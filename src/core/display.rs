use core::str;
use std::cmp::min;
use std::io::Write;
use std::mem::size_of;
use std::net::TcpStream;
use std::sync::{Arc, Mutex, RwLock};
use std::sync::atomic::AtomicBool;
use std::thread::{self, ThreadId};
use std::{time, io};

use atomic::Ordering;
use murmurhash64::murmur_hash64a;

use crate::core::geometry::{Bounds2i, Point2i};
use crate::core::pbrt::Float;
use crossbeam_utils::thread::Scope;
use crate::core::film::Pixel;

const TILE_SIZE: i32 = 128;

type Pixels<'a> = Arc<&'a RwLock<Vec<Pixel>>>;
type Function = fn(Bounds2i, Pixels, usize, &mut Vec<Vec<Float>>);

struct DisplayDirective;

#[allow(unused)]
impl DisplayDirective {
    pub const OPEN_IMAGE: u8 = 0;
    pub const RELOAD_IMAGE: u8 = 1;
    pub const CLOSE_IMAGE: u8 = 2;
    pub const UPDATE_IMAGE: u8 = 3;
    pub const CREATE_IMAGE: u8 = 4;
}

struct DisplayItem<'a> {
    title: String,
    resolution: Point2i,
    get_tile_values: Function,
    channel_buffers: Vec<ImageChannelBuffer>,
    vec: Pixels<'a>,
    opened_image: bool,
}

impl<'a> DisplayItem<'a> {
    pub fn new(
        base_title: &str,
        resolution: Point2i,
        channel_names: Vec<String>,
        vec: Pixels<'a>,
        get_tile_values: Function,
    ) -> DisplayItem<'a> {
        let title = format!("{} {:?}", base_title, get_thread_id());

        let mut channel_buffers = Vec::new();
        for channel_name in channel_names {
            channel_buffers.push(ImageChannelBuffer::new(
                channel_name,
                title.clone(),
            ));
        }

        DisplayItem {
            title,
            resolution,
            get_tile_values,
            channel_buffers,
            vec,
            opened_image: false,
        }
    }

    pub fn display(&mut self, ipc_channel: &mut TcpStream) -> bool {
        // Open image in tev if not done already
        if !self.opened_image {
            let opened = self.send_open_image(ipc_channel);
            if let Err(err) = opened {
                println!("Could not send open image packet: {}", err);
                return false;
            }
            self.opened_image = true;
        }

        // Create image buffer
        let size = self.channel_buffers.len();
        let inner_size = (TILE_SIZE * TILE_SIZE) as usize;
        let mut display_values: Vec<Vec<Float>> = Vec::with_capacity(size);
        for _ in 0..size { display_values.push(Vec::with_capacity(inner_size)); }

        let mut tile_index = 0;
        let mut y = 0;
        while y < self.resolution.y {
            let mut x = 0;
            while x < self.resolution.x {
                let height = min(y + TILE_SIZE, self.resolution.y) - y;
                let width = min(x + TILE_SIZE, self.resolution.x) - x;

                for channel_buffer in &mut self.channel_buffers {
                    channel_buffer.buffer.resize(channel_buffer.channel_values_offset, 0);
                    channel_buffer.set_tile_bounds(x, y, width, height);
                }

                let b = Bounds2i {
                    p_min: Point2i { x, y },
                    p_max: Point2i {
                        x: x + width,
                        y: y + height,
                    },
                };

                // empty out display old values
                for values in display_values.iter_mut() { values.clear(); }
                (self.get_tile_values)(b, self.vec.clone(), self.resolution.x as usize, &mut display_values);
                // Send the RGB buffers if they differ from the last sent version
                for (channel_buffer, values) in self.channel_buffers.iter_mut().zip(&display_values) {
                    // Insert image values into channel buffer
                    let values: Vec<u8> = values.iter().flat_map(|x| x.to_le_bytes()).collect();
                    channel_buffer.buffer.extend_from_slice(&values);
                }

                for channel_buffer in self.channel_buffers.iter_mut() {
                    let sent = channel_buffer.send_if_changed(ipc_channel, tile_index);
                    if !sent {
                        // self.opened_image = false;
                        return false;
                    }
                }

                x += TILE_SIZE;
                tile_index += 1;
            }

            y += TILE_SIZE;
        }

        true
    }

    fn send_open_image(&self, ipc_channel: &mut TcpStream) -> std::io::Result<usize> {
        // Create "open the image" message buffer
        let mut buffer = Vec::with_capacity(1024);
        create_open_packet(&self.title, &self.resolution, &mut buffer);

        // Send off the buffer
        ipc_channel.write(&buffer)
    }
}

struct ImageChannelBuffer {
    buffer: Vec<u8>,
    channel_values_offset: usize,
    tile_hashes: Vec<u64>,
    // Should probably be a map?
}

impl ImageChannelBuffer {
    fn new(channel_name: String, title: String) -> Self {
        let buffer_size =
            TILE_SIZE as usize * TILE_SIZE as usize * size_of::<Float>() + title.capacity() + 32;
        let mut buffer = Vec::with_capacity(buffer_size);

        buffer.extend_from_slice(&0_i32.to_le_bytes()); // reserve space for message length. Maybe we should check somewhere, that the message length will fit in 32 bits
        buffer.extend_from_slice(&DisplayDirective::UPDATE_IMAGE.to_le_bytes());
        buffer.push(0); // grab focus
        buffer.extend_from_slice(title.as_bytes());
        buffer.push(0); // Null terminate string
        buffer.extend_from_slice(channel_name.as_bytes());
        buffer.push(0); // Null terminate string

        // The original implementation says we have a problem with the offset now not being float-aligned
        // As far as I can tell this should be u8 aligned
        // Also this is apparently not a problem on x86...
        let channel_values_offset = buffer.len();

        let tile_hashes: Vec<u64> = Vec::new();

        ImageChannelBuffer {
            buffer,
            channel_values_offset,
            tile_hashes,
        }
    }

    fn set_tile_bounds(&mut self, x: i32, y: i32, width: i32, height: i32) {
        self.buffer.extend_from_slice(&x.to_le_bytes());
        self.buffer.extend_from_slice(&y.to_le_bytes());
        self.buffer.extend_from_slice(&width.to_le_bytes());
        self.buffer.extend_from_slice(&height.to_le_bytes());
    }

    fn send_if_changed(&mut self, ipc_channel: &mut TcpStream, tile_index: usize) -> bool {
        let hash = murmur_hash64a(&self.buffer[self.channel_values_offset..], 0);
        if let Some(&tile_hash) = self.tile_hashes.get(tile_index) {
            if tile_hash == hash { return true; }
        }

        let message_length = self.buffer.len();
        let message_length_bytes = (message_length as i32).to_ne_bytes();
        self.buffer.splice(..4, message_length_bytes);

        let sent = ipc_channel.write(&self.buffer);
        if let Err(err) = sent {
            dbg!(err);
            return false;
        }

        self.tile_hashes.insert(tile_index, hash);
        true
    }
}

fn get_thread_id() -> ThreadId {
    thread::current().id()
}

fn create_open_packet(title: &str, resolution: &Point2i, buf: &mut Vec<u8>) {
    buf.extend_from_slice(&0_i32.to_le_bytes()); // Make space for message length
    buf.push(DisplayDirective::CREATE_IMAGE); // Create Image
    buf.push(1); // Grab focus
    buf.extend_from_slice(title.as_bytes());
    buf.push(0); // Null terminate string
    buf.extend_from_slice(&resolution.x.to_le_bytes());
    buf.extend_from_slice(&resolution.y.to_le_bytes());

    let channels = ["R", "G", "B"];
    let size = channels.len() as i32;
    buf.extend_from_slice(&size.to_le_bytes());
    for channel in channels {
        buf.extend_from_slice(channel.as_bytes());
        buf.push(0); // Null terminate string
    }

    let message_length = buf.len();
    let message_length_bytes = (message_length as i32).to_le_bytes();
    buf.splice(..4, message_length_bytes);
}

pub struct Preview<'a> {
    exit_thread: Arc<AtomicBool>,
    dynamic_items: Arc<Mutex<Vec<DisplayItem<'a>>>>,
    // update_thread: Option<JoinHandle<()>>,
    dynamic_channel: TcpStream,
}

impl<'a> Preview<'a> {
    pub fn connect_to_display_server(host: &str) -> io::Result<Preview> {
        let exit_thread = Arc::new(AtomicBool::new(false));
        let dynamic_channel = TcpStream::connect(host)?;
        let dynamic_items: Arc<Mutex<Vec<DisplayItem>>> = Arc::new(Mutex::new(Vec::new()));

        Ok(Preview {
            exit_thread,
            dynamic_items,
            dynamic_channel,
        })
    }

    pub fn disconnect_from_display_server(&mut self) {
        self.exit_thread.store(true, Ordering::Relaxed);
    }

    pub fn display_dynamic(&mut self, title: &str, resolution: Point2i, channel_names: Vec<String>,
                           scope: &Scope<'a>, arc: Arc<&'a RwLock<Vec<Pixel>>>, get_tile_values: Function) {
        dbg!("display_dynamic");

        let cloned_exit_thread = self.exit_thread.clone();
        let cloned_dynamic_items = self.dynamic_items.clone();
        let mut cloned_channel = self.dynamic_channel.try_clone().unwrap();
        // self.update_thread = Some(
        scope.spawn(move |_| {
            update_dynamic_items(cloned_exit_thread, &mut cloned_channel, cloned_dynamic_items);
        });
        // );
        let mut display_items = self.dynamic_items.lock().unwrap();
        display_items.push(DisplayItem::new(title, resolution, channel_names, arc, get_tile_values))
    }

    #[allow(unused)]
    fn display_static(&mut self, title: &str, resolution: Point2i, vec: Pixels, channel_names: Vec<String>,
                      get_tile_values: Function) {
        let mut item = DisplayItem::new(title, resolution, channel_names, vec, get_tile_values);
        if !item.display( &mut self.dynamic_channel) {
            println!("Unable to display static content {}", title);
        }
    }
}

fn update_dynamic_items(exit_thread: Arc<AtomicBool>, channel: &mut TcpStream, items: Arc<Mutex<Vec<DisplayItem>>>) {
    while !exit_thread.load(Ordering::Relaxed) {
        thread::sleep(time::Duration::from_millis(250));

        let mut items = items.lock().unwrap();
        for item in items.iter_mut() {
            item.display( channel);
        }
    }

    let mut items = items.lock().unwrap();
    for item in items.iter_mut() {
        item.display(channel);
    }

    items.clear();
}


#[cfg(test)]
mod test {
    use std::thread;
    use std::time;

    use crate::core::display::Preview;
    use crate::core::geometry::{Bounds2i, Point2i};
    use crate::core::pbrt::Float;
    use std::sync::{Arc, RwLock, Mutex};
    use std::time::Duration;
    use crossbeam_utils::thread::Scope;
    use crate::core::film::Pixel;


    #[test]
    /// Manual test for tev remote
    fn display_remote() {
        let address = "127.0.0.1:14158";

        let display = Preview::connect_to_display_server(address);
        assert!(display.is_ok());
        let mut display = display.unwrap();
        let resolution = Point2i { x: 200, y: 200 };

        let mut image: Vec<Pixel> = Vec::with_capacity(resolution.x as usize);
        for x in 0..resolution.x {
            for y in 0..resolution.y {
                let color = (x * y) as Float / (resolution.x * resolution.y - 1) as Float;
                let mut pixel = Pixel::default();
                pixel.xyz = [color; 3];
                image.push(pixel);
            }
        }

        let data = &RwLock::new(image);
        let arc = Arc::new(data);

        let get_values = move |b: Bounds2i, arc: Arc<&RwLock<Vec<Pixel>>>, width: usize, values: &mut Vec<Vec<Float>>| {
            for col in b.p_min.y as usize.. b.p_max.y as usize {
                for row in b.p_min.x as usize .. b.p_max.x as usize {
                    let v = {
                        let clone = arc.read().unwrap();
                        clone[col * width + row].xyz
                    };

                    for i in 0..3 {
                        values[i].push(v[i]);
                    }
                }
            }
        };

        crossbeam::scope(|scope| {
            display.display_dynamic("Test", resolution,
                                    vec!["R".to_string(), "G".to_string(), "B".to_string()], scope, arc.clone(), get_values);

            thread::sleep(time::Duration::from_millis(1000));
            for cols in 0..resolution.x as usize {
                for rows in 0..resolution.y as usize {
                    let mut arc = arc.write().unwrap();
                    arc[cols * resolution.x as usize + rows] = Pixel::default();
                }
            }
            thread::sleep(time::Duration::from_millis(1000));
            display.disconnect_from_display_server();
        }).unwrap();
    }

    #[test]
    fn mutate_data_while_sharing() {
        let num = (0, 0);

        let arc = Arc::new(Mutex::new(num));
        let clone = arc.clone();
        let clone2 = arc.clone();

        crossbeam::scope(|scope| {
            mutate(scope, clone);
            for _ in 0..100 {
                {
                    let mut arc = arc.lock().unwrap();
                    arc.0 += 1;
                }
                thread::sleep(Duration::from_millis(3));
            }
        }).unwrap();

        let num = clone2.lock().unwrap();
        println!("{:?}", num);
        println!("Done.");
    }

    fn mutate(scope: &Scope, clone: Arc<Mutex<(i32, i32)>>) {
        scope.spawn(move |_| {
            for _ in 0..10 {
                {
                    let mut clone = clone.lock().unwrap();
                    println!("{:?}", clone.0);
                    clone.1 += 1;
                }
                thread::sleep(Duration::from_millis(2));
            }
        });
    }
}
