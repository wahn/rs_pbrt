use core::str;
use std::net::TcpStream;
use std::sync::atomic::AtomicBool;
use std::sync::{Arc, Mutex, RwLock};
use std::thread::{self, ThreadId};
use std::time;

use atomic::Ordering;
use crossbeam_utils::thread::Scope;
use tev_client::{PacketCreateImage, PacketUpdateImage, TevClient, TevError};

use crate::core::geometry::{Bounds2i, Point2i};
use crate::core::pbrt::Float;

const TILE_SIZE: i32 = 128;

type Pixels<'a, T> = Arc<&'a RwLock<Vec<T>>>;
type Function<T> = fn(Bounds2i, Pixels<T>, usize, &mut Vec<Vec<Float>>);

struct DisplayItem<'a, T: Send + Sync> {
    title: String,
    resolution: Point2i,
    get_tile_values: Function<T>,
    // channel_buffers: Vec<ImageChannelBuffer>,
    channel_names: Vec<String>,
    vec: Pixels<'a, T>,
    opened_image: bool,
}

impl<'a, T: Send + Sync> DisplayItem<'a, T> {
    pub fn new(
        base_title: &str,
        resolution: Point2i,
        channel_names: Vec<String>,
        vec: Pixels<'a, T>,
        get_tile_values: Function<T>,
    ) -> DisplayItem<'a, T> {
        let title = format!("{} {:?}", base_title, get_thread_id());

        DisplayItem {
            title,
            resolution,
            get_tile_values,
            channel_names,
            vec,
            opened_image: false,
        }
    }

    pub fn display_with_tev_client(&mut self, client: &mut TevClient) -> bool {
        // Open image if not opened already
        if !self.opened_image {
            if self.send_create_image_packet(client).is_ok() {
                self.opened_image = true;
            } else {
                return false;
            }
        }

        // Create image buffer
        let size = self.channel_names.len();
        let inner_size = (TILE_SIZE * TILE_SIZE) as usize;
        let mut display_values: Vec<Vec<Float>> = Vec::with_capacity(size);
        for _ in 0..size {
            display_values.push(Vec::with_capacity(inner_size));
        }

        // Bounds for Tile
        let bounds = Bounds2i {
            p_min: Point2i { x: 0, y: 0 },
            p_max: self.resolution,
        };

        // Retrieve image values for Tile with bounds
        (self.get_tile_values)(
            bounds,
            self.vec.clone(),
            self.resolution.x as usize,
            &mut display_values,
        );

        // debug_assert!(!display_values.iter().all(|x| x.iter().all(|y| *y == 0.0)),
        //               "All display values are zero");

        let channel0 = display_values[0].iter();
        let channel1 = display_values[1].iter();
        let channel2 = display_values[2].iter();
        let mut data: Vec<f32> = vec![];
        data.push(1.0);
        for ((x, y), z) in channel0.zip(channel1).zip(channel2) {
            data.push(*x);
            data.push(*y);
            data.push(*z);
        }

        let packet = PacketUpdateImage {
            image_name: &self.title,
            grab_focus: false,
            channel_names: &["R", "G", "B"],
            channel_offsets: &[1, 2, 3],
            channel_strides: &[3, 3, 3],
            x: 0,
            y: 0,
            width: self.resolution.x as u32,
            height: self.resolution.y as u32,
            data: &data,
        };

        client.send(packet).expect("TODO: panic message");

        true
    }

    fn send_create_image_packet(&self, client: &mut TevClient) -> std::io::Result<()> {
        client.send(PacketCreateImage {
            image_name: &self.title,
            grab_focus: false,
            width: self.resolution.x as u32,
            height: self.resolution.y as u32,
            channel_names: &["R", "G", "B"],
        })
    }
}

/*impl ImageChannelBuffer {
   fn send_if_changed(&mut self, ipc_channel: &mut TcpStream, tile_index: usize) -> bool {
        let hash = murmur_hash64a(&self.buffer[self.channel_values_offset..], 0);
        if let Some(&tile_hash) = self.tile_hashes.get(tile_index) {
            if tile_hash == hash {
                return true;
            }
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
}*/

fn get_thread_id() -> ThreadId {
    thread::current().id()
}

pub struct Preview<'a, T: Send + Sync> {
    pub exit_thread: Arc<AtomicBool>,
    dynamic_items: Arc<Mutex<Vec<DisplayItem<'a, T>>>>,
    dynamic_channel: TcpStream,
}

impl<'a, T: Send + Sync> Preview<'a, T> {
    pub fn connect_to_display_server(host: &str) -> Result<Preview<'a, T>, TevError> {
        let exit_thread = Arc::new(AtomicBool::new(false));
        let dynamic_channel = TcpStream::connect(host)?;
        let dynamic_items: Arc<Mutex<Vec<DisplayItem<T>>>> = Arc::new(Mutex::new(Vec::new()));

        Ok(Preview {
            exit_thread,
            dynamic_items,
            dynamic_channel,
        })
    }

    pub fn disconnect_from_display_server(&mut self) {
        dbg!("Disconnecting from Tev");
        self.exit_thread.store(true, Ordering::Relaxed);
    }

    pub fn display_dynamic(
        self,
        title: &str,
        resolution: Point2i,
        channel_names: Vec<String>,
        scope: &Scope<'a>,
        arc: Arc<&'a RwLock<Vec<T>>>,
        get_tile_values: Function<T>,
    ) {
        let cloned_exit_thread = self.exit_thread.clone();
        let cloned_dynamic_items = self.dynamic_items.clone();
        let mut cloned_channel = self.dynamic_channel.try_clone().unwrap();
        scope.spawn(move |_| {
            update_dynamic_items(
                cloned_exit_thread,
                &mut cloned_channel,
                cloned_dynamic_items,
            );
        });

        let cloned_items = self.dynamic_items.clone();
        let mut display_items = cloned_items.lock().unwrap();
        display_items.push(DisplayItem::new(
            title,
            resolution,
            channel_names,
            arc,
            get_tile_values,
        ));
    }

    #[allow(unused)]
    fn display_static(
        &mut self,
        title: &str,
        resolution: Point2i,
        vec: Pixels<T>,
        channel_names: Vec<String>,
        get_tile_values: Function<T>,
    ) {
        let mut item = DisplayItem::new(title, resolution, channel_names, vec, get_tile_values);
        let mut client = TevClient::wrap(self.dynamic_channel.try_clone().unwrap());
        if !item.display_with_tev_client(&mut client) {
            println!("Unable to display static content {}", title);
        }
    }
}

fn update_dynamic_items<T: Send + Sync>(
    exit_thread: Arc<AtomicBool>,
    channel: &mut TcpStream,
    items: Arc<Mutex<Vec<DisplayItem<T>>>>,
) {
    let mut client = TevClient::wrap(channel.try_clone().expect("Could not clone TCP Channel"));
    while !exit_thread.load(Ordering::Relaxed) {
        thread::sleep(time::Duration::from_millis(250));

        let mut items = items.lock().expect("Could not lock");
        for item in items.iter_mut() {
            item.display_with_tev_client(&mut client);
        }
    }

    let mut items = items.lock().expect("Could not lock");
    for item in items.iter_mut() {
        item.display_with_tev_client(&mut client);
    }

    items.clear();
}

#[cfg(test)]
mod test {
    use std::sync::atomic::Ordering::Relaxed;
    use std::sync::{Arc, Mutex, RwLock};
    use std::thread;
    use std::time;
    use std::time::Duration;

    use crossbeam_utils::thread::Scope;

    use crate::core::display::Preview;
    use crate::core::film::Pixel;
    use crate::core::geometry::{Bounds2i, Point2i};
    use crate::core::pbrt::Float;

    #[ignore]
    #[test]
    /// Manual test for tev remote
    fn display_remote() {
        let address = "127.0.0.1:14158";

        let display = Preview::connect_to_display_server(address);
        // Do not fail the test if tev is not running
        if display.is_err() {
            return;
        }
        let display = display.unwrap();
        let exit_thread = display.exit_thread.clone();
        let resolution = Point2i { x: 200, y: 200 };

        let mut image: Vec<Pixel> = Vec::with_capacity(resolution.x as usize);
        for x in 0..resolution.x {
            for y in 0..resolution.y {
                let mut pixel = Pixel::default();
                pixel.xyz = [
                    x as Float / resolution.x as Float,
                    y as Float / resolution.y as Float,
                    0.0,
                ];
                image.push(pixel);
            }
        }

        let data = &RwLock::new(image);
        let arc = Arc::new(data);

        let get_values = move |b: Bounds2i,
                               arc: Arc<&RwLock<Vec<Pixel>>>,
                               width: usize,
                               values: &mut Vec<Vec<Float>>| {
            for col in b.p_min.y as usize..b.p_max.y as usize {
                for row in b.p_min.x as usize..b.p_max.x as usize {
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

        // let display = Arc::new(Mutex::new(display));
        crossbeam::scope(|scope| {
            display.display_dynamic(
                "rs_pbrt",
                resolution,
                vec!["R".to_string(), "G".to_string(), "B".to_string()],
                scope,
                arc.clone(),
                get_values,
            );

            thread::sleep(time::Duration::from_millis(1000));
            for cols in 0..resolution.x as usize {
                for rows in 0..resolution.y as usize {
                    let mut arc = arc.write().unwrap();
                    arc[cols * resolution.x as usize + rows] = Pixel::default();
                }
            }
            thread::sleep(time::Duration::from_millis(1000));
            exit_thread.store(true, Relaxed);
        })
        .unwrap();
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
        })
        .unwrap();

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
