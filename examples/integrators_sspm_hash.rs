use pbrt::core::geometry::Point3i;

fn hash(p: &Point3i, hash_size: i32) -> usize {
    (((p.x * 73856093) ^ (p.y * 19349663) ^ (p.z * 83492791)) as u32 % hash_size as u32) as usize
}

fn main() {
    let photon_grid_index: Point3i = Point3i {
        x: 225,
        y: 0,
        z: 267,
    };
    let hash_size: usize = 700000;
    let h: usize = hash(&photon_grid_index, hash_size as i32);
    assert!(h < hash_size, "hash({:?}, {:?})", photon_grid_index, hash_size);
    println!("hash({:?}, {:?}) = {}", photon_grid_index, hash_size, h);
}
