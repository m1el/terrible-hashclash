use rand::{Rng, SeedableRng};

// use std::hash::{Hash, Hasher};
// use std::collections::hash_map::{DefaultHasher, HashMap, Entry};
use std::collections::hash_map::{HashMap, Entry};
use std::sync::atomic::{AtomicUsize, AtomicU64, AtomicU32, Ordering};
use std::sync::{Arc, Mutex};
use std::convert::TryInto;
use std::thread;

mod xoroshiro;
mod md5;
use xoroshiro::Xoroshiro128Plus;
use md5::compress;


fn next_point(state: [u32; 3], mask: u32) -> [u32; 3] {
    let mut ihv = if state[0] < state[1] { IHV1 } else { IHV2 };
    let mut data = [0_u32; 16];
    /*
    let mut hex_state = [0_u8; 4 * 3 * 2];
    fn hex(x: u32) -> u8 {
        assert!(x <= 15, "YOU LIED");
        (x as u8) + if x <= 9 { b'0' } else { b'A' - 10 }
    }
    for (si, sb) in state.iter().cloned().enumerate() {
        for ii in 0..4 {
            let byte = sb >> ii * 8;
            let lo = hex(byte & 0xf);
            let hi = hex((byte & 0xf0) >> 4);
            hex_state[si * 8 + ii * 2] = lo;
            hex_state[si * 8 + ii * 2 + 1] = hi;
        }
    }
    let hex_state = unsafe { core::mem::transmute::<_, [u32; 6]>(hex_state) };

    data[10..].copy_from_slice(&hex_state);
    */
    // println!("ihv={:x?} data={:x?}", ihv, data);
    data[13..].copy_from_slice(&state);
    compress(&mut ihv, &data);
    let mut rv: [u32; 3] = ihv[..3].try_into().unwrap();
    rv[2] &= mask;
    return rv;
}

fn find_collision(mut a: [u32; 3], len_a: u64, mut b: [u32; 3], len_b: u64, mask: u32) -> ([u32; 3], [u32; 3]) {
    let len = if len_a > len_b {
        for ii in 0..(len_a - len_b) {
            a = next_point(a, mask);
        }
        len_b
    } else {
        for ii in 0..(len_b - len_a) {
            b = next_point(b, mask);
        }
        len_a
    };
    loop {
        let pa = a;
        let pb = b;
        a = next_point(a, mask);
        b = next_point(b, mask);
        if a == b {
            return (pa, pb);
        }
    }
}

/*
fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}
*/

static SYNC_POINT: AtomicUsize = AtomicUsize::new(0);

static TRAILS: AtomicU64 = AtomicU64::new(0);
static COLLISIONS: AtomicU64 = AtomicU64::new(0);
static ROBIN_HOOD: AtomicU64 = AtomicU64::new(0);
static LOOPS: AtomicU64 = AtomicU64::new(0);
static HASH_COUNT: AtomicU64 = AtomicU64::new(0);
static CONTENTIONS: AtomicU64 = AtomicU64::new(0);

static HASH_MASK: AtomicU32 = AtomicU32::new(0);
static TRAIL_MASK: AtomicU32 = AtomicU32::new(0);
const IHV1: [u32; 4] = [ 0x266e2670, 0x9b8a1b87, 0x923fd523, 0x8c4fcf12 ];
const IHV2: [u32; 4] = [ 0xa0da787b, 0xb3cb406c, 0xfe644118, 0xd7c59003 ];

fn main() {
    let trails = HashMap::<[u32; 3], Vec<([u32; 3], u64)>>::new();
    let trails = Arc::new(Mutex::new(trails));
    let thread_count = num_cpus::get() as u64;
    SYNC_POINT.store(1, Ordering::SeqCst);

    let mut workers = Vec::new();
    for worker in 0..thread_count {
        let trails = trails.clone();
        let mut worker_rng = Xoroshiro128Plus::from_seed([0x1337133713371337 + worker, 0x1337133713371337 + worker * 2]);
        workers.push(thread::spawn(move || loop {
            // wait to start work
            while SYNC_POINT.load(Ordering::SeqCst) != 0 {
                std::hint::spin_loop();
            }

            let hash_mask = HASH_MASK.load(Ordering::SeqCst);
            let trail_mask = TRAIL_MASK.load(Ordering::SeqCst);

            'outer: loop {

                let start = [worker_rng.next_u32(), worker_rng.next_u32(), worker_rng.next_u32()];
                // let mut value = start;
                let mut trail_steps = 0;
                let mut value = start;
                // Generate a trail.
                // The loop will return a distinguishing point.
                let hash = loop {
                    if SYNC_POINT.load(Ordering::Relaxed) != 0 { break 'outer }
                    trail_steps += 1;

                    let hash = next_point(value, hash_mask);
                    // How many hashes to put in the same trail. This is 1/theta.
                    if hash[0] & trail_mask == 0 {
                        TRAILS.fetch_add(1, Ordering::Relaxed);
                        HASH_COUNT.fetch_add(trail_steps, Ordering::Relaxed);
                        break hash;
                    }
                    // give up on this trail
                    if trail_steps > 20 * (trail_mask as u64) {
                        LOOPS.fetch_add(1, Ordering::Relaxed);
                        HASH_COUNT.fetch_add(trail_steps, Ordering::Relaxed);
                        continue 'outer;
                    }
                    value = hash;
                };

                let mut trails = if let Ok(lock) = trails.try_lock() {
                    lock
                } else {
                    CONTENTIONS.fetch_add(1, Ordering::Relaxed);
                    trails.lock().unwrap()
                };

                let new_trail = (start, trail_steps);
                let mut collided = false;
                match trails.entry(hash) {
                    Entry::Occupied(mut entry) => {
                        let trails = entry.get_mut();
                        for (old_start, old_len) in trails.iter() {
                            if hash == *old_start {
                                ROBIN_HOOD.fetch_add(1, Ordering::Relaxed);
                                continue;
                            }
                            let (ca, cb) = find_collision(*old_start, *old_len, start, trail_steps, hash_mask);
                            let ta = ca[0] < ca[1];
                            let tb = cb[0] < cb[1];
                            if ta == tb {
                                println!("found self-collision :(");
                            } else {
                                println!("found a good collision! ta={} tb={} ca={:x?} cb={:x?}", ta, tb, ca, cb);
                                const MAX_COLLISIONS: u64 = 1;
                                if COLLISIONS.fetch_add(1, Ordering::Relaxed) >= MAX_COLLISIONS - 1 {
                                    let _ = SYNC_POINT.compare_exchange(0, 1, Ordering::SeqCst, Ordering::SeqCst);
                                }
                            }
                            break;
                        }
                        trails.push(new_trail);
                        collided = true;
                    }
                    Entry::Vacant(entry) => {
                        entry.insert(vec![new_trail]);
                    }
                };
                if collided {
                    trails.clear();
                }
            }

            // Notify that current thread has finished its work
            SYNC_POINT.fetch_add(1, Ordering::SeqCst);
        }));
    }


    println!("hash_bits,trail_bits,hashes_per_coll,time_to_coll,loops_per_coll,contentions_per_sec,robin_hood_per_coll");

    let memory_bits = 44;
    // let hash_bits = 64;
    for hash_bits in (64_usize..=96).step_by(4) {
    let min_trail_bits = hash_bits.saturating_sub(memory_bits).max(6);
    let max_trail_bits = (hash_bits / 2).min(32);
    // let trail_bits = 12;
    for trail_bits in min_trail_bits..=max_trail_bits {

    let hash_mask = if hash_bits == 64 { 0 } else { (!0_u32) >> (96 - hash_bits) };
    let trail_mask = (1_u32 << trail_bits) - 1;
    println!("hash_mask={:016x}, trail_mask={:016x}", hash_mask, trail_mask);
    HASH_MASK.store(hash_mask, Ordering::SeqCst);
    TRAIL_MASK.store(trail_mask, Ordering::SeqCst);

    let start = std::time::Instant::now();
    SYNC_POINT.store(0, Ordering::SeqCst);

    // wait for all workers to finish their work
    while SYNC_POINT.load(Ordering::SeqCst) < workers.len() + 1 {
        std::hint::spin_loop();
    }
    let stat_time = start.elapsed().as_millis();

    let loops = LOOPS.load(Ordering::Relaxed);
    let contentions = CONTENTIONS.load(Ordering::Relaxed);
    let collisions = COLLISIONS.load(Ordering::Relaxed);
    let hashes = HASH_COUNT.load(Ordering::Relaxed);
    let robin_hood = ROBIN_HOOD.load(Ordering::Relaxed);
    println!("{:2},{:2},{:4.0},{:2.5},{:4},{:4},{:4}",
             hash_bits, trail_bits,
             (hashes as f64) / (collisions as f64),
             (stat_time as f64 / 1000.0) / (collisions as f64),
             (loops as f64) / (collisions as f64),
             (contentions as f64) / (stat_time as f64 / 1000.0),
             (robin_hood as f64) / (collisions as f64));

    trails.lock().unwrap().clear();
    LOOPS.store(0, Ordering::SeqCst);
    CONTENTIONS.store(0, Ordering::SeqCst);
    COLLISIONS.store(0, Ordering::SeqCst);
    HASH_COUNT.store(0, Ordering::SeqCst);
    ROBIN_HOOD.store(0, Ordering::SeqCst);
    TRAILS.store(0, Ordering::SeqCst);
    }
    }
}
