use rand::{Rng, SeedableRng};

use std::hash::{Hash, Hasher};
use std::collections::hash_map::{DefaultHasher, HashMap, Entry};
use std::sync::atomic::{AtomicUsize, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;

mod xoroshiro;
use xoroshiro::Xoroshiro128Plus;

fn calculate_hash<T: Hash>(t: &T) -> u64 {
    let mut s = DefaultHasher::new();
    t.hash(&mut s);
    s.finish()
}

static SYNC_POINT: AtomicUsize = AtomicUsize::new(0);

static TRAILS: AtomicU64 = AtomicU64::new(0);
static COLLISIONS: AtomicU64 = AtomicU64::new(0);
static ROBIN_HOOD: AtomicU64 = AtomicU64::new(0);
static LOOPS: AtomicU64 = AtomicU64::new(0);
static HASH_COUNT: AtomicU64 = AtomicU64::new(0);
static CONTENTIONS: AtomicU64 = AtomicU64::new(0);

static HASH_MASK: AtomicU64 = AtomicU64::new(0);
static TRAIL_MASK: AtomicU64 = AtomicU64::new(0);

fn main() {

    let trails = HashMap::<u64, Vec<(u64, u64)>>::new();
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
                let start = worker_rng.next_u64();
                // let mut value = start;
                let mut trail_steps = 0;
                let mut value = start;
                // Generate a trail.
                // The loop will return a distinguishing point.
                let hash = loop {
                    if SYNC_POINT.load(Ordering::Relaxed) != 0 { break 'outer }
                    trail_steps += 1;

                    // Limit hash size for faster collision testing.
                    // How many hashes to put in the same trail. This is 1/theta.
                    let hash_full = calculate_hash(&value);
                    let hash = hash_full & hash_mask;
                    if hash & trail_mask == 0 {
                        TRAILS.fetch_add(1, Ordering::Relaxed);
                        HASH_COUNT.fetch_add(trail_steps, Ordering::Relaxed);
                        break hash;
                    }
                    // give up on this trail
                    if trail_steps > 20 * trail_mask {
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
                        for (old_start, _) in trails.iter() {
                            if hash == *old_start {
                                ROBIN_HOOD.fetch_add(1, Ordering::Relaxed);
                                continue;
                            }
                            const MAX_COLLISIONS: u64 = 1000;
                            if COLLISIONS.fetch_add(1, Ordering::Relaxed) >= MAX_COLLISIONS - 1 {
                                SYNC_POINT.compare_exchange(0, 1, Ordering::SeqCst, Ordering::SeqCst);
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
    for hash_bits in (48_usize..=64).step_by(4) {
    let min_trail_bits = hash_bits.saturating_sub(memory_bits).max(6);
    let max_trail_bits = (hash_bits / 2).min(32);
    // let trail_bits = 12;
    for trail_bits in min_trail_bits..=max_trail_bits {

    let hash_mask = (!0) >> (64 - hash_bits);
    let trail_mask = (1 << trail_bits) - 1;
    // println!("hash_mask={:016x}, trail_mask={:016x}", hash_mask, trail_mask);
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
