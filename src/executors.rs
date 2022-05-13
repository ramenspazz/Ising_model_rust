use executors::*;
use executors::crossbeam_workstealing_pool;
use std::sync::{Arc, RwLock};
use std::sync::mpsc::channel;

struct lattice {
    test_vec: Arc<RwLock<Vec<i32>>>,
}

impl lattice {
    fn new() -> Self {
        Self {
            test_vec: Arc::new(RwLock::new(vec![])),
        }
    }

    fn push(&mut self, item: i32) {
        match self.test_vec.write().as_mut() {
            Ok(value) => value.push(item),
            Err(_) => return,
        }
    }
}

fn thread_func(tx: std::sync::mpsc::Sender<i32>, vec_obj: &Arc<RwLock<Vec<i32>>>) {
    // get a ref to the test_vec and then sum it
    let mut psum = 0;
    match vec_obj.read().as_deref() {
        Ok(value) => {
            for item in value {
                psum += item;
            }
        },
        Err(_) => return, 
    }
    tx.send(psum).unwrap();
}

fn main() {
    let n_workers = 1;
    let n_jobs = 1;
    let pool = crossbeam_workstealing_pool::small_pool(n_workers);
    let mut test_lat = lattice::new();

    for i in 1..5 {
        test_lat.push(i);
    }

    let (tx, rx) = channel();
    // let (finished_tx, finished_rx) = channel();

    for _ in 0..n_jobs {
        let shared_data = test_lat.test_vec.clone();
        let tx = tx.clone();
        pool.execute(move || thread_func(tx, &shared_data));
    }

    println!("{}", rx.iter().take(n_jobs).fold(0, |a, b| a + b));
    println!("Ok");
}