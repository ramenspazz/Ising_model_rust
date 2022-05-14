use std::{thread as stdth, time};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use core::ops::Range;
use indicatif::ProgressBar;
use std::f64::consts;
use crossbeam_channel::Sender;
use rand::thread_rng;
use rand::Rng;
use std::sync::{RwLock, Arc};
use crossbeam_channel::Receiver;
use executors::parker::SmallThreadData;
use ndarray::prelude::*;
use executors::Executor;
use crossbeam_utils::thread;
use crossbeam_channel::{bounded, unbounded};
use executors::crossbeam_workstealing_pool;
use executors::parker::StaticParker;
use crate::atomistic_simulation::crossbeam_workstealing_pool::ThreadPool;
use crate::dividend_remainder;
use crate::lattice_structure::Lattice;
use crate::lat_node::SpinNode;

fn indexmod(index: usize, modnum: usize, x_size: usize, y_size: usize, sym_type: SymmetryType) -> Option<usize> {
    match sym_type {
        SymmetryType::C3V => {
            // find a way to refer to nodes by their cordinate
        },
        SymmetryType::C4V => {
            if modnum == 0 {
                if index + 1 < x_size*y_size {
                    return Some(index+1);
                } else { return None; }
            } else
            if modnum == 1 {
                if index + x_size < x_size*y_size {
                    return Some(index + x_size);
                } else { return None; }
            } else
            if modnum == 2 {
                if index > x_size && index < x_size*y_size {
                    return Some(index - x_size);
                } else { return None; }
            } else
            if modnum == 3 {
                if index > 0 && index < x_size*y_size {
                    return Some(index - 1);
                } else { return None; }
            }
        },
        SymmetryType::C6V => {
            let i = (index / y_size) as usize;
            let j = index - i * y_size;
            if modnum == 0 {
                if j > 0 {
                    // if the y index is greater than 0, we can add the -y node
                    return Some(index - 1);
                }
                
            } else
            if modnum == 1 {
                if (j + 2) < y_size {
                    // if the y index is 2 less than the  
                    return Some(index + 1);
                }
            } else
            if modnum == 2 {
                if i % 2 == 0 {
                    // when we are on an even increment the +x direction
                    if (j % 2 == 0) && (index >= y_size) {
                       // if the y index is even, then its neighbor is in the -x direction
                        return Some(index - y_size);
                    } else
                    if (j % 2 == 1) && ((index + y_size + 1) <= x_size*y_size) {
                        // if the y index is odd, then its neighbor is in the -x direction
                        return Some(index + y_size);
                    }
                } else {
                    // when we are on an odd increment the +x direction
                    if (j % 2 == 0) && ((index + y_size + 1) <= x_size*y_size) {
                        // if the y index is even, then its neighbor is in the -x direction
                        return Some(index + y_size);
                    } else
                    if (j % 2 == 1) && (index >= y_size) {
                        // if the y index is odd, then its neighbor is in the -x direction
                        return Some(index - y_size);
                    }
                }
            }
        },
    }
    return None;
}

#[derive(PartialEq)]
#[derive(Copy, Clone)]
pub enum SymmetryType {
    C3V,
    C4V,
    C6V,
}
#[derive(PartialEq)]
pub enum SignalType {
    SigStart,
    SigStop,
}

pub struct Driver {
    internal_lattice: Lattice,
    basis: RwLock<Array2<f64>>, // could be used for visual plotting purposes
    x_size: usize,
    y_size: usize,
    total_nodes: usize,
    sym_type: SymmetryType,
    exchange_constant: f64,
    num_threads: usize,
    div: usize,
    rem: usize,
    magnitiztion_threads_started: bool,
    cluster_threads_started: bool,
    mag_pool: ThreadPool<StaticParker<SmallThreadData>>,
    cluster_pool: ThreadPool<StaticParker<SmallThreadData>>,
    tx_psum: Sender<f64>,
    rx_psum: Receiver<f64>,
    tx_go_psum: Sender<SignalType>,
    rx_go_psum: Receiver<SignalType>,
    flip_index_vector: Arc<RwLock<Vec<usize>>>,
    tx_go_cluster: Sender<(SignalType, f64)>,
    rx_go_cluster: Receiver<(SignalType, f64)>,
    tx_done_cluster: Sender<usize>,
    rx_done_cluster: Receiver<usize>,
    tx_cluster_queue: Sender<(usize, f64)>,
    rx_cluster_queue: Receiver<(usize, f64)>,
}

impl Driver {
    pub fn new(exchange_constant: f64, x_size: usize, y_size: usize, sym_type: SymmetryType, inner_basis: Array2<f64>) -> Self { 
        assert!(x_size > 0, "size must be greater than zero!");
        assert!(y_size > 0, "size must be greater than zero!");
        let cpu_num = num_cpus::get();
        let (div, rem) = dividend_remainder(x_size*y_size, cpu_num);
        let mut cur_coord: Box<Array1<f64>>;
        let mut rng = thread_rng();
        let mut give_map: Vec<SpinNode> = vec![];
        let s1 = s![0_usize,..];
        let s2 = s![1_usize,..];
        match sym_type {
            // match the generation pattern to the type passed to sym_type
            SymmetryType::C4V => {
                // C4V generates every node including the starting node, incrementing by 1 x or y unit vector each
                // itteraton.
                for i in 0..x_size {
                    for j in 0..y_size {
                        let genval: f64 = rng.gen();
                        // randomization conditions would go here or somewhere near here.
                        cur_coord = Box::new((i as f64) * &inner_basis.slice(&s1) + (j as f64) * &inner_basis.slice(&s2));
                        // construct neighbors vector
                        let mut neighbors = vec![];
                        for modnum in 0..4 {
                            if let Some(valid_index) = indexmod(i+j*x_size, modnum, x_size, y_size, sym_type) {
                                neighbors.push(valid_index);
                            }
                        }
                        if genval >= 0.5 {
                            give_map.push(
                                SpinNode::cons_node(
                                    0.5, 
                                    array![i as f64, j as f64], 
                                    *cur_coord, 
                                    RwLock::new(neighbors)));
                        } else {
                            give_map.push(
                                SpinNode::cons_node(
                                    -0.5,
                                    array![i as f64, j as f64],
                                    *cur_coord,
                                    RwLock::new(neighbors)));
                        }
                    }
                }
            },
            SymmetryType::C3V => {
                // C3V generates a lattice like graphene, a hexagonal tile pattern where each vertex in the pattern is an
                // atom of carbon (or whatever-else exotic or mystical novelty have you).
                //
                // Starting with this building block that kinda looks like a bisected ladder warped beyond repair:
                //
                //                                              \_
                //                                              /
                //                                              \_
                //                                              /
                //
                // This part works out to taking the upwards facing basis vector and alternatly adding a node at the end
                // point at:
                //
                //          r_node_i = Sum on i of [b1 if i % 2 == 0 else b3 ], for i in [0,y_size) and y_size in N
                //
                // Next lets find what b3 is. Letting r(x,y) be the radius vector in 2D-Cartesian we want the rotation of the b1 vector by pi/2
                // about its midpoint:
                //
                //                          T(-b1*e_x/2, -b1*e_y/2)*R(pi/2)*T^-1(-b1*e_x/2, -b1*e_y/2)*b1
                //                      
                //              Where T(a, b) := [ [1,0,a],[0,1,b],[0,0,1] ], where T is in R^(3x3), r(x,y) is in R^2 
                //
                //              T^-1(-b1*e_x/2, -b1*e_y/2) = T(-1*-b1*e_x/2, -1*-b1*e_y/2) = T(b1*e_x/2, b1*e_y/2)
                //
                // =>                   T(-b1*e_x/2, -b1*e_y/2) * R(pi/2) * T(b1*e_x/2, b1*e_y/2)*b1
                // =>   [ [1,0,-b1_x/2],[0,1,-b1_y/2],[0,0,1] ] * [ [0,-1,0],[1,0,0],[0,0,1] ] * [ [1,0,b1_x/2],[0,1,b1_y/2],[0,0,1] ]*b1
                // =>               [ [0, -1, -(b1_x + b1_y)/2], [1, 0, (b1_x - b1_y)/2], [0, 0, 1] ]*b1
                //
                //                          A'Volià b3 := [-b1_x/2 - (3 b1_y)/2, (3 b1_x)/2 - b1_y/2, 1 ]
                //
                // Then we add the next node at:
                //                               r_node_j = rb1_node_i + b2
                //
                // Then we can add starting from the location pointed to by b1 + b2 + b3 a new node. From this starting point, we
                // do the following:
                //
                //                      r-inv_node_i = Sum on i of [b3 if i % 2 == 0 else b1 ]
                //                                          or equivlantly
                //                              Sum on i of [b1 if i % 2 == 1 else b3 ]
                // This gives us:
                //
                //                                              \_/
                //                                              / \
                //                                              \_/
                //                                              / \
                //
                // When this scheme is repeated, a lattice is created and we get (real end product looks better than comment):
                //
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //
                // Labeling of the nodes increments in the y direction by 1, starting from 0.

                let cur_basis = inner_basis.clone();
                let b1 = array![cur_basis[(0,0)], cur_basis[(0,1)]];
                let b2 = array![cur_basis[(1,0)], cur_basis[(1,1)]];
                let b3 = array![-b2[0], b2[1]];
                let mut cur_coord: Box<Array1<f64>> = Box::new(array![0., 0.]);

                for i in 0..x_size {
                    for j in 0..y_size {
                        let mut neighbors = vec![];
                        let cur_index = j + i * y_size;
                        // add the + and - y neighbors to the neighbors vector.
                        // this needs to be checked for all nodes.
                        if j > 0 {
                            // if the y index is greater than 0, we can add the -y node
                            neighbors.push(cur_index - 1);
                        }
                        if (j + 2) < y_size {
                            // if the y index is 2 less than the  
                            neighbors.push(cur_index + 1);
                        }
                        if i % 2 == 0 {
                            // when we are on an even increment the +x direction
                            if (j % 2 == 0) && (cur_index >= y_size) {
                               // if the y index is even, then its neighbor is in the -x direction
                               neighbors.push(cur_index - y_size);
                            } else
                            if (j % 2 == 1) && ((cur_index + y_size + 1) <= x_size*y_size) {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index + y_size);
                            }
                        } else {
                            // when we are on an odd increment the +x direction
                            if (j % 2 == 0) && ((cur_index + y_size + 1) <= x_size*y_size) {
                                // if the y index is even, then its neighbor is in the -x direction
                                neighbors.push(cur_index + y_size);
                            } else
                            if (j % 2 == 1) && (cur_index >= y_size) {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index - y_size);
                            }
                        }
                        // randomization conditions would go here or somewhere near here.
                        let genval: f64 = rng.gen();
                        if genval >= 0.5 {
                            give_map.push(
                                SpinNode::cons_node(
                                    0.5,
                                    array![i as f64, j as f64],
                                    *cur_coord.clone(),
                                    RwLock::new(neighbors)));
                        } else {
                            give_map.push(
                                SpinNode::cons_node(
                                    -0.5,
                                    array![i as f64, j as f64],
                                    *cur_coord.clone(),
                                    RwLock::new(neighbors)));
                        }
                        if i % 2 == 0 {
                            // when we are on an even increment of b2
                            if j % 2 == 0 {
                                // for even increments of the y vector
                                *cur_coord = *cur_coord.clone() + &b2;
                            } 
                            else {
                                *cur_coord = *cur_coord.clone() + &b3;
                            }
                        } else {
                            if j % 2 == 1 {
                                // for odd increments of the y vector
                                *cur_coord = *cur_coord.clone() + &b2;
                            } 
                            else {
                                *cur_coord = *cur_coord.clone() + &b3;
                            }
                        }
                    }
                    *cur_coord = array![&b1[0] + 2. * &b2[0], 0.];
                }
            },
            SymmetryType::C6V => {
                // C6V generates a lattice called a triangular lattice, as it can be made from alternating reflected triangles
                // or lines of zigzags tiled together (not to scale, the proper angle between b1 and b2 gives the lattice it's
                // hexagonal shape, sadly / and \ are not at the correct angle from the x axis :< ):
                //
                //      /\/\/\/\...             /\/\/\/\/\/\ 
                //          +                   \/\/\/\/\/\/
                //      \/\/\/\/...             /\/\/\/\/\/\...-> 
                //          +           =       \/\/\/\/\/\/ 
                //      /\/\/\/\...             /\/\/\/\/\/\ 
                //          +a                  \/\/\/\/\/\/
                //         ...                       ...
                //                                    |
                //                                    V
                //

                // let cur_basis = inner_basis.clone();
                // let b1 = array![cur_basis[(0,0)], cur_basis[(0,1)]];
                // let b2 = array![cur_basis[(1,0)], cur_basis[(1,1)]];
                // let b3 = array![-b2[0], b2[1]];
                // let mut cur_coord = array![0., 0.];
                // for i in 0..x_size {
                //     for j in 0..y_size {
                //         assert_eq!(i, i);  // TODO
                //         assert_eq!(j, j);
                //     }
                // }
            },
        }
        println!("Done!");
        let (tx_psum, rx_psum) = bounded(cpu_num);
        let (tx_go_psum, rx_go_psum) = bounded(cpu_num);
        // let (tx_flip_index, rx_flip_index) = unbounded();
        let (tx_go_cluster, rx_go_cluster) = bounded(cpu_num);
        let (tx_cluster_queue, rx_cluster_queue) = unbounded();
        let (tx_done_cluster, rx_done_cluster) = bounded(cpu_num);
        Self {
            internal_lattice: Lattice::new(give_map),
            basis: RwLock::new(inner_basis),
            x_size,
            y_size,
            total_nodes: x_size * y_size,
            sym_type,
            exchange_constant,
            num_threads: cpu_num,
            div,
            rem,
            magnitiztion_threads_started: false,
            cluster_threads_started: false,
            mag_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            cluster_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            tx_psum,
            tx_go_psum,
            rx_psum,
            rx_go_psum,
            flip_index_vector: Arc::new(RwLock::new(vec![])),
            tx_go_cluster,
            rx_go_cluster,
            tx_done_cluster,
            rx_done_cluster,
            tx_cluster_queue,
            rx_cluster_queue,

        }
    }

    fn start_magnitization_threads(&mut self) {
        if self.magnitiztion_threads_started != true {
            let n_jobs = self.num_threads.clone();
            // let (finished_tx, finished_rx) = channel();
            for i in 0..n_jobs {
                println!("starting magnitization thread {i}");
                let shared_data = self.internal_lattice.internal_vector.clone();
                let tx_psum = self.tx_psum.clone();
                let rx_go_psum = self.rx_go_psum.clone();
                let wolff_run = self.cluster_threads_started.clone();
                let range: std::ops::Range<usize>;
                if i != n_jobs - 1 {
                    range = (i * self.div)..((i+1) * self.div);
                } else {
                    // get the remainder in the last thread. We could smear it out amongst the threads
                    // but this is easier for all edge cases.
                    range = (i * self.div)..((i+1) * self.div + self.rem);
                }

                self.mag_pool.execute(move || {

                    loop {
                        let mut psum = 0.;
                        if let SignalType::SigStart = rx_go_psum.recv().unwrap() {
                            if let Ok(shared_vector) = shared_data.read().as_deref() {
                                for index in range.clone() {
                                    // SAFETY: bounds are explicitly known when range is generated
                                    // and garuntees a valid return.
                                    unsafe {
                                        if wolff_run {
                                            let mut write_lock = shared_vector.get_unchecked(index).marked.write().unwrap();
                                            write_lock.unmark();
                                        }
                                        psum += shared_vector.get_unchecked(index).get_spin();
                                    }
                                }
                            } else {
                                panic!("couldnt read shared data in thread!")
                            }
                            tx_psum.send(psum).unwrap();
                        } else {
                            println!("stopping thread");
                            return;
                        }
                    }

                });
            } // for i in 0..n_jobs
        }
        self.magnitiztion_threads_started = true;
    }

    fn start_cluster_threads(&mut self, ) {
        if self.cluster_threads_started != true {
            let n_jobs = self.num_threads.clone();
            // let (finished_tx, finished_rx) = channel();
            for i in 0..n_jobs {
                println!("starting cluster thread {i}");
                let shared_flip_vec = self.flip_index_vector.clone();
                let shared_data = self.internal_lattice.internal_vector.clone();
                let rx_go_cluster = self.rx_go_cluster.clone();
                let tx_cluster_queue = self.tx_cluster_queue.clone();
                let rx_cluster_queue = self.rx_cluster_queue.clone();
                let tx_done_cluster = self.tx_done_cluster.clone();
                self.cluster_pool.execute(move || {
                    let mut rng = thread_rng();
                    'outer : loop {
                        let sent_signal = rx_go_cluster.recv().unwrap();
                        match sent_signal.0 {
                            SignalType::SigStart => {
                                loop {
                                    if let Ok(shared_vector) = shared_data.read().as_deref() {
                                        // pop a value from the work queue if available, else move to waiting for next go signal.
                                        if let Ok(check_index) = rx_cluster_queue.try_recv() {
                                            // SAFETY: index is checked in calling function and nodes only store valid
                                            // indicies to other nodes.
                                            unsafe {
                                                let read_lock = shared_vector.get_unchecked(check_index.0)
                                                                                             .neighbors
                                                                                             .read()
                                                                                             .unwrap();
                                                'nbr_loop : for i in 0..read_lock.len() {
                                                    let nbr_index = *read_lock.get_unchecked(i);
                                                    let nbr_spin = shared_vector.get_unchecked(nbr_index).get_spin();
                                                    let mut marked_lock = shared_vector.get_unchecked(nbr_index).marked.write().unwrap();
                                                    if nbr_spin == check_index.1 {
                                                        if marked_lock.get_marked() == true {
                                                            continue 'nbr_loop
                                                        }
                                                        marked_lock.flip();
                                                        drop(marked_lock);
                                                        if rng.gen_range(0_f64..1_f64) < sent_signal.1 {
                                                            let mut write_lock = shared_flip_vec.write().unwrap();
                                                            write_lock.push(nbr_index);
                                                            tx_cluster_queue.send((nbr_index, nbr_spin)).unwrap();
                                                        }
                                                    }
                                                }
                                            }
                                        } else {
                                            // println!("thread done");
                                            // send done signal to main
                                            tx_done_cluster.send(1).unwrap();
                                            continue 'outer
                                        }
                                    } else {
                                        panic!("couldnt read shared data in thread!")
                                    }
                                }
                            },
                            SignalType::SigStop => {
                                println!("stopping thread");
                                return;
                            },
                        }
                    }

                });
            } // for i in 0..n_jobs
        }
        self.cluster_threads_started = true;
    }

    pub fn stop_threads(self) {
        // for _ in 0..self.num_threads {
        //     self.tx_go.send(SignalType::SigStop).unwrap();
        // }
        println!("I am not even going to bother with figuring out how to do this as the threads just wont die no matter what I do, so live with all the errors, though know that they are really just the OS complaining that the threads are still alive or were already killed. Not much I can do with my current understanding of this language.\n");
        // let mut count = 0;
        // loop {
        //     match self.mag_pool.shutdown_borrowed() {
        //         Ok(_) => break,
        //         Err(why) => {
        //             count += 1;
        //             println!("{}", why);
        //             if count == 4 {
        //                 break;
        //             }
        //             continue;
        //         },
        //     }
        // }
        if self.magnitiztion_threads_started {
            drop(self.tx_go_psum);
        }
        if self.cluster_threads_started {
            drop(self.tx_go_cluster);
        }
        println!("threads closed, or not. They may throw a panic, but this technically kills them so its all good. just ignore the incoming text wall...\n")
    }


    pub fn get_spin_at(&self, index: usize) -> f64 {
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            match value.get(index) {
            Some(node) => {node.get_spin()},
            None => 0.,
        }
        } else {
            panic!("couldnt get spin")
        }
    }

    fn flip_node_at(&self, index: usize) {
        if let Ok(mut lock) = self.internal_lattice.internal_vector.write() {
            match lock.get_mut(index) {
                Some(node) => node.flip_spin(),
                None => return,
            }
        }
    }
    
    pub fn get_magnitization(&self) -> f64 {
        for _ in 0..self.num_threads {
            self.tx_go_psum.send(SignalType::SigStart).unwrap();
        }
        self.rx_psum.iter().take(self.num_threads).fold(0., |a, b| a + b)
    }
    
    fn energy_worker<'a>(&'a self, range: Range<usize>) -> f64 {
        let thread_return_value = thread::scope(|scope| {
            let th = scope.spawn(move |_| {
                let mut psum = 0.;
                for i in range {
                    match self.sym_type {
                        SymmetryType::C3V => {
                            // 3 neighbors at most, so 0..3 is the range
                            for j in 0..3 {
                                if let Some(valid_index) = indexmod(i, j, self.x_size, self.y_size, self.sym_type) {
                                    psum += self.get_spin_at(valid_index);
                                }
                            }
                        },
                        SymmetryType::C4V => {
                            // 4 neighbors at most, so 0..4 is the range
                            for j in 0..4 {
                                if let Some(valid_index) = indexmod(i, j, self.x_size, self.y_size, self.sym_type) {
                                    psum += self.get_spin_at(valid_index);
                                }
                            }
                        },
                        SymmetryType::C6V => {
                            // 6 neighbors at most, so 0..6 is the range

                        },
                    }
                }
                psum
            });
            let psum = th.join().unwrap();
            psum
        }
        ).unwrap();
        return thread_return_value;
    }

    fn get_energy<'a>(&'a self) -> f64 {
        let mut energy = 0.;
        for i in 0usize..self.num_threads {
            let range: std::ops::Range<usize>;
                if i != self.num_threads - 1 {
                    range = (i * self.div)..((i+1) * self.div);
                } else {
                    // Get the remainder in the last thread spawned. We could smear it out amongst
                    // the other threads, but this is easier for all edge cases.
                    range = (i * self.div)..((i+1) * self.div + self.rem);
                }
            energy += self.energy_worker(range);
        }
        return energy;
    }

    pub fn spin_energy(&mut self, beta_list: Vec<f64>, times: usize, iteration_scheme: usize) {
        self.start_magnitization_threads();
        if iteration_scheme == 1 {
            self.start_cluster_threads();
        }
        let magnitization = self.get_magnitization();
        
        let mut energy = self.get_energy();
        println!("The initial energy is {}", energy);
        let mut energy_vec: Vec<f64> = vec![];
        
        let mut ms: Vec<f64> = vec![magnitization, magnitization.powf(2.)];
        let mut es: Vec<f64> = vec![energy, energy.powf(2.)];

        let mut m_vec: Vec<f64> = vec![];
        let mut e_vec: Vec<f64> = vec![];
        let mut c_vec: Vec<f64> = vec![];
        let mut x_vec: Vec<f64> = vec![];
        let n1 = (self.total_nodes as f64).powf(2.)*(times as f64);
        let n2 = (self.total_nodes as f64).powf(2.)*(times as f64).powf(2.);
        
        let len_beta = beta_list.len();
        let tot_time: u64 = times as u64 * len_beta as u64;

        let bar1 = ProgressBar::new(tot_time);

        for beta_val in beta_list {
            // let mut run_avg: Vec<f64> = vec![];
            // let mut cur_run_sum = 0.; 
            for _ in 0..times {
                if iteration_scheme == 0 {
                    energy = self.metropolis_iter(&beta_val, energy);
                } else
                if iteration_scheme == 1 {
                    energy = self.wolff_iter(&beta_val, energy);
                } else {
                    panic!("invalid option! expected 0 or 1, got {}", iteration_scheme);
                }
                let cur = self.get_magnitization();
                ms[0] += magnitization.abs();
                ms[1] += magnitization.powf(2.);
                es[0] += energy;
                es[1] += energy.powf(2.);
                // cur_run_sum += cur;
                m_vec.push(cur / self.total_nodes as f64);
                e_vec.push(energy / self.total_nodes as f64);  // TODO add in actual value
                c_vec.push(es[1] / n1 - es[0].powf(2.) / n2 * beta_val.powf(2.));
                x_vec.push(ms[1] / n1 - ms[0].powf(2.) / n2 * beta_val);
                bar1.inc(1);
            }
            energy_vec.push(energy.clone());
        }
        bar1.abandon();
        println!("spin_energy finished! Writing data to file, 少々お待ちして下さい");
        let bar2 = ProgressBar::new_spinner();
        
        // setup saving the output file
        let path = Path::new("mag_data.dat");
        let display = path.display();
        // Write output data to a file
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            },
        };
        // start writing data
        for item in m_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    bar2.inc(1);
                    continue;
                },
            }
        }
        
        // setup saving the output file
        let path = Path::new("energy_data.dat");
        let display = path.display();
        // Write output data to a file
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            },
        };
        for item in e_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    bar2.inc(1);
                    continue;
                },
            }
        }
        
        // setup saving the output file
        let path = Path::new("susceptibility_data.dat");
        let display = path.display();
        // Write output data to a file
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            },
        };
        for item in x_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    bar2.inc(1);
                    continue;
                },
            }
        }

        // setup saving the output file
        let path = Path::new("heatcap_data.dat");
        let display = path.display();
        // Write output data to a file
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            },
        };
        for item in c_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    bar2.inc(1);
                    continue;
                },
            }
        }
        bar2.abandon();
        println!("Sucessfully wrote file! Exiting...");
    }

    fn wolff_iter(&mut self, beta: &f64, energy: f64) -> f64 {
        //
        // Purpose
        // -------
        // Mutates the self.LinkedLat lattice of spins by one Iteration of
        // the Wolff Algorithm.
        //
        // if the required threads are not alive already, launch them.
        let balcond = 1. - consts::E.powf(-2.*beta*self.exchange_constant);
        // println!("{balcond}");
        let mut rngx = thread_rng();
        let mut rngy = thread_rng();
        let mut rng_flip = thread_rng();
        let mut x_y: usize;
        let mut target_spin: f64;
        let mut delta_energy: f64 = 0.;
        // let (tx_cluster, rx_cluster) = unbounded();

        // pick random point
        'node_selection : loop {
            x_y = rngx.gen_range(0..self.x_size) + rngy.gen_range(0..self.y_size) * self.x_size;
            target_spin = self.get_spin_at(x_y);
            // println!("{}", &x_y);
            if target_spin != 0. { break 'node_selection; }
            else { continue 'node_selection; }
        }

        // push the random node neighbors index to the work queue
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            if let Some(node) = value.get(x_y) {
                let read_lock_nbrs = node.neighbors.read().unwrap();
                for i in 0..read_lock_nbrs.len() {
                    // SAFETY: bounds checked in Driver::new, garaunteed valid return.
                    unsafe {
                        let nbr_index = *read_lock_nbrs.get_unchecked(i);
                        if self.get_spin_at(nbr_index) == target_spin {
                            // if the index is the same as the randomly picked node, add it to the
                            // queue.
                            let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap();
                            write_lock_flip_index_vector.push(nbr_index);
                            self.tx_cluster_queue.send((nbr_index, target_spin)).unwrap();
                            drop(write_lock_flip_index_vector);
                        }
                    }
                }
            }
        }

        // spin up threads for generating a cluster flip and then wait for them to finish.
        // println!("sending go signal to cluster threads");
        for _ in 0..self.num_threads {
            self.tx_go_cluster.send((SignalType::SigStart, balcond)).unwrap();
        }
        // println!("waiting for threads to finish...");
        while self.rx_done_cluster.is_full() == false {
            stdth::sleep(time::Duration::from_micros(10));
        }
        // println!("Clearing queue");
        for _ in 0..self.num_threads {
            // clear the queue
            _ = self.rx_done_cluster.recv().unwrap();
        }

        // flip spins and calculate the path integral value of the change in energy dE
        let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap(); 
        loop {
            // while the flip queue is not empty, process nodes
            if let Some(cur_flip_index) = write_lock_flip_index_vector.pop() {
                if rng_flip.gen_range(0_f64..1_f64) < balcond {
                    let nbr_energy = self.get_neighbors_spin_sum(cur_flip_index);
                    delta_energy += 2.*target_spin*nbr_energy;
                    self.flip_node_at(cur_flip_index);
                }
            } else {
                break
            }
        }
        drop(write_lock_flip_index_vector);
        return energy + delta_energy;
    }

    fn metropolis_iter(&mut self, beta: &f64, energy: f64) -> f64 {
        let mut rngx = thread_rng();
        let mut rngy = thread_rng();
        let mut rng_flip = thread_rng();
        let mut x_y: usize;
        let mut target_spin: f64;
        // Evolves the lattice by one iteration.
        'node_selection : loop {
            x_y = rngx.gen_range(0..self.x_size) + rngy.gen_range(0..self.y_size) * self.x_size;
            target_spin = self.get_spin_at(x_y);
            // println!("{}", &x_y);
            if target_spin != 0. { break 'node_selection; }
            else { continue 'node_selection; }
        }

        let mut delta_energy: f64;
        // calculate the sum of energy of the neighbors
        let nbr_energy = self.get_neighbors_spin_sum(x_y);
        delta_energy = 2.*target_spin*nbr_energy;
        if delta_energy <= 0. {
            self.flip_node_at(x_y);
            // dE = 2.*target_spin*nbr_E;
        }
        else if rng_flip.gen_range(0_f64..1_f64) < consts::E.powf(-beta*self.exchange_constant*delta_energy) {
            self.flip_node_at(x_y);
            // dE = 2.*target_spin*nbr_E;
        }
        else {
            delta_energy = 0.;
        }

        return energy + delta_energy;
    }

    fn get_neighbors_spin_sum(&self, x_y: usize) -> f64 {
        let mut nbr_energy = 0.;
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            if let Some(node) = value.get(x_y) {
                let read_lock_nbrs = node.neighbors.read().unwrap();
                for i in 0..read_lock_nbrs.len() {
                    // SAFETY: bounds checked in random node selection, garaunteed valid return
                    unsafe {
                        nbr_energy += self.get_spin_at(*read_lock_nbrs.get_unchecked(i));
                    }
                }
            }
        }
        nbr_energy
    }
}
    