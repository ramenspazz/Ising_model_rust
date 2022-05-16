use std::{thread as stdth, time};
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
// use core::ops::Range;
use indicatif::ProgressBar;
use std::f64::consts;
use rand::thread_rng;
use rand::Rng;
use std::sync::{RwLock, Arc};
use executors::parker::SmallThreadData;
use ndarray::prelude::*;
use executors::Executor;
use executors::crossbeam_workstealing_pool;
use executors::parker::StaticParker;
use crate::atomistic_simulation::crossbeam_workstealing_pool::ThreadPool;
use crate::dividend_remainder;
use crate::SignalContainer::*;
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
    energy_threads_started: bool,
    mag_pool: ThreadPool<StaticParker<SmallThreadData>>,
    energy_pool: ThreadPool<StaticParker<SmallThreadData>>,
    cluster_pool: ThreadPool<StaticParker<SmallThreadData>>,
    magnitization_psum_signaler: SignalContainer<f64>,
    magnitization_go_stop_signaler: SignalContainer<SignalType>,
    energy_psum_signaler: SignalContainer<f64>,
    energy_go_stop_signaler: SignalContainer<SignalType>,
    cluster_go_stop_signaler: SignalContainer<(SignalType, f64)>,
    cluster_done_signaler: SignalContainer<usize>,
    cluster_queue_signaler: SignalContainer<(usize, f64)>,
    flip_index_vector: Arc<RwLock<Vec<usize>>>,
}

impl Driver {
    pub fn new(exchange_constant: f64, x_size: usize, y_size: usize, sym_type: SymmetryType, inner_basis: Array2<f64>, spin_up_chance: f64, spin_unit: f64) -> Self { 
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
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        // randomization conditions would go here or somewhere near here.
                        cur_coord = Box::new((i as f64) * &inner_basis.slice(&s1) + (j as f64) * &inner_basis.slice(&s2));
                        // construct neighbors vector
                        let mut neighbors = vec![];
                        for modnum in 0..4 {
                            if let Some(valid_index) = indexmod(i+j*x_size, modnum, x_size, y_size, sym_type) {
                                neighbors.push(valid_index);
                            }
                        }
                        if genval >= spin_up_chance {
                            give_map.push(
                                SpinNode::cons_node(
                                    spin_unit, 
                                    array![i as f64, j as f64], 
                                    *cur_coord, 
                                    RwLock::new(neighbors)));
                        } else {
                            give_map.push(
                                SpinNode::cons_node(
                                    -spin_unit,
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
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        if genval >= spin_up_chance {
                            give_map.push(
                                SpinNode::cons_node(
                                    spin_unit,
                                    array![i as f64, j as f64],
                                    *cur_coord.clone(),
                                    RwLock::new(neighbors)));
                        } else {
                            give_map.push(
                                SpinNode::cons_node(
                                    -spin_unit,
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

        let psum_signaler = SignalContainer::new(cpu_num);
        let psum_go_signaler = SignalContainer::new(cpu_num);
        let energy_signlaer = SignalContainer::new(cpu_num);
        let energy_go_signaler = SignalContainer::new(cpu_num);
        let cluster_go_signaler = SignalContainer::new(cpu_num);
        let cluster_done_signaler = SignalContainer::new(cpu_num);
        let cluster_queue_signaler = SignalContainer::new(0);

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
            energy_threads_started: false,
            mag_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            energy_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            cluster_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            magnitization_psum_signaler: psum_signaler,
            magnitization_go_stop_signaler: psum_go_signaler,
            energy_psum_signaler: energy_signlaer,
            energy_go_stop_signaler: energy_go_signaler,
            cluster_go_stop_signaler: cluster_go_signaler,
            cluster_done_signaler,
            cluster_queue_signaler,
            flip_index_vector: Arc::new(RwLock::new(vec![])),
        }
    }

    pub fn save_state(&self, fname: &str) {
        self.internal_lattice.export_state_to_file(fname);
    }

    pub fn load_state(&mut self, fname: &str) {
        self.internal_lattice.load_state_from_file(fname);
    }

    fn start_energy_threads(&mut self) {
        if self.energy_threads_started != true {
            let n_jobs = self.num_threads.clone();

            for i in 0..n_jobs {
                let shared_data = self.internal_lattice.internal_vector.clone();
                let energy_psum_signaler = self.energy_psum_signaler.clone();
                let energy_go_stop_signaler = self.energy_go_stop_signaler.clone();

                let range: std::ops::Range<usize>;
                if i != self.num_threads - 1 {
                    println!("starting energy thread {i} with range {} to {}", (i * self.div), ((i+1) * self.div - 1));
                    range = (i * self.div)..((i+1) * self.div - 1);
                } else {
                    // Get the remainder in the last thread spawned. We could smear it out amongst
                    // the other threads, but this is easier for all edge cases.
                    println!("starting energy thread {i} with range {} to {}", (i * self.div), ((i+1) * self.div + self.rem - 1));
                    range = (i * self.div)..((i+1) * self.div + self.rem - 1);
                }

                self.energy_pool.execute(move || {

                    loop {
                        let mut energy_psum = 0.;
                        if let SignalType::SigStart = energy_go_stop_signaler.recv().unwrap() {
                            if let Ok(shared_vector) = shared_data.read().as_deref() {
                                unsafe { for index in range.clone() {
                                    let read_lock_nbrs = shared_vector.get_unchecked(index).neighbors.read().unwrap();
                                    let target_spin = shared_vector.get_unchecked(index).get_spin();
                                    for ith_nbr in 0..read_lock_nbrs.len() {
                                        // SAFETY: bounds are explicitly known when range is generated
                                        // and garuntees a valid return.
                                        let nbr_index = read_lock_nbrs.get_unchecked(ith_nbr);
                                        energy_psum += target_spin * shared_vector.get_unchecked(*nbr_index).get_spin();
                                    }
                                }}
                            } else {
                                panic!("couldnt read shared data in thread!")
                            }
                            match energy_psum_signaler.send(-energy_psum) {
                                Ok(_) => continue,
                                Err(_) => return,
                            }
                        } else {
                            println!("stopping thread");
                            return;
                        }
                    }

                });
            } // for i in 0..n_jobs
        }
        self.energy_threads_started = true;
    }

    fn start_magnitization_threads(&mut self) {
        if self.magnitiztion_threads_started != true {
            let n_jobs = self.num_threads.clone();

            for i in 0..n_jobs {
                let shared_data = self.internal_lattice.internal_vector.clone();
                let magnitization_psum_signaler = self.magnitization_psum_signaler.clone();
                let magnitization_go_stop_signaler = self.magnitization_go_stop_signaler.clone();

                let wolff_run = self.cluster_threads_started.clone();
                let range: std::ops::Range<usize>;

                if i != self.num_threads - 1 {
                    println!("starting magnitization thread {i} with range {} to {}", (i * self.div), ((i+1) * self.div - 1));
                    range = (i * self.div)..((i+1) * self.div - 1);
                } else {
                    // Get the remainder in the last thread spawned. We could smear it out amongst
                    // the other threads, but this is easier for all edge cases.
                    println!("starting magnitization thread {i} with range {} to {}", (i * self.div), ((i+1) * self.div + self.rem - 1));
                    range = (i * self.div)..((i+1) * self.div + self.rem - 1);
                }

                self.mag_pool.execute(move || {

                    loop {
                        let mut psum = 0.;
                        if let Ok(SignalType::SigStart) = magnitization_go_stop_signaler.recv() {
                            if let Ok(shared_vector) = shared_data.read().as_deref() {
                                for index in range.clone() {
                                    // SAFETY: bounds are explicitly known when range is generated
                                    // and garuntees a valid return.
                                unsafe { if wolff_run == true {
                                        let mut write_lock = shared_vector.get_unchecked(index).marked.write().unwrap();
                                        write_lock.unmark();
                                    }
                                    psum += shared_vector.get_unchecked(index).get_spin(); }
                                }
                            } else {
                                panic!("couldnt read shared data in thread!")
                            }
                            match magnitization_psum_signaler.send(psum) {
                                Ok(_) => continue,
                                Err(_) => return,
                            }
                        } else {
                            println!("Stop stignal received, closing thread.");
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
            for i in 0..n_jobs {
                println!("starting cluster thread {i}");
                let shared_flip_vec = self.flip_index_vector.clone();
                let shared_data = self.internal_lattice.internal_vector.clone();
                let cluster_go_stop_signaler = self.cluster_go_stop_signaler.clone();
                let cluster_queue_signaler = self.cluster_queue_signaler.clone();
                let cluster_done_signaler = self.cluster_done_signaler.clone();

                self.cluster_pool.execute(move || {
                    let mut rng = thread_rng();
                    'outer : loop {
                        if let Ok(sent_signal) = cluster_go_stop_signaler.recv() {
                            if let SignalType::SigStart = sent_signal.0 {
                                loop {
                                    // pop a value from the work queue if available, else move to waiting for next go signal.
                                    if let (Ok(read_lock_internal_vector), Ok(check_index)) = (shared_data.read().as_deref(), cluster_queue_signaler.try_recv()) {
                                        // SAFETY: index is checked in calling function and nodes only store valid
                                        // indicies to other nodes.
                                    unsafe { let read_lock_shared_vector = read_lock_internal_vector.get_unchecked(check_index.0)
                                                                .neighbors
                                                                .read()
                                                                .unwrap();
                                        'nbr_loop : for i in 0..read_lock_shared_vector.len() {
                                            let nbr_index = *read_lock_shared_vector.get_unchecked(i);
                                            let nbr_spin = read_lock_internal_vector.get_unchecked(nbr_index).get_spin();
                                            // if the two spins are the same
                                            if nbr_spin == check_index.1 {
                                                // check that the node was not previously processed and marked
                                                if let Ok(read_lock_marked) = read_lock_internal_vector.get_unchecked(nbr_index).marked.try_read() {
                                                    if read_lock_marked.get_marked() == true {
                                                        drop(read_lock_marked);
                                                        continue 'nbr_loop
                                                    }
                                                } else { continue 'nbr_loop }
                                                // if the node was not processed, mark it now
                                                if let Ok(mut write_lock_marked) = read_lock_internal_vector.get_unchecked(nbr_index).marked.try_write() {
                                                    write_lock_marked.mark();
                                                } else { continue 'nbr_loop }
                                                if rng.gen_range(0_f64..1_f64) < sent_signal.1 {
                                                    let mut write_lock = shared_flip_vec.write().unwrap();
                                                    write_lock.push(nbr_index);
                                                    cluster_queue_signaler.send((nbr_index, check_index.1)).unwrap();
                                                }
                                            }
                                        }}
                                    } else {
                                        // send done signal to main
                                        cluster_done_signaler.send(1).unwrap();
                                        continue 'outer
                                    }
                                }
                            } else {
                                println!("Stop stignal received, closing thread.");
                                return;
                            }
                        }
                    }

                });
            } // for i in 0..n_jobs
        }
        self.cluster_threads_started = true;
    }

    pub fn stop_threads(self) {
        if self.magnitiztion_threads_started {
            self.magnitization_go_stop_signaler.close_channel();
        }
        if self.energy_threads_started {
            self.energy_go_stop_signaler.close_channel();
        }
        if self.cluster_threads_started {
            self.cluster_go_stop_signaler.close_channel();
        }
        println!("Threads closed!\n")
    }

    pub fn get_spin_at(&self, index: usize) -> f64 {
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            match value.get(index) {
            Some(node) => {node.get_spin()},
            None => 0.,
        }
        } else {
            panic!("Couldnt get spin")
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

    pub fn get_spin_energy(&self) -> (f64, f64) {
        return (self.get_magnitization(), self.get_energy());
    }
    
    pub fn get_magnitization(&self) -> f64 {
        let mut mag = 0.;
        for _ in 0..self.num_threads {
            self.magnitization_go_stop_signaler.send(SignalType::SigStart).unwrap();
        }
        for _ in 0..self.num_threads {
            if let Ok(psum) = self.magnitization_psum_signaler.recv() {
                mag += psum;
            } else { panic!("Channel is already closed, something went wrong!") }
        }
        mag
    }

    pub fn get_energy(&self) -> f64 {
        let mut energy = 0.;
        for _ in 0..self.num_threads {
            self.energy_go_stop_signaler.send(SignalType::SigStart).unwrap();
        }
        for _ in 0..self.num_threads {
            if let Ok(psum) = self.energy_psum_signaler.recv() {
                energy += psum;
            } else { panic!("Channel is already closed, something went wrong!") }
        }
        energy
    }
    

    fn save_files(&self, m_vec: Vec<f64>, e_vec: Vec<f64>, x_vec: Vec<f64>, c_vec: Vec<f64>) {
        println!("Writing data to file, 少々お待ちして下さい");
        let file_progress = ProgressBar::new_spinner();
        
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
                    file_progress.inc(1);
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
                    file_progress.inc(1);
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
                    file_progress.inc(1);
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
                    file_progress.inc(1);
                    continue;
                },
            }
        }
        file_progress.abandon();
        println!("Sucessfully wrote files!");
    }

    pub fn spin_energy(&mut self, beta_list: Vec<f64>, times: usize, iteration_scheme: usize) {
        self.start_magnitization_threads();
        self.start_energy_threads();
        if iteration_scheme == 1 {
            self.start_cluster_threads();
        }
        
        let (initial_magnitization, initial_energy) = self.get_spin_energy();
        println!("\nThe initial energy is {}, and the initial magnitization is {}.\n", initial_energy, initial_magnitization);
        
        let mut energy: f64;
        let mut magnitization: f64;
        let mut ms: Vec<f64> ;
        let mut es: Vec<f64>;
        
        let mut m_vec: Vec<f64> = vec![];
        let mut e_vec: Vec<f64> = vec![];
        let mut c_vec: Vec<f64> = vec![];
        let mut x_vec: Vec<f64> = vec![];
        let n1 = (self.total_nodes as f64)*(times as f64);
        let n2 = (self.total_nodes as f64)*(times as f64).powf(2.);
        
        let len_beta = beta_list.len();
        let tot_time: u64 = times as u64 * len_beta as u64;
        
        let bar1 = ProgressBar::new(tot_time);
        
        // load the initial state at tbe begining of each new beta
        for beta_val in beta_list {
            match self.sym_type {
                SymmetryType::C3V => {
                    self.internal_lattice.load_state_from_file("c3v.dat");
                },
                SymmetryType::C4V => {
                    self.internal_lattice.load_state_from_file("c4v.dat");
                },
                SymmetryType::C6V => {
                    self.internal_lattice.load_state_from_file("c6v.dat");
                },
                
            }
            
            magnitization = initial_magnitization;
            energy = initial_energy;
            let (mut d_energy, mut d_mag);
            ms = vec![0., 0.];
            es = vec![0., 0.];

            for _ in 0..times {
                // preform an iteration scheme of metropolis or wolff
                if iteration_scheme == 0 {
                    (d_energy, d_mag) = self.metropolis_iter(&beta_val);
                } else
                if iteration_scheme == 1 {
                    (d_energy, d_mag) = self.wolff_iter(&beta_val);
                } else {
                    panic!("invalid option! expected 0 or 1, got {}", iteration_scheme);
                }
                // (magnitization, energy) = self.get_spin_energy();
                // Beginning the data collection after a few iterations gives better
                // overall graphs becuase the system has had time to relax
                // if cur_t > (0.375 * times as f64) as usize {
                // }
                magnitization += d_mag;
                energy += d_energy;
                ms[0] += magnitization;
                ms[1] += magnitization.powf(2.);
                es[0] += energy;
                es[1] += energy.powf(2.);
                bar1.inc(1);
            }
            m_vec.push(ms[0] / n1);
            e_vec.push(es[0] / n1);
            c_vec.push((es[1] / n1 - es[0].powf(2.) / n2) * beta_val.powf(2.));
            x_vec.push((ms[1] / n1 - ms[0].powf(2.) / n2) * beta_val);
        } // for beta_val in beta_list end
        bar1.abandon();
        println!("spin_energy finished!");
        
        // save files
        if true {
            self.save_files(m_vec, e_vec, x_vec, c_vec);
        }
    }
    
    /// Purpose
    /// -------
    /// Mutates the self.LinkedLat lattice of spins by one Iteration of
    /// the Wolff Algorithm.
    ///
    /// if the required threads are not alive already, launch them.
    fn wolff_iter(&mut self, beta: &f64) -> (f64, f64) {
        let balance_condition = 1. - consts::E.powf(-20.*beta*self.exchange_constant);
        let mut rngspin = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;
        let mut delta_energy: f64 = 0.;
        let mut delta_mag: f64 = 0.;

        // select a random node
        'node_selection : loop {
            target_index = rngspin.gen_range(0..(self.x_size*self.y_size - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. { break 'node_selection; }
            else { continue 'node_selection; }
        }

        // add target node to flip list
        // SAFETY: target_index was verified in get_spin_at already
        unsafe {
            let mut write_lock_internal_vector = self.internal_lattice.internal_vector.write().unwrap();
            let target_node = write_lock_internal_vector.get_unchecked_mut(target_index);
            target_node.marked
                .write()
                .unwrap()
                .mark();
            drop(write_lock_internal_vector);
        }

        // push the random node neighbors index to the work queue
        if let Ok(read_lock_internal_vector) = self.internal_lattice.internal_vector.read() {
            if let Some(target_node) = read_lock_internal_vector.get(target_index) {
                let read_lock_nbrs = target_node.neighbors.read().unwrap();
                let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap();
                for i in 0..read_lock_nbrs.len() {
                // SAFETY: bounds checked in Driver::new, garaunteed valid return.
                unsafe { let nbr_index = *read_lock_nbrs.get_unchecked(i);
                    if self.get_spin_at(nbr_index) == target_spin {
                        // if the spin is the same as the randomly picked node, add it to the
                        // queue.
                        if i == 0 {
                            write_lock_flip_index_vector.push(target_index);
                        }
                        write_lock_flip_index_vector.push(nbr_index);
                        self.cluster_queue_signaler.send((nbr_index, target_spin)).unwrap();
                    }}
                }
                drop(write_lock_flip_index_vector);
            }
        }

        // spin up threads for generating a cluster flip and then wait for them to finish.
        // println!("sending go signal to cluster threads");
        for _ in 0..self.num_threads {
            self.cluster_go_stop_signaler.send((SignalType::SigStart, balance_condition)).unwrap();
        }
        // println!("waiting for threads to finish...");
        while self.cluster_done_signaler.is_full() == false {
            stdth::sleep(time::Duration::from_micros(1));
        }
        // println!("Clearing queue");
        for _ in 0..self.num_threads {
            // clear the queue
            _ = self.cluster_done_signaler.recv().unwrap();
        }
        
        // flip spins and calculate the path integral value of the change in energy dE
        let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap(); 
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();
        loop {
            // while the flip queue is not empty, process nodes
            if let Some(cur_flip_index) = write_lock_flip_index_vector.pop() {
        // SAFETY: cur_flip_index will always point to a valid node
        unsafe {let cur_target = read_lock_internal_vector.get_unchecked(cur_flip_index);
                let cur_target_nbrs_lock = cur_target.neighbors.read().unwrap();
                for flip_node_nbrs_index in cur_target_nbrs_lock.iter() {
                    delta_energy += cur_target.get_spin() * self.get_neighbors_spin_sum(*flip_node_nbrs_index);
                }}
                delta_energy += -target_spin * self.get_neighbors_spin_sum(cur_flip_index);
                delta_mag += -target_spin;
                } else {
                    break
                }
            }
            drop(write_lock_flip_index_vector);
            println!("{}, {}", delta_energy, delta_mag);
            return (delta_energy, delta_mag);
    }

    /// Evolves the lattice by one iteration using the metropolis-hastings scheme.
    fn metropolis_iter(&mut self, beta: &f64) -> (f64, f64) {
        let mut rngspin = thread_rng();
        let mut rng_flip = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;

        // select a random node
        'node_selection : loop {
            target_index = rngspin.gen_range(0..(self.x_size*self.y_size - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. { break 'node_selection; }
            else { continue 'node_selection; }
        }

        let mut delta_energy: f64 = 0.;
        let delta_mag: f64;
        // calculate dE of the proposed spin flip at x_y
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();
        // SAFETY: x_y was checked already to be a valid node
        unsafe {let cur_target = read_lock_internal_vector.get_unchecked(target_index);
            let cur_target_nbrs_lock = cur_target.neighbors.read().unwrap();
            for flip_node_nbrs_index in cur_target_nbrs_lock.iter() {
                delta_energy += cur_target.get_spin() * self.get_neighbors_spin_sum(*flip_node_nbrs_index);
            }
            delta_energy += -target_spin * self.get_neighbors_spin_sum(target_index);
        }
        drop(read_lock_internal_vector);

                
        if delta_energy > 0. && (rng_flip.gen_range(0_f64..1_f64) < consts::E.powf(-beta*self.exchange_constant*delta_energy)) {
            self.flip_node_at(target_index);
            delta_mag = -target_spin;
        }
        else if delta_energy <= 0. {
            self.flip_node_at(target_index);
            delta_mag = -target_spin;
        }
        else {
            delta_energy = 0.;
            delta_mag = 0.;
        }
        return (delta_energy, delta_mag);
    }

    fn get_neighbors_spin_sum(&self, x_y: usize) -> f64 {
        let mut nbr_energy = 0.;
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            if let Some(node) = value.get(x_y) {
                // let target_spin = node.get_spin();
                let read_lock_nbrs = node.neighbors.read().unwrap();
                for ith_nbr_index in 0..read_lock_nbrs.len() {
                    // SAFETY: bounds checked in random node selection, garaunteed valid return
                    unsafe { nbr_energy += self.get_spin_at(*read_lock_nbrs.get_unchecked(ith_nbr_index)); }
                }
            }
        }
        return nbr_energy;
    }
}
    