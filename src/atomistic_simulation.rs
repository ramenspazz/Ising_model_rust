use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use std::{thread as stdth, time};
// use core::ops::Range;
use crate::{dividend_remainder, get_input_as_usize};
use crate::lat_node::*;
use crate::lattice_structure::Lattice;
use crate::signal_container::*;
use indicatif::ProgressBar;
use ndarray::prelude::*;
use rand::thread_rng;
use rand::Rng;
use rayon::prelude::*;
use std::f64::consts;
use std::sync::{Arc, RwLock};

fn indexmod(
    index: usize,
    modnum: usize,
    x_size: usize,
    y_size: usize,
    sym_type: SymmetryType,
) -> Option<usize> {
    let i = (index % x_size) as usize;
    let j = (index / x_size) as usize;
    match sym_type {
        SymmetryType::C3V => {
            if modnum == 0 {
                if j > 0 {
                    // if the y index is greater than 0, we can add the -y node
                    return Some(index - 1);
                }
            } else if modnum == 1 {
                if (j + 2) < y_size {
                    // if the y index is 2 less than the
                    return Some(index + 1);
                }
            } else if modnum == 2 {
                if i % 2 == 0 {
                    // when we are on an even increment the +x direction
                    if (j % 2 == 0) && (index >= y_size) {
                        // if the y index is even, then its neighbor is in the -x direction
                        return Some(index - y_size);
                    } else if (j % 2 == 1) && ((index + y_size + 1) <= x_size * y_size) {
                        // if the y index is odd, then its neighbor is in the -x direction
                        return Some(index + y_size);
                    }
                } else {
                    // when we are on an odd increment the +x direction
                    if (j % 2 == 0) && ((index + y_size + 1) <= x_size * y_size) {
                        // if the y index is even, then its neighbor is in the -x direction
                        return Some(index + y_size);
                    } else if (j % 2 == 1) && (index >= y_size) {
                        // if the y index is odd, then its neighbor is in the -x direction
                        return Some(index - y_size);
                    }
                }
            }
        }
        SymmetryType::C4V => {
            // we are using rc cola indexing (row-column) such that 'i' is the row, and 'j' is the column.
            if modnum == 0 {
                // +x
                if (i + 1) < x_size {
                    return Some(index + 1);
                } else {
                    return None;
                }
            } else if modnum == 1 {
                // +y
                if j + 1 < y_size {
                    // println!("index {} + xsize < ysize", &index);
                    return Some(index + x_size);
                } else {
                    return None;
                }
            } else if modnum == 2 {
                // -x
                if i > 0 {
                    // println!("index {} - 1 > 0", &index);
                    return Some(index - 1);
                } else {
                    return None;
                }
            } else if modnum == 3 {
                // -y
                if j > 0 {
                    return Some(index - x_size);
                } else {
                    return None;
                }
            }
        }
        SymmetryType::C6V => {}
    }
    None
}

#[derive(PartialEq, Copy, Clone)]
pub enum SymmetryType {
    C3V,
    C4V,
    C6V,
}

pub struct Driver {
    internal_lattice: Lattice,
    basis: RwLock<Array2<f64>>, // could be used for visual plotting purposes
    spin_unit: f64,
    x_size: usize,
    y_size: usize,
    total_nodes: usize,
    sym_type: SymmetryType,
    exchange_constant: f64,
    num_threads: usize,
    div: usize,
    rem: usize,
    cluster_threads_started: bool,
    energy_threads_started: bool,
    // mag_pool: ThreadPool<StaticParker<SmallThreadData>>,
    energy_pool: rayon::ThreadPool,
    cluster_pool: rayon::ThreadPool,
    // magnitization_psum_signaler: SignalContainer<f64>,
    energy_psum_signaler: SignalContainer<f64>,
    cluster_go_stop_signaler: SignalContainer<(SignalType, f64, f64)>,
    cluster_done_signaler: SignalContainer<usize>,
    cluster_queue_signaler: SignalContainer<(usize, f64)>,
    flip_index_vector: Arc<RwLock<Vec<usize>>>,
    touched_index_vec: Arc<RwLock<Vec<usize>>>,
    fname: String,
}

impl Driver {
    pub fn new(
        exchange_constant: f64,
        x_size: usize,
        y_size: usize,
        sym_type: SymmetryType,
        inner_basis: Array2<f64>,
        spin_up_chance: f64,
        spin_unit: f64,
        fname: String,
    ) -> Self {
        assert!(x_size > 0, "size must be greater than zero!");
        assert!(y_size > 0, "size must be greater than zero!");
        let cpu_num = num_cpus::get();
        let (div, rem) = dividend_remainder(x_size * y_size, cpu_num);
        let mut cur_coord: Box<Array1<f64>>;
        let mut rng = thread_rng();
        let mut give_map: Vec<SpinNode> = vec![];
        let s1 = s![0_usize, ..];
        let s2 = s![1_usize, ..];
        match sym_type {
            // match the generation pattern to the type passed to sym_type
            SymmetryType::C4V => {
                // C4V generates every node including the starting node, incrementing by 1 x or y unit vector each
                // itteraton.
                for j in 0..y_size {
                    for i in 0..x_size {
                        let cur_index = i + j * x_size;
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        // randomization conditions would go here or somewhere near here.
                        cur_coord = Box::new(
                            (i as f64) * &inner_basis.slice(&s1)
                                + (j as f64) * &inner_basis.slice(&s2),
                        );
                        // construct neighbors vector
                        let mut neighbors = vec![];
                        for modnum in 0..4 {
                            let result = indexmod(cur_index, modnum, x_size, y_size, sym_type);
                            // if result.is_none() {
                            //     println!("got an invalid result for index ({}, {}) with modnum = {}", &i, &j, &modnum);
                            // }
                            if let Some(valid_index) = result {
                                neighbors.push(valid_index);
                            }
                        }
                        if genval >= spin_up_chance {
                            give_map.push(SpinNode::cons_node(
                                spin_unit,
                                array![i as f64, j as f64],
                                *cur_coord,
                                RwLock::new(neighbors),
                            ));
                        } else {
                            give_map.push(SpinNode::cons_node(
                                -spin_unit,
                                array![i as f64, j as f64],
                                *cur_coord,
                                RwLock::new(neighbors),
                            ));
                        }
                    }
                }
            }
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
                let b1 = array![cur_basis[(0, 0)], cur_basis[(0, 1)]];
                let b2 = array![cur_basis[(1, 0)], cur_basis[(1, 1)]];
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
                            } else if (j % 2 == 1) && ((cur_index + y_size + 1) <= x_size * y_size)
                            {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index + y_size);
                            }
                        } else {
                            // when we are on an odd increment the +x direction
                            if (j % 2 == 0) && ((cur_index + y_size + 1) <= x_size * y_size) {
                                // if the y index is even, then its neighbor is in the -x direction
                                neighbors.push(cur_index + y_size);
                            } else if (j % 2 == 1) && (cur_index >= y_size) {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index - y_size);
                            }
                        }
                        // randomization conditions would go here or somewhere near here.
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        if genval >= spin_up_chance {
                            give_map.push(SpinNode::cons_node(
                                spin_unit,
                                array![i as f64, j as f64],
                                *cur_coord.clone(),
                                RwLock::new(neighbors),
                            ));
                        } else {
                            give_map.push(SpinNode::cons_node(
                                -spin_unit,
                                array![i as f64, j as f64],
                                *cur_coord.clone(),
                                RwLock::new(neighbors),
                            ));
                        }
                        if i % 2 == 0 {
                            // when we are on an even increment of b2
                            if j % 2 == 0 {
                                // for even increments of the y vector
                                *cur_coord = *cur_coord.clone() + &b2;
                            } else {
                                *cur_coord = *cur_coord.clone() + &b3;
                            }
                        } else if j % 2 == 1 {
                            // for odd increments of the y vector
                            *cur_coord = *cur_coord.clone() + &b2;
                        } else {
                            *cur_coord = *cur_coord.clone() + &b3;
                        }
                    }
                    *cur_coord = array![b1[0] + 2. * &b2[0], 0.];
                }
            }
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
            }
        }
        println!("Done!");

        // let psum_signaler = SignalContainer::new(cpu_num);
        let energy_signlaer = SignalContainer::new(cpu_num);
        let cluster_go_signaler = SignalContainer::new(cpu_num);
        let cluster_done_signaler = SignalContainer::new(cpu_num);
        let cluster_queue_signaler = SignalContainer::new(0);
        let cluster_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cpu_num)
            .build()
            .unwrap();
        let energy_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cpu_num)
            .build()
            .unwrap();

        Self {
            internal_lattice: Lattice::new(give_map),
            basis: RwLock::new(inner_basis),
            spin_unit,
            x_size,
            y_size,
            total_nodes: x_size * y_size,
            sym_type,
            exchange_constant,
            num_threads: cpu_num,
            div,
            rem,
            cluster_threads_started: false,
            energy_threads_started: false,
            energy_pool,
            cluster_pool,
            energy_psum_signaler: energy_signlaer,
            cluster_go_stop_signaler: cluster_go_signaler,
            cluster_done_signaler,
            cluster_queue_signaler,
            flip_index_vector: Arc::new(RwLock::new(vec![])),
            touched_index_vec: Arc::new(RwLock::new(vec![])),
            fname,
        }
    }

    pub fn save_state(&self, fname: &str) {
        self.internal_lattice.export_state_to_file(fname);
    }

    pub fn load_state(&mut self, fname: &str) {
        self.internal_lattice
            .load_state_from_file(fname.to_string());
    }

    pub fn get_energy(&mut self) -> f64 {
        if !self.energy_threads_started {
            let n_jobs = self.num_threads;
            for i in 0..n_jobs {
                let this_thread_num = i;
                let shared_internal_vector = self.internal_lattice.internal_vector.clone();
                let energy_psum_signaler = self.energy_psum_signaler.clone();

                // start threads here
                let div = self.div;
                let rem = self.rem;

                self.energy_pool.spawn(move || {
                    let mut energy_psum = 0.;
                    unsafe {
                        let range = if this_thread_num == n_jobs - 1 {
                            (this_thread_num * div)..((this_thread_num + 1) * div + rem)
                        } else {
                            (this_thread_num * div)..((this_thread_num + 1) * div)
                        };
                        if let Ok(read_lock_internal_vector) = shared_internal_vector.read() {
                            for cur_node_index in range {
                                // iterate through all nodes in this threads range
                                let cur_node =
                                    read_lock_internal_vector.get_unchecked(cur_node_index);
                                let read_lock_neighbors_of_cur_node = cur_node
                                    .neighbors
                                    .read()
                                    .expect("Couldnt read from internal_vector!");
                                let mut sum_of_cur_node_nbrs_spins = 0.;
                                for nbr_index in read_lock_neighbors_of_cur_node.iter() {
                                    // if cur_node_index == 1 {
                                    //     // print the neighbors indices of the node
                                    //     println!(
                                    //         "index {} has neighbor {}",
                                    //         &cur_node_index, *nbr_index
                                    //     );
                                    // }
                                    // println!("node {} has {} neighbors", &cur_node_index, &read_lock_neighbors_of_cur_node.len());
                                    // iterate over the current nodes neighbors
                                    let cur_nodes_neighbor_ref =
                                        read_lock_internal_vector.get_unchecked(*nbr_index);
                                    sum_of_cur_node_nbrs_spins += cur_nodes_neighbor_ref.get_spin();
                                }
                                // println!(
                                //     "sum of spins for node {} with spin of {} is {}, giving an energy of {}",
                                //     &cur_node_index, &cur_node.get_spin(), &sum_of_cur_node_nbrs_spins,
                                //     sum_of_cur_node_nbrs_spins * cur_node.get_spin()
                                // );
                                // multiply the sum of spins by the current spin value
                                energy_psum += sum_of_cur_node_nbrs_spins * cur_node.get_spin();
                            }
                        }
                    }
                    energy_psum_signaler.send(energy_psum).unwrap();
                });
            } // for i in 0..n_jobs
        }
        // self.energy_threads_started = true;
        let mut energy = 0.;
        for _ in 0..self.num_threads {
            if let Ok(psum) = self.energy_psum_signaler.recv() {
                energy += psum;
            } else {
                panic!("Channel is already closed, something went wrong!")
            }
        }
        energy
    }

    pub fn get_magnitization(&self) -> f64 {
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();
        read_lock_internal_vector
            .par_iter()
            .map(|i| i.get_spin())
            .sum()
    }

    fn start_cluster_threads(&mut self) {
        if !self.cluster_threads_started {
            let n_jobs = self.num_threads;
            for i in 0..n_jobs {
                println!("starting cluster thread {i}");
                let shared_flip_vec = self.flip_index_vector.clone();
                let shared_touched_vec = self.touched_index_vec.clone();
                let shared_data = self.internal_lattice.internal_vector.clone();
                let cluster_go_stop_signaler = self.cluster_go_stop_signaler.clone();
                let cluster_queue_signaler = self.cluster_queue_signaler.clone();
                let cluster_done_signaler = self.cluster_done_signaler.clone();

                self.cluster_pool.spawn(move || {
                    let mut rng = thread_rng();
                    'outer: loop {
                        if let Ok((sent_signal, balance_condition, original_target_spin)) =
                            cluster_go_stop_signaler.recv()
                        {
                            if SignalType::SigStart == sent_signal {
                                'inner: loop {
                                    if let (
                                        Ok(read_lock_internal_vector),
                                        Ok((proposed_nodes_index, proposed_nodes_spin)),
                                    ) = (
                                        shared_data.read().as_deref(),
                                        cluster_queue_signaler.try_recv(),
                                    ) {
                                        unsafe {
                                            let proposed_node = read_lock_internal_vector
                                                .get_unchecked(proposed_nodes_index);
                                            if proposed_node.get_status() != StateValue::Unmarked
                                                || proposed_nodes_spin != original_target_spin
                                            {
                                                proposed_node.marked.write().unwrap().set_marked();
                                                shared_touched_vec.write().unwrap().push(proposed_nodes_index);
                                                continue 'inner;
                                            }
                                            // indicate that a thread is processing this node
                                            if let Ok(mut write_lock_proposed_node) =
                                                proposed_node.marked.write()
                                            {
                                                write_lock_proposed_node.set_processing();
                                                shared_touched_vec.write().unwrap().push(proposed_nodes_index);
                                            } else {
                                                continue 'inner;
                                            }
                                            let num_neighbors =
                                                proposed_node.neighbors.read().unwrap().len();
                                            'nbr_cycle: for i in 0..num_neighbors {
                                                let nbr_index = *proposed_node
                                                    .neighbors
                                                    .read()
                                                    .unwrap()
                                                    .get_unchecked(i);
                                                let nbr_spin = read_lock_internal_vector
                                                    .get_unchecked(nbr_index)
                                                    .get_spin();
                                                if nbr_spin != original_target_spin {
                                                    continue 'nbr_cycle;
                                                }
                                                if rng.gen_range(0_f64..1_f64) < balance_condition {
                                                    // get write lock to accepted node
                                                    if let Ok(mut write_lock) =
                                                        shared_flip_vec.try_write()
                                                    {
                                                        if read_lock_internal_vector
                                                            .get_unchecked(nbr_index)
                                                            .get_status()
                                                            != StateValue::Pushed
                                                        {
                                                            let mut write_lock_nbr =
                                                                read_lock_internal_vector
                                                                    .get_unchecked(nbr_index)
                                                                    .marked
                                                                    .write()
                                                                    .unwrap();
                                                            write_lock_nbr.set_pushed();
                                                            drop(write_lock_nbr);
                                                            write_lock.push(nbr_index);
                                                            cluster_queue_signaler
                                                                .send((
                                                                    nbr_index,
                                                                    original_target_spin,
                                                                ))
                                                                .unwrap();
                                                        }
                                                    } else {
                                                        continue 'nbr_cycle;
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        // send done signal to main
                                        cluster_done_signaler.send(1).unwrap();
                                        continue 'outer;
                                    }
                                }
                            } else if SignalType::SigStop == sent_signal {
                                println!("Stop stignal received, closing thread.");
                                return;
                            } else {
                                panic!("An unexpected event occured!")
                            }
                        }
                    }
                });
            } // for i in 0..n_jobs
        }
        self.cluster_threads_started = true;
    }

    pub fn stop_threads(self) {
        if self.cluster_threads_started {
            self.cluster_go_stop_signaler.close_channel();
        }
        println!("Threads closed!\n")
    }

    pub fn get_spin_at(&self, index: usize) -> f64 {
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            match value.get(index) {
                Some(node) => node.get_spin(),
                None => 0.,
            }
        } else {
            panic!("Couldnt get spin")
        }
    }

    /// Flips the spin pointed to by the parameter `index`. This function requires a write lock on `internal_lattice.internal_vector`
    pub fn flip_node_at(&self, index: usize) {
        if let Ok(mut lock) = self.internal_lattice.internal_vector.write() {
            match lock.get_mut(index) {
                Some(node) => node.flip_spin(),
                None => (),
            }
        }
    }

    fn save_files(&self, m_vec: Vec<f64>, e_vec: Vec<f64>, c_vec: Vec<f64>, x_vec: Vec<f64>) {
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
            }
        };
        // start writing data
        for item in m_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    file_progress.inc(1);
                    continue;
                }
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
            }
        };
        for item in e_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    file_progress.inc(1);
                    continue;
                }
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
            }
        };
        for item in x_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    file_progress.inc(1);
                    continue;
                }
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
            }
        };
        for item in c_vec.into_iter() {
            let write_str = format!("{}\n", item);
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    file_progress.inc(1);
                    continue;
                }
            }
        }
        file_progress.abandon();
        println!("Sucessfully wrote files!");
    }

    pub fn spin_energy(&mut self, beta_list: Vec<f64>, times: usize, iteration_scheme: usize, ignore_n_runs: usize) {
        // if self.magnitiztion_threads_started != true { self.start_magnitization_threads(); }
        // if self.energy_threads_started != true { self.start_energy_threads(); }
        assert!(ignore_n_runs < times);
        if iteration_scheme == 1 {
            self.start_cluster_threads();
        }

        let (temp_magnitization, initial_energy) = (self.get_magnitization(), self.get_energy());
        let initial_magnitization = temp_magnitization.abs();
        // println!("\nThe initial energy is {}, and the initial magnitization is {}.\n", initial_energy, initial_magnitization);

        let mut energy: f64;
        let mut magnitization: f64;
        let mut ms: Vec<f64>;
        let mut es: Vec<f64>;

        let mut m_vec: Vec<f64> = vec![];
        let mut e_vec: Vec<f64> = vec![];
        let mut c_vec: Vec<f64> = vec![];
        let mut x_vec: Vec<f64> = vec![];
        let n1 = (self.total_nodes as f64) * (times as f64);
        let n2 = (self.total_nodes as f64) * (times as f64).powf(2.);

        let len_beta = beta_list.len();
        let tot_time: u64 = times as u64 * len_beta as u64;

        let bar1 = ProgressBar::new(tot_time);

        // load the initial state at tbe begining of each new beta
        for beta_val in beta_list {
            self.internal_lattice
                .load_state_from_file(self.fname.to_string());

            magnitization = initial_magnitization;
            energy = initial_energy;
            let (mut d_energy, mut d_mag);
            // let mut cur_mag = 0.;
            // let mut cur_energy = 0.;
            ms = vec![0., 0.];
            es = vec![0., 0.];

            for cur_t in 0..times {
                // preform an iteration scheme of metropolis or wolff
                if iteration_scheme == 0 {
                    (d_energy, d_mag) = self.metropolis_iter(&beta_val);
                    if cur_t > ignore_n_runs {
                        magnitization = self.get_magnitization();
                        energy += d_energy;
                    }
                } else if iteration_scheme == 1 {
                    (d_energy, d_mag) = self.wolff_iter(&beta_val);
                    if cur_t > ignore_n_runs {
                        magnitization += d_mag;
                        energy += d_energy;
                    }
                } else {
                    panic!("invalid option! expected 0 or 1, got {}", iteration_scheme);
                }
                // Beginning the data collection after a few iterations gives better
                // overall graphs becuase the system has had some time to relax and the
                // resulting data doesnt wildly fluctuate at the start of the resulting
                // data plot.
                if cur_t > ignore_n_runs {
                    ms[0] += magnitization.abs();
                    ms[1] += magnitization.powf(2.);
                    es[0] += energy;
                    es[1] += energy.powf(2.);
                }
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
            self.save_files(m_vec, e_vec, c_vec, x_vec);
        }
    }

    /// Purpose
    /// -------
    /// Mutates the `self.LinkedLat` lattice of spins by one Iteration of the Wolff Algorithm.
    fn wolff_iter(&mut self, beta: &f64) -> (f64, f64) {
        // self.unmark_all_nodes();
        let balance_condition = 1. - consts::E.powf(-2. * beta * self.exchange_constant);
        let mut rng_spin_select = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;
        let mut delta_energy: f64 = 0.;
        let mut delta_mag: f64 = 0.;

        // select a random node
        'node_selection: loop {
            target_index = rng_spin_select.gen_range(0..(self.x_size * self.y_size - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }

        // push the random node neighbors indicies to the work queue
        if let Ok(read_lock_internal_vector) = self.internal_lattice.internal_vector.read() {
            if let Some(target_node) = read_lock_internal_vector.get(target_index) {
                let read_lock_nbrs = target_node.neighbors.read().unwrap();
                for i in 0..read_lock_nbrs.len() {
                    // SAFETY: bounds checked in Driver::new, garaunteed valid return.
                    unsafe {
                        let nbr_index = *read_lock_nbrs.get_unchecked(i);
                        let nbr_spin = self.get_spin_at(nbr_index);
                        if nbr_spin == target_spin {
                            // if the spin is the same as the randomly picked node, add it to the
                            // queue.
                            self.cluster_queue_signaler
                                .send((nbr_index, nbr_spin))
                                .expect("Couldnt send message through cluster_queue_signaler!");
                        }
                    }
                }
            }
        }

        // spin up threads for generating a cluster flip
        for _ in 0..self.num_threads {
            self.cluster_go_stop_signaler
                .send((SignalType::SigStart, balance_condition, target_spin))
                .unwrap();
        }
        // waiting for threads to finish
        while !self.cluster_done_signaler.is_full() {
            stdth::sleep(time::Duration::from_micros(1));
        }
        // Clear queue
        for _ in 0..self.num_threads {
            // clear the queue
            _ = self.cluster_done_signaler.recv().unwrap();
        }
        // now we flip the resulting cluster and calculate the change in energy and magnitization
        delta_mag += -1. * target_spin * self.flip_index_vector.read().unwrap().len() as f64;
        // get a write lock to flip_index_vecor
        unsafe {
            // println!("indices of nodes to be flipped is/are: ");
            // for item in write_lock_flip_index_vector.iter() {
                //     println!("{}", item);
                // }
                // println!("\n");
                // get_input_as_usize(None);
            'flip: loop {
                let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap();
                let opt_cur_flip_node_index= write_lock_flip_index_vector.pop();
                if let Some(cur_flip_node_index) = opt_cur_flip_node_index {
                    let mut write_lock_internal_vector =
                    self.internal_lattice.internal_vector.write().unwrap();
                    // preform the accepted spin flip
                    write_lock_internal_vector
                    .get_unchecked_mut(cur_flip_node_index)
                    .flip_spin();
                } else { break 'flip; }
                drop(write_lock_flip_index_vector);
                // flip already happened, so take the negative of the change in energy
                delta_energy -= self.calculate_energy_change_of_suggested_flip(opt_cur_flip_node_index.unwrap());
            }
            // unmark all touched nodes, I cant really think of a better way to do this at the moment.
            let mut write_lock_touched_index_vector = self.touched_index_vec.write().unwrap();
            if let Some(touched_index) = write_lock_touched_index_vector.pop() {
                let mut write_lock_internal_vector =
                    self.internal_lattice.internal_vector.write().unwrap();
                write_lock_internal_vector
                    .get_unchecked_mut(touched_index)
                    .marked
                        .write()
                        .unwrap()
                        .set_unmark();
            }
        }

        return (delta_energy, delta_mag)
    }

    /// Evolves the lattice by one iteration using the metropolis-hastings scheme.
    fn metropolis_iter(&mut self, beta: &f64) -> (f64, f64) {
        let mut rngspin = thread_rng();
        let mut rng_flip = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;

        // select a random node
        'node_selection: loop {
            target_index = rngspin.gen_range(0..(self.x_size * self.y_size - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }

        // calculate dE of the proposed spin flip at x_y
        let mut delta_energy: f64 = self.calculate_energy_change_of_suggested_flip(target_index);
        let delta_mag: f64;

        // flip node if accepted
        if (delta_energy > 0.
            && (rng_flip.gen_range(0_f64..1_f64)
                < consts::E.powf(-beta * self.exchange_constant * delta_energy)))
            || delta_energy <= 0.
        {
            self.flip_node_at(target_index);
            delta_mag = -target_spin;
        } else {
            delta_energy = 0.;
            delta_mag = 0.;
        }

        // println!("using get_energy={}, using calculate_energy_change_of_suggested_flip={}", delta_energy, dE);

        (delta_energy, delta_mag)
    }

    fn get_neighbors_spin_sum(&self, x_y: usize) -> f64 {
        let mut nbr_energy = 0.;
        if let Ok(value) = self.internal_lattice.internal_vector.read() {
            if let Some(node) = value.get(x_y) {
                // let target_spin = node.get_spin();
                let read_lock_nbrs = node.neighbors.read().unwrap();
                for ith_nbr_index in 0..read_lock_nbrs.len() {
                    // SAFETY: bounds checked in random node selection, garaunteed valid return
                    unsafe {
                        nbr_energy +=
                            self.get_spin_at(*read_lock_nbrs.get_unchecked(ith_nbr_index));
                    }
                }
            }
        }
        nbr_energy
    }

    fn unmark_all_nodes(&mut self) {
        let mut write_lock = self.internal_lattice.internal_vector.write().unwrap();
        for node in write_lock.iter_mut() {
            node.marked.write().unwrap().set_unmark();
        }
    }

    pub fn calculate_energy_change_of_suggested_flip(&self, target_index: usize) -> f64 {
        let mut energy_i = 0.;
        let mut energy_f = 0.;
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();
        // SAFETY: x_y was checked already to be a valid node
        unsafe {
            let cur_target = read_lock_internal_vector.get_unchecked(target_index);
            let cur_target_nbrs_lock = cur_target.neighbors.read().unwrap();
            // iterate through the neighbors of the suggested node to flip
            for nbrs_index_of_flip_node in cur_target_nbrs_lock.iter() {
                let target_of_cur_target = read_lock_internal_vector.get_unchecked(*nbrs_index_of_flip_node);
                let target_of_cur_target_nbrs_lock = target_of_cur_target.neighbors.read().unwrap();
                // cycle through the neighbors of the neighbors of the suggested node to flip (wordy, yeah)
                for nbrs_of_nbrs_index_of_flip_node in target_of_cur_target_nbrs_lock.iter() {
                    let nbr_of_nbrs_of_flip_node = read_lock_internal_vector.get_unchecked(*nbrs_of_nbrs_index_of_flip_node);
                    // if nbrs_of_nbrs_index_of_flip_node is the same as the target_index, calculate the E_i and E_f
                    if *nbrs_of_nbrs_index_of_flip_node == target_index {
                        // get the energy using a flipped value of spin for nbr_of_nbrs_of_flip_node
                        energy_f += -nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                    }
                    else {
                        energy_f += nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                    }
                    // get the regular energy
                    energy_i += nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                }
            }   
        }
        return (energy_f - energy_i) / self.spin_unit
    }
}
