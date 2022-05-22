use ordered_float::OrderedFloat;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
//use std::time::Duration;
use crate::lat_node::*;
use crate::lattice_structure::Lattice;
use crate::signal_container::*;
use crate::{dividend_remainder, sim_params};
use indicatif::ProgressBar;
use ndarray::prelude::*;
use rand::thread_rng;
use rand::Rng;
use rayon::prelude::*;
use std::f64::consts;
use std::sync::{Arc, RwLock};
use std::{thread as stdth, time};

fn mean(data: &[f64]) -> Option<f64> {
    let sum = data.par_iter().sum::<f64>() as f64;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f64),
        _ => None,
    }
}

fn std_deviation(data: &[f64]) -> Option<f64> {
    match (mean(data), data.len()) {
        (Some(data_mean), count) if count > 0 => {
            let variance = data
                .par_iter()
                .map(|value| {
                    let diff = data_mean - (*value as f64);

                    diff * diff
                })
                .sum::<f64>()
                / count as f64;

            Some(variance.sqrt())
        }
        _ => None,
    }
}

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
                    return Some(index + x_size);
                } else {
                    return None;
                }
            } else if modnum == 2 {
                // -x
                if i > 0 {
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
    parameters: sim_params::SimulationParameters,
    internal_lattice: Lattice,
    num_threads: usize,
    div: usize,
    rem: usize,
    cluster_threads_started: bool,
    energy_threads_started: bool,
    energy_pool: rayon::ThreadPool,
    cluster_pool: rayon::ThreadPool,
    energy_psum_signaler: SignalContainer<f64>,
    cluster_go_stop_signaler: SignalContainer<(SignalType, f64, f64)>,
    cluster_done_signaler: SignalContainer<(f64, f64)>,
    cluster_queue_signaler: SignalContainer<usize>,
    touched_index_vec: Arc<RwLock<Vec<usize>>>,
}

impl Driver {
    pub fn new(parameters: sim_params::SimulationParameters) -> Self {
        let num_threads = num_cpus::get();
        let (div, rem) = dividend_remainder(parameters.num_nodes(), num_threads);
        let mut cur_coord: Box<Array1<f64>>;
        let mut rng = thread_rng();
        let mut give_map: Vec<SpinNode> = vec![];
        let s1 = s![0_usize, ..];
        let s2 = s![1_usize, ..];
        match parameters.get_symtype() {
            // match the generation pattern to the type passed to sym_type
            SymmetryType::C4V => {
                // C4V generates every node including the starting node, incrementing by 1 x or y unit vector each
                // itteraton.
                for j in 0..parameters.get_ysize() {
                    for i in 0..parameters.get_xsize() {
                        let cur_index = i + j * parameters.get_xsize();
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        // randomization conditions would go here or somewhere near here.
                        cur_coord = Box::new(
                            (i as f64) * &parameters.get_basis().slice(&s1)
                                + (j as f64) * &parameters.get_basis().slice(&s2),
                        );
                        // construct neighbors vector
                        let mut neighbors = vec![];
                        for modnum in 0..4 {
                            let result = indexmod(
                                cur_index,
                                modnum,
                                parameters.get_xsize(),
                                parameters.get_ysize(),
                                parameters.get_symtype(),
                            );

                            if let Some(valid_index) = result {
                                neighbors.push(valid_index);
                            }
                        }
                        if genval >= parameters.get_spin_up_chance() {
                            give_map.push(SpinNode::cons_node(
                                parameters.get_spin_unit(),
                                array![i as f64, j as f64],
                                *cur_coord,
                                RwLock::new(neighbors),
                            ));
                        } else {
                            give_map.push(SpinNode::cons_node(
                                -parameters.get_spin_unit(),
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

                let cur_basis = parameters.get_basis();
                let b1 = array![cur_basis[(0, 0)], cur_basis[(0, 1)]];
                let b2 = array![cur_basis[(1, 0)], cur_basis[(1, 1)]];
                let b3 = array![-b2[0], b2[1]];
                let mut cur_coord: Box<Array1<f64>> = Box::new(array![0., 0.]);

                for i in 0..parameters.get_xsize() {
                    for j in 0..parameters.get_ysize() {
                        let mut neighbors = vec![];
                        let cur_index = j + i * parameters.get_ysize();
                        // add the + and - y neighbors to the neighbors vector.
                        // this needs to be checked for all nodes.
                        if j > 0 {
                            // if the y index is greater than 0, we can add the -y node
                            neighbors.push(cur_index - 1);
                        }
                        if (j + 2) < parameters.get_ysize() {
                            // if the y index is 2 less than the
                            neighbors.push(cur_index + 1);
                        }
                        if i % 2 == 0 {
                            // when we are on an even increment the +x direction
                            if (j % 2 == 0) && (cur_index >= parameters.get_ysize()) {
                                // if the y index is even, then its neighbor is in the -x direction
                                neighbors.push(cur_index - parameters.get_ysize());
                            } else if (j % 2 == 1)
                                && ((cur_index + parameters.get_ysize() + 1)
                                    <= parameters.num_nodes())
                            {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index + parameters.get_ysize());
                            }
                        } else {
                            // when we are on an odd increment the +x direction
                            if (j % 2 == 0)
                                && ((cur_index + parameters.get_ysize() + 1)
                                    <= parameters.num_nodes())
                            {
                                // if the y index is even, then its neighbor is in the -x direction
                                neighbors.push(cur_index + parameters.get_ysize());
                            } else if (j % 2 == 1) && (cur_index >= parameters.get_ysize()) {
                                // if the y index is odd, then its neighbor is in the -x direction
                                neighbors.push(cur_index - parameters.get_ysize());
                            }
                        }
                        // randomization conditions would go here or somewhere near here.
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
                        if genval >= parameters.get_spin_up_chance() {
                            give_map.push(SpinNode::cons_node(
                                parameters.get_spin_unit(),
                                array![i as f64, j as f64],
                                *cur_coord.clone(),
                                RwLock::new(neighbors),
                            ));
                        } else {
                            give_map.push(SpinNode::cons_node(
                                -parameters.get_spin_unit(),
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
                todo!();
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
        let energy_signlaer = SignalContainer::new(num_threads);
        let cluster_go_signaler = SignalContainer::new(num_threads);
        let cluster_done_signaler = SignalContainer::new(num_threads);
        let cluster_queue_signaler = SignalContainer::new(0);
        let cluster_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();
        let energy_pool = rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap();

        Self {
            parameters,
            internal_lattice: Lattice::new(give_map),
            num_threads,
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
            touched_index_vec: Arc::new(RwLock::new(vec![])),
        }
    }

    pub fn get_parameters(&self) -> sim_params::SimulationParameters {
        self.parameters.clone()
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
                let ext_mag_field = self.parameters.get_b_field();

                self.energy_pool.spawn(move || {
                    let mut energy_psum = 0.;
                    let mut mag_field_psum = 0.;
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
                                if ext_mag_field != 0. {
                                    mag_field_psum -= ext_mag_field * cur_node.get_spin();
                                }
                                for nbr_index in read_lock_neighbors_of_cur_node.iter() {
                                    // iterate over the current nodes neighbors
                                    let cur_nodes_neighbor_ref =
                                        read_lock_internal_vector.get_unchecked(*nbr_index);
                                    sum_of_cur_node_nbrs_spins += cur_nodes_neighbor_ref.get_spin();
                                }
                                // multiply the sum of spins by the current spin value
                                energy_psum -= sum_of_cur_node_nbrs_spins * cur_node.get_spin();
                            }
                        }
                    }
                    energy_psum_signaler
                        .send(energy_psum + mag_field_psum)
                        .unwrap();
                });
            } // for i in 0..n_jobs
        }
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
                let shared_touched_vec = self.touched_index_vec.clone();
                let shared_data = self.internal_lattice.internal_vector.clone();
                let cluster_go_stop_signaler = self.cluster_go_stop_signaler.clone();
                let cluster_queue_signaler = self.cluster_queue_signaler.clone();
                let cluster_done_signaler = self.cluster_done_signaler.clone();
                let ext_mag_field = self.parameters.get_b_field();
                let spin_unit = self.parameters.get_spin_unit();

                self.cluster_pool.spawn(move || {
                    let mut rng = thread_rng();
                    unsafe {
                        'main_loop: loop {
                            if let Ok(start_signal) = cluster_go_stop_signaler.recv() {
                                if start_signal.0 == SignalType::SigStop {
                                    break 'main_loop;
                                }
                                let mut delta_mag = 0.;
                                let mut delta_energy = 0.;
                                let (balance_condition, target_spin) =
                                    (start_signal.1, start_signal.2);
                                'cluster_loop: loop {
                                    if let Ok(node_index) = cluster_queue_signaler.try_recv() {
                                        let read_lock_internal_vec = shared_data.read().unwrap();
                                        if read_lock_internal_vec
                                            .get_unchecked(node_index)
                                            .get_status()
                                            != StateValue::Unmarked
                                        {
                                            continue 'cluster_loop;
                                        } else {
                                            let mut write_lock_touched_index_vector =
                                                shared_touched_vec.write().unwrap();
                                            write_lock_touched_index_vector.push(node_index);
                                            drop(write_lock_touched_index_vector);
                                            read_lock_internal_vec
                                                .get_unchecked(node_index)
                                                .marked
                                                .write()
                                                .unwrap()
                                                .set_marked();
                                            if target_spin
                                                != read_lock_internal_vec
                                                    .get_unchecked(node_index)
                                                    .get_spin()
                                            {
                                                continue 'cluster_loop;
                                            }
                                        }
                                        drop(read_lock_internal_vec);
                                        if rng.gen_range(0_f64..1_f64) < balance_condition {
                                            let mut energy_i = 0.;
                                            let mut energy_f = 0.;
                                            let mut mag_energy_i = 0.;
                                            let mut mag_energy_f = 0.;
                                            // iterate through the neighbors of the suggested node to flip
                                            for nbrs_index_of_flip_node in shared_data
                                                .read()
                                                .unwrap()
                                                .get_unchecked(node_index)
                                                .neighbors
                                                .read()
                                                .unwrap()
                                                .iter()
                                            {
                                                let read_lock_internal_vec =
                                                    shared_data.read().unwrap();
                                                let target_of_cur_target = read_lock_internal_vec
                                                    .get_unchecked(*nbrs_index_of_flip_node);
                                                let target_of_cur_target_nbrs_lock =
                                                    target_of_cur_target.neighbors.read().unwrap();
                                                // cycle through the neighbors of the neighbors of the suggested node to flip (wordy, yeah)
                                                for nbrs_of_nbrs_index_of_flip_node in
                                                    target_of_cur_target_nbrs_lock.iter()
                                                {
                                                    let read_lock_internal_vec =
                                                        shared_data.read().unwrap();
                                                    let nbr_of_nbrs_of_flip_node =
                                                        read_lock_internal_vec.get_unchecked(
                                                            *nbrs_of_nbrs_index_of_flip_node,
                                                        );
                                                    // if nbrs_of_nbrs_index_of_flip_node is the same as the target_index, calculate the E_i and E_f
                                                    if *nbrs_of_nbrs_index_of_flip_node
                                                        == node_index
                                                    {
                                                        // get the energy using a flipped value of spin for nbr_of_nbrs_of_flip_node
                                                        energy_f -= -nbr_of_nbrs_of_flip_node
                                                            .get_spin()
                                                            * target_of_cur_target.get_spin();
                                                        mag_energy_f -= -ext_mag_field
                                                            * target_of_cur_target.get_spin();
                                                    } else {
                                                        energy_f -= nbr_of_nbrs_of_flip_node
                                                            .get_spin()
                                                            * target_of_cur_target.get_spin();
                                                        mag_energy_i -= ext_mag_field
                                                            * target_of_cur_target.get_spin();
                                                    }
                                                    // get the regular energy
                                                    energy_i -= nbr_of_nbrs_of_flip_node.get_spin()
                                                        * target_of_cur_target.get_spin();
                                                }
                                            }
                                            delta_energy += (energy_f - energy_i) / spin_unit
                                                + (mag_energy_f - mag_energy_i);
                                            delta_mag -= target_spin;
                                            let mut write_lock_internal_vector =
                                                shared_data.write().unwrap();
                                            // println!("flipped a spin");
                                            write_lock_internal_vector
                                                .get_unchecked_mut(node_index)
                                                .flip_spin();
                                            drop(write_lock_internal_vector);
                                            let read_lock_internal_vec =
                                                shared_data.read().unwrap();
                                            let nbs_lock_of_target = read_lock_internal_vec
                                                .get_unchecked(node_index)
                                                .neighbors
                                                .read()
                                                .unwrap();
                                            for nbr in nbs_lock_of_target.iter() {
                                                let read_lock_nbr =
                                                    read_lock_internal_vec.get_unchecked(*nbr);
                                                if read_lock_nbr.get_spin() == target_spin {
                                                    read_lock_internal_vec
                                                        .get_unchecked(*nbr)
                                                        .marked
                                                        .write()
                                                        .unwrap()
                                                        .set_pushed();
                                                    cluster_queue_signaler.send(*nbr).unwrap();
                                                }
                                            }
                                        }
                                    } else {
                                        cluster_done_signaler
                                            .send((delta_mag, delta_energy))
                                            .unwrap();
                                        continue 'main_loop;
                                    } // if the stack is empty, wait at the start of the main loop for the next go signal
                                }
                            } else {
                                break 'main_loop;
                            } // break the loop and exit the thread of the channel is closed
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

    fn save_files(
        &self,
        m_vec: Vec<f64>,
        e_vec: Vec<f64>,
        c_vec: Vec<f64>,
        x_vec: Vec<f64>,
        beta_vec: Vec<f64>,
    ) {
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

        // setup saving the output file
        let path = Path::new("beta.dat");
        let display = path.display();
        // Write output data to a file
        let mut file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            }
        };
        for item in beta_vec.into_iter() {
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

    pub fn spin_energy(
        &mut self,
        beta_vec: Vec<f64>,
        times: usize,
        iteration_scheme: usize,
        ignore_n_runs: usize,
        anneal: bool,
    ) {
        if !anneal {
            assert!(ignore_n_runs < times);
        }
        if iteration_scheme == 1 {
            self.start_cluster_threads();
        }

        let (temp_magnitization, initial_energy) = (self.get_magnitization(), self.get_energy());
        let initial_magnitization = temp_magnitization.abs();
        println!(
            "\nThe initial energy is {}, and the initial magnitization is {}.\n",
            initial_energy, initial_magnitization
        );

        let mut m_vec: Vec<f64> = vec![];
        let mut e_vec: Vec<f64> = vec![];
        let mut c_vec: Vec<f64> = vec![];
        let mut x_vec: Vec<f64> = vec![];
        let n = self.parameters.num_nodes() as f64;

        let len_beta = beta_vec.len();
        let tot_time: u64 = times as u64 * len_beta as u64;

        let bar1 = ProgressBar::new(tot_time);
        // load the initial state at tbe begining of each new beta
        for beta_val in beta_vec.iter() {
            self.internal_lattice
                .load_state_from_file(self.parameters.get_fname().to_string());

            let mut magnitization: f64 = initial_magnitization;
            let mut energy: f64 = initial_energy;
            let (mut d_energy, mut d_mag);
            let mut mag_vec = vec![];
            let mut energy_vec = vec![];

            for cur_t in 0..times {
                // preform an iteration scheme of metropolis or wolff
                if iteration_scheme == 0 {
                    (d_energy, d_mag) = self.metropolis_iter(beta_val);
                } else if iteration_scheme == 1 {
                    (d_energy, d_mag) = self.wolff_iter_multithreaded(beta_val);
                } else {
                    panic!("invalid option! expected 0 or 1, got {}", iteration_scheme);
                }
                // Beginning the data collection after a few iterations gives better
                // overall graphs becuase the system has had some time to relax and the
                // resulting data doesnt wildly fluctuate at the start of the resulting
                // data plot.
                if cur_t + 1 > ignore_n_runs && !anneal {
                    // // preform a sum and sum of squares for statistics later
                    magnitization += d_mag;
                    energy += d_energy;
                    mag_vec.push(magnitization / n);
                    energy_vec.push(energy / n);
                }
                bar1.inc(1);
            }
            if !anneal {
                m_vec.push(mean(&mag_vec).unwrap());
                e_vec.push(mean(&energy_vec).unwrap());
                c_vec.push(std_deviation(&energy_vec).unwrap());
                x_vec.push(std_deviation(&mag_vec).unwrap());
            }
        } // for beta_val in beta_list end
        bar1.abandon();
        println!("spin_energy finished!");

        // save files
        if anneal != true {
            self.save_files(m_vec, e_vec, c_vec, x_vec, beta_vec);
        }
    }

    /// Purpose
    /// -------
    /// Mutates the `self.LinkedLat` lattice of spins by one Iteration of the Wolff Algorithm.
    fn wolff_iter_multithreaded(&mut self, beta: &f64) -> (f64, f64) {
        // self.unmark_all_nodes();
        let balance_condition = 1. - consts::E.powf(-2. * beta * self.parameters.get_j());
        let mut rng_spin_select = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;
        let mut delta_energy: f64 = 0.;
        let mut delta_mag: f64 = 0.;

        // select a random node
        'node_selection: loop {
            target_index = rng_spin_select
                .gen_range(0..(self.parameters.get_xsize() * self.parameters.get_ysize() - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }

        if rng_spin_select.gen_range(0_f64..1_f64) < balance_condition {
            delta_energy += self.calculate_energy_change_of_suggested_flip(target_index);
            delta_mag += -target_spin;
            self.flip_node_at(target_index);

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
                                    .send(nbr_index)
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
                // clear the queue and total changes
                let (tmp_delta_mag, tmp_delta_energy) = self.cluster_done_signaler.recv().unwrap();
                delta_mag += tmp_delta_mag;
                delta_energy += tmp_delta_energy;
            }

            unsafe {
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
                    drop(write_lock_internal_vector);
                    if self
                        .internal_lattice
                        .internal_vector
                        .read()
                        .unwrap()
                        .get_unchecked(touched_index)
                        .get_spin()
                        != -target_spin
                    {
                        delta_energy +=
                            self.calculate_energy_change_of_suggested_flip(target_index);
                        self.internal_lattice
                            .internal_vector
                            .write()
                            .unwrap()
                            .get_unchecked_mut(touched_index)
                            .flip_spin();
                    }
                }
            }
        }

        (delta_energy, delta_mag)
    }

    /// Evolves the lattice by one iteration using the metropolis-hastings scheme.
    fn metropolis_iter(&mut self, beta: &f64) -> (f64, f64) {
        let mut rngspin = thread_rng();
        let mut rng_flip = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;

        // select a random node
        'node_selection: loop {
            target_index = rngspin
                .gen_range(0..(self.parameters.get_xsize() * self.parameters.get_ysize() - 1));
            target_spin = self.get_spin_at(target_index);
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }
        // get the change in energy of the suggested spin flip
        let mut delta_energy: f64 = self.calculate_energy_change_of_suggested_flip(target_index);
        let balance_condition = consts::E.powf(-beta * self.parameters.get_j() * delta_energy);

        let delta_mag: f64;

        // flip node if accepted
        match OrderedFloat(delta_energy).cmp(&OrderedFloat(0_f64)) {
            std::cmp::Ordering::Less | std::cmp::Ordering::Equal => {
                self.flip_node_at(target_index);
                delta_mag = -target_spin;
            }
            std::cmp::Ordering::Greater => {
                if rng_flip.gen_range(0_f64..1_f64) < balance_condition {
                    self.flip_node_at(target_index);
                    delta_mag = -target_spin;
                } else {
                    delta_energy = 0.;
                    delta_mag = 0.;
                }
            }
        }

        (delta_energy, delta_mag)
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
        let mut mag_energy_i = 0.;
        let mut mag_energy_f = 0.;
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();

        // SAFETY: x_y was checked already to be a valid node
        unsafe {
            let cur_target = read_lock_internal_vector.get_unchecked(target_index);
            let cur_target_nbrs_lock = cur_target.neighbors.read().unwrap();
            // iterate through the neighbors of the suggested node to flip
            for nbrs_index_of_flip_node in cur_target_nbrs_lock.iter() {
                let target_of_cur_target =
                    read_lock_internal_vector.get_unchecked(*nbrs_index_of_flip_node);
                let target_of_cur_target_nbrs_lock = target_of_cur_target.neighbors.read().unwrap();
                // cycle through the neighbors of the neighbors of the suggested node to flip (wordy, yeah)
                for nbrs_of_nbrs_index_of_flip_node in target_of_cur_target_nbrs_lock.iter() {
                    let nbr_of_nbrs_of_flip_node =
                        read_lock_internal_vector.get_unchecked(*nbrs_of_nbrs_index_of_flip_node);
                    // if nbrs_of_nbrs_index_of_flip_node is the same as the target_index, calculate the E_i and E_f
                    if *nbrs_of_nbrs_index_of_flip_node == target_index {
                        // get the energy using a flipped value of spin for nbr_of_nbrs_of_flip_node
                        energy_f -=
                            -nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                        mag_energy_f -=
                            -self.parameters.get_b_field() * target_of_cur_target.get_spin();
                    } else {
                        energy_f -=
                            nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                        mag_energy_i -=
                            self.parameters.get_b_field() * target_of_cur_target.get_spin();
                    }
                    // get the regular energy
                    energy_i -=
                        nbr_of_nbrs_of_flip_node.get_spin() * target_of_cur_target.get_spin();
                }
            }
        }
        (energy_f - energy_i) / self.parameters.get_spin_unit() + (mag_energy_f - mag_energy_i)
    }
}
