use ordered_float::OrderedFloat;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
//use std::time::Duration;
use crate::{lat_node::*, thread_workers};
use crate::lattice_structure::Lattice;
use crate::signal_container::*;
use crate::{dividend_remainder, sim_params};
use indicatif::ProgressBar;
use ndarray::prelude::*;
use rand::thread_rng;
use rand::Rng;
use rayon::prelude::*;
use std::f64::consts;
use std::sync::RwLock;

fn mean(data: &[f64]) -> Option<f64> {
    let sum = data.par_iter().sum::<f64>() as f64;
    let count = data.len();

    match count {
        positive if positive > 0 => Some(sum / count as f64),
        _ => None,
    }
}

/// Purpose
/// -------
/// Multithreaded implementation of the sample standard deviation of a array/vector of float 64 numbers.
/// Uses the multihreaded implementation of `mean` as well.
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
                / ((count - 1) as f64);

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
                    // if the y index is 2 less than the max y heigth
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
    energy_psum_signaler: SignalContainer<f64>,
    cluster_done_signaler: SignalContainer<Array1<f64>>,
    cluster_queue_signaler: SignalContainer<usize>,
    // touched_index_vec: Arc<RwLock<Vec<usize>>>,
}

impl Driver {
    pub fn new(parameters: sim_params::SimulationParameters) -> Self {
    // faster this way :(, I still have a lot to learn in rust, but it does work either way at least...
        let num_threads = 1;//if num_cpus::get() > 1 {num_cpus::get() - 1} else {1};
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
                        let genval: f64 = rng.gen_range(0_f64..1_f64);
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
        let cluster_done_signaler: SignalContainer<ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>>> = SignalContainer::new(num_threads);
        let cluster_queue_signaler: SignalContainer<usize> = SignalContainer::new(0);

        Self {
            parameters,
            internal_lattice: Lattice::new(give_map),
            num_threads,
            div,
            rem,
            energy_psum_signaler: energy_signlaer,
            cluster_done_signaler,
            cluster_queue_signaler,
            // touched_index_vec: Arc::new(RwLock::new(vec![])),
        }
    }

    pub fn get_parameters(&self) -> sim_params::SimulationParameters {
        self.parameters.clone()
    }

    pub fn update_j(&mut self, j_val: f64) {
        self.parameters.set_j(j_val);
    }

    pub fn update_b(&mut self, b_val: f64) {
        self.parameters.set_b(b_val);
    }

    pub fn save_state(&self, fname: &str) {
        self.internal_lattice.export_state_to_file(fname);
    }

    pub fn load_state(&mut self, fname: &str) {
        self.internal_lattice
            .load_state_from_file(fname.to_string());
    }
    
    fn get_cluster(&mut self, target_spin: f64, balance_condition: f64) -> Array1<f64> {
        let n_jobs = self.num_threads;
        let mut result = Array1::from(vec![0., 0.]); 
        rayon::scope(|s| {
            for i in 0..n_jobs {
                // println!("starting cluster thread {i}");
                // let shared_touched_vec = self.touched_index_vec.clone();
                let shared_data = self.internal_lattice.internal_vector.clone();
                let cluster_queue_signaler = self.cluster_queue_signaler.clone();
                let cluster_done_signaler = self.cluster_done_signaler.clone();
                let ext_mag_field = self.parameters.get_b_field();
                let spin_unit = self.parameters.get_spin_unit();
                s.spawn(move |_| {
                    cluster_done_signaler
                        .send(thread_workers::thread_cluster_worker(
                                cluster_queue_signaler,
                                shared_data,
                                // shared_touched_vec,
                                target_spin,
                                balance_condition,
                                ext_mag_field,
                                spin_unit))
                        .unwrap();
                });
            } // for i in 0..n_jobs
        }); // rayon::scope(|s|

        for _ in 0..n_jobs {
            if let Ok(thread_result) = self.cluster_done_signaler.recv() {
                result = result + thread_result;
            } else {
                panic!("could not receive value from cluster_done_signaler!");
            }
        }
        // println!("{:?}", &result);
        result
    }

    pub fn get_energy(&mut self) -> f64 {
        let n_jobs = self.num_threads;
        rayon::scope(|s| {
            for i in 0..n_jobs {
                let this_thread_num = i;
                let shared_internal_vector = self.internal_lattice.internal_vector.clone();
                let energy_psum_signaler = self.energy_psum_signaler.clone();

                // start threads here
                let div = self.div;
                let rem = self.rem;
                let ext_mag_field = self.parameters.get_b_field();

                s.spawn(move |_| 
                {
                    thread_workers::thread_energy_worker(
                        n_jobs,
                        this_thread_num,
                        div,
                        rem,
                        shared_internal_vector,
                        ext_mag_field,
                        energy_psum_signaler);
                });
            } // for i in 0..n_jobs

        });
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

    pub fn get_magnetization(&self) -> f64 {
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();
        read_lock_internal_vector
            .par_iter()
            .map(|i| i.get_spin())
            .sum()
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
        
        let n = self.parameters.num_nodes() as f64;
        let spin_unit = self.parameters.get_spin_unit();
        let (temp_magnetization, initial_energy) = (self.get_magnetization()  / (n * spin_unit), self.get_energy() / (n * spin_unit));
        let initial_magnetization = temp_magnetization.abs();
        println!(
            "\nThe initial energy is {}, and the initial magnetization is {}.\n",
            initial_energy, initial_magnetization
        );

        let mut m_vec: Vec<f64> = vec![];
        let mut e_vec: Vec<f64> = vec![];
        let mut c_vec: Vec<f64> = vec![];
        let mut x_vec: Vec<f64> = vec![];

        let len_beta = beta_vec.len();
        let tot_time: u64 = times as u64 * len_beta as u64;

        let bar1 = ProgressBar::new(tot_time);
        // load the initial state at tbe begining of each new beta
        for beta_val in beta_vec.iter() {
            let balance_condition = if iteration_scheme == 1 {
                1. - consts::E.powf(-2. * beta_val * self.parameters.get_j() * self.parameters.get_spin_unit().powf(2.))
            } else {
                0.
            };
            self.internal_lattice
                // .load_state_from_file("minimum_energy_state_anneal.dat".to_string());
                .load_state_from_file(self.parameters.get_fname().to_string());
            
            let mut magnetization: f64 = initial_magnetization;
            let mut energy: f64 = initial_energy;
            let mut deltas: ArrayBase<ndarray::OwnedRepr<f64>, Dim<[usize; 1]>> = array![0., 0.];
            let mut mag_vec = vec![];
            let mut energy_vec = vec![];

            for cur_t in 0..times {
                // preform an iteration scheme of metropolis or wolff
                if iteration_scheme == 0 {
                    deltas = self.metropolis_iter(beta_val);
                } else if iteration_scheme == 1 {
                    deltas = self.wolff_iter_multithreaded(balance_condition);
                } else {
                    panic!("invalid option! expected 0 or 1, got {}", iteration_scheme);
                }
                // Beginning the data collection after a few iterations gives better
                // overall graphs becuase the system has had some time to relax and the
                // resulting data doesnt wildly fluctuate at the start of the resulting
                // data plot.
                if cur_t + 1 > ignore_n_runs && !anneal {
                    // // preform a sum and sum of squares for statistics later
                    energy += deltas[0]  / (n * spin_unit);
                    energy_vec.push(energy);
                    magnetization += deltas[1]  / (n * spin_unit);
                    mag_vec.push(magnetization);
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
        if anneal == false {
            println!("spin_energy finished!");
            self.save_files(m_vec, e_vec, c_vec, x_vec, beta_vec);
        }
        else { println!("Annealing finished!"); }
    }

    /// Purpose
    /// -------
    /// Mutates the `self.LinkedLat` lattice of spins by one Iteration of the Wolff Algorithm.
    fn wolff_iter_multithreaded(&mut self, balance_condition: f64) -> Array1<f64> {
        // self.unmark_all_nodes();
        let mut rng_spin_select = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;
        let mut deltas = array![0., 0.];

        // select a random node
        'node_selection: loop {
            target_index = rng_spin_select
                .gen_range(0..(self.parameters.get_xsize() * self.parameters.get_ysize() - 1));
            target_spin = {
                let ref this = self;
                if let Ok(value) = this.internal_lattice.internal_vector.read() {
                    match value.get(target_index) {
                        Some(node) => node.get_spin(),
                        None => 0.,
                    }
                } else {
                    panic!("Couldnt get spin")
                }
            };
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }

        if rng_spin_select.gen_range(0_f64..1_f64) < balance_condition {
            deltas = deltas + array![self.calculate_energy_change_of_suggested_flip(target_index), -target_spin];
            self.flip_node_at(target_index);

            // push the random node neighbors indicies to the work queue
            if let Ok(read_lock_internal_vector) = self.internal_lattice.internal_vector.read() {
                if let Some(target_node) = read_lock_internal_vector.get(target_index) {
                    let read_lock_nbrs = target_node.neighbors.read().unwrap();
                    for i in 0..read_lock_nbrs.len() {
                        // SAFETY: bounds checked in Driver::new, garaunteed valid return.
                        unsafe {
                            let nbr_index = *read_lock_nbrs.get_unchecked(i);
                            let nbr_spin = {
                                let ref this = self;
                                if let Ok(value) = this.internal_lattice.internal_vector.read() {
                                    match value.get(nbr_index) {
                                        Some(node) => node.get_spin(),
                                        None => 0.,
                                    }
                                } else {
                                    panic!("Couldnt get spin")
                                }
                            };
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
            deltas = deltas + self.get_cluster(target_spin, balance_condition);

        }
        deltas
    }

    /// Evolves the lattice by one iteration using the metropolis-hastings scheme.
    fn metropolis_iter(&mut self, beta: &f64) -> Array1<f64> {
        let mut rngspin = thread_rng();
        let mut rng_flip = thread_rng();
        let mut target_index: usize;
        let mut target_spin: f64;
        let mut deltas = array![0., 0.];

        // select a random node
        'node_selection: loop {
            target_index = rngspin
                .gen_range(0..(self.parameters.get_xsize() * self.parameters.get_ysize() - 1));
            target_spin = {
                let ref this = self;
                if let Ok(value) = this.internal_lattice.internal_vector.read() {
                    match value.get(target_index) {
                        Some(node) => node.get_spin(),
                        None => 0.,
                    }
                } else {
                    panic!("Couldnt get spin")
                }
            };
            if target_spin != 0. {
                break 'node_selection;
            } else {
                continue 'node_selection;
            }
        }
        // get the change in energy of the suggested spin flip
        let mut delta_energy: f64 = self.calculate_energy_change_of_suggested_flip(target_index);
        let balance_condition = consts::E.powf(-beta * self.parameters.get_j() * delta_energy);

        // flip node if accepted
        match OrderedFloat(delta_energy).cmp(&OrderedFloat(0_f64)) {
            std::cmp::Ordering::Less | std::cmp::Ordering::Equal => {
                self.flip_node_at(target_index);
                deltas = array![delta_energy, -target_spin];
            }
            std::cmp::Ordering::Greater => {
                if rng_flip.gen_range(0_f64..1_f64) < balance_condition {
                    self.flip_node_at(target_index);
                    deltas = array![delta_energy, -target_spin];
                }
            }
        }

        deltas
    }

    // fn unmark_all_nodes(&mut self) {
    //     let mut write_lock = self.internal_lattice.internal_vector.write().unwrap();
    //     for node in write_lock.iter_mut() {
    //         node.marked.write().unwrap().set_unmark();
    //     }
    // }

    pub fn calculate_energy_change_of_suggested_flip(&self, target_index: usize) -> f64 {
        let mut energy_i = 0.;
        let mut energy_f = 0.;
        let mut mag_energy_i = 0.;
        let mut mag_energy_f = 0.;
        let read_lock_internal_vector = self.internal_lattice.internal_vector.read().unwrap();

        // SAFETY: x_y was checked already to be a valid node
        if let Some(cur_target) = read_lock_internal_vector.get(target_index) {
            let cur_target_nbrs_lock = cur_target.neighbors.read().unwrap();
            // iterate through the neighbors of the suggested node to flip
            for nbrs_index_of_flip_node in cur_target_nbrs_lock.iter() {
                let target_of_cur_target =
                    read_lock_internal_vector
                        .get(*nbrs_index_of_flip_node)
                        .unwrap();
                let target_of_cur_target_nbrs_lock = target_of_cur_target.neighbors.read().unwrap();
                // cycle through the neighbors of the neighbors of the suggested node to flip (wordy, yeah)
                for nbrs_of_nbrs_index_of_flip_node in target_of_cur_target_nbrs_lock.iter() {
                    let nbr_of_nbrs_of_flip_node =
                        read_lock_internal_vector
                            .get(*nbrs_of_nbrs_index_of_flip_node)
                            .unwrap();
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
