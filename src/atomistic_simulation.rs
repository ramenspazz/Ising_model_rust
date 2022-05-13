use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use core::ops::Range;
use indicatif::ProgressBar;
use std::f64::consts;
use crossbeam_channel::Sender;
use rand::thread_rng;
use rand::Rng;
use std::sync::Arc;
use std::sync::RwLock;
use crossbeam_channel::Receiver;
use executors::parker::SmallThreadData;
use ndarray::prelude::*;
use executors::Executor;
use crossbeam_utils::thread;
use crossbeam_channel::bounded;
use executors::crossbeam_workstealing_pool;
use executors::parker::StaticParker;
use crate::atomistic_simulation::crossbeam_workstealing_pool::ThreadPool;
use crate::DividendRemainder;
use crate::lattice_structure::lattice;
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
    internal_lattice: lattice,
    basis: RwLock<Array2<f64>>,
    x_size: usize,
    y_size: usize,
    total_nodes: usize,
    sym_type: SymmetryType,
    num_threads: usize,
    div: usize,
    rem: usize,
    threads_started: bool,
    mag_pool: ThreadPool<StaticParker<SmallThreadData>>,
    energy_pool: ThreadPool<StaticParker<SmallThreadData>>,
    tx_psum: Sender<f64>,
    rx_psum: Receiver<f64>,
    tx_go_psum: Sender<SignalType>,
    rx_go_psum: Receiver<SignalType>,
}

impl Driver {
    pub fn new(x_size: usize, y_size: usize, sym_type: SymmetryType, inner_basis: Array2<f64>) -> Self { 
        assert!(x_size > 0, "size must be greater than zero!");
        assert!(y_size > 0, "size must be greater than zero!");
        let cpu_num = num_cpus::get();
        let (div, rem) = DividendRemainder(x_size*y_size, cpu_num);
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
                        if genval >= 0.5 {
                            give_map.push(SpinNode::cons_node(0.5, array![i as f64, j as f64], *cur_coord));
                        } else {
                            give_map.push(SpinNode::cons_node(-0.5, array![i as f64, j as f64], *cur_coord));
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
                // When this is repeated, a lattice is created and we get (real end product looks better than comment):
                //
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //                                         \_/ \_/ \_/ \
                //                                         / \_/ \_/ \_/
                //
                let cur_basis = inner_basis.clone();
                let b1 = array![cur_basis[(0,0)], cur_basis[(0,1)]];
                let b2 = array![cur_basis[(1,0)], cur_basis[(1,1)]];
                let b3 = array![-b2[0], b2[1]];
                let mut cur_coord: Box<Array1<f64>> = Box::new(array![0., 0.]);
                for i in 0..x_size {
                    for j in 0..y_size {
                        // randomization conditions would go here or somewhere near here.
                        let genval: f64 = rng.gen();
                        if genval >= 0.5 {
                            give_map.push(SpinNode::cons_node(0.5, array![i as f64, j as f64], *cur_coord.clone()));
                        } else {
                            give_map.push(SpinNode::cons_node(-0.5, array![i as f64, j as f64], *cur_coord.clone()));
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
                let cur_basis = inner_basis.clone();
                let b1 = array![cur_basis[(0,0)], cur_basis[(0,1)]];
                let b2 = array![cur_basis[(1,0)], cur_basis[(1,1)]];
                let b3 = array![-b2[0], b2[1]];
                let mut cur_coord = array![0., 0.];
                for i in 0..x_size {
                    for j in 0..y_size {
                        assert_eq!(i, i);  // TODO
                        assert_eq!(j, j);
                    }
                }
            },
        }
        println!("Done!");
        let (tx_psum, rx_psum) = bounded(cpu_num);
        let (tx_go_psum, rx_go_psum) = bounded(cpu_num);
        // let (tx_energy, rx_energy) = bounded(cpu_num);
        // let (tx_go_energy, rx_go_energy) = bounded(cpu_num);
        Self {
            internal_lattice: lattice::new(give_map),
            basis: RwLock::new(inner_basis),
            x_size,
            y_size,
            total_nodes: x_size * y_size,
            sym_type,
            num_threads: cpu_num,
            div,
            rem,
            threads_started: false,
            mag_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            energy_pool: crossbeam_workstealing_pool::small_pool(cpu_num),
            tx_psum,
            tx_go_psum,
            rx_psum,
            rx_go_psum,

        }
    }

    fn start_threads(&mut self) {
        if self.threads_started != true {
            let n_jobs = self.num_threads.clone();
            // let (finished_tx, finished_rx) = channel();
            for i in 0..n_jobs {
                println!("starting thread {i}");
                let shared_data = self.internal_lattice.internal_vector.clone();
                let tx_psum = self.tx_psum.clone();
                let rx_go_psum = self.rx_go_psum.clone();
                let range = 0..2;

                self.mag_pool.execute(move || {

                    'outer : loop {
                        let mut psum = 0.;
                        match rx_go_psum.recv().unwrap() {
                            SignalType::SigStart => {
                                match shared_data.read().as_deref() {
                                    Ok(value) => {
                                        for item in value {
                                            psum += item.get_spin();
                                        }
                                    },
                                    Err(_) => panic!("couldnt read shared data in thread!"), 
                                }
                                tx_psum.send(psum).unwrap();
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
        self.threads_started = true;
    }

    fn stop_threads(&self) {
        // for _ in 0..self.num_threads {
        //     self.tx_go.send(SignalType::SigStop).unwrap();
        // }
        println!("I am not even going to bother with figuring out how to do this as the threads just wont die no matter what I do, so live with all the errors, though know that they are really just the OS complaining that the threads are still alive or were already killed. Not much I can do with my current understanding of this language.\n");
        let mut count = 0;
        loop {
            match self.mag_pool.shutdown_borrowed() {
                Ok(_) => break,
                Err(why) => {
                    count += 1;
                    println!("{}", why);
                    if count == 4 {
                        break;
                    }
                    continue;
                },
            }
        }
        println!("threads closed, or not. They may throw a panic, but this technically kills them so its all good. just ignore the incoming text wall...\n")
    }


    pub fn get_spin_at(&self, index: usize) -> f64 {
        match self.internal_lattice.internal_vector.read() {
            Ok(value) => match value.get(index) {
                Some(node) => node.get_spin(),
                None => 0.,
            },
            Err(_) => panic!("couldnt get spin"),
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

                        },
                        SymmetryType::C4V => {
                            for j in 0..4 {
                                if let Some(valid_index) = indexmod(i, j, self.x_size, self.y_size, self.sym_type) {
                                    psum += self.get_spin_at(valid_index);
                                }
                            }
                        },
                        SymmetryType::C6V => {

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
            let range = (i * self.div)..((i+1) * self.div);
            energy += self.energy_worker(range);
        }
        return energy;
    }

    pub fn spin_energy(&mut self, beta_list: Vec<f64>, times: usize) {
        self.start_threads();
        let magnitization = self.get_magnitization();
        let return_vec = Box::new(vec![1.,2.]);
        
        let mut energy = self.get_energy();
        println!("The initial energy is {}", energy);
        let mut energy_vec: Vec<f64> = vec![];
        
        // stfu rust no one cares about your idiomatic naming conventions.
        // I will continue to name these with a capital CammelCaseYouShmuckLol
        // cuz #physics.
        let mut Ms: Vec<f64> = vec![magnitization, magnitization.powf(2.)];
        let mut Es: Vec<f64> = vec![energy, energy.powf(2.)];

        let mut M: Vec<f64> = vec![];
        let mut E: Vec<f64> = vec![];
        let mut C: Vec<f64> = vec![];
        let mut X: Vec<f64> = vec![];
        
        let len_beta = beta_list.len();
        let tot_time: u64 = times as u64 * len_beta as u64;

        println!("Creating bar1 for spin-energy calculation");
        let bar1 = ProgressBar::new(tot_time);

        for beta_val in beta_list {
            // let mut run_avg: Vec<f64> = vec![];
            // let mut cur_run_sum = 0.; 
            for i in 0..times {
                energy = self.metropolis_iter(&beta_val, energy);
                let cur = self.get_magnitization();
                Ms[0] += magnitization.abs();
                Ms[1] += magnitization.powf(2.);
                Es[0] += energy;
                Es[1] += energy.powf(2.);
                // cur_run_sum += cur;
                M.push(cur / self.total_nodes as f64);
                E.push(energy / self.total_nodes as f64);  // TODO add in actual value
                C.push(Es[1] / ((self.total_nodes as f64 -1.)*(times as f64).powf(2.)) - 
                Es[0].powf(2.) / (self.total_nodes as f64 -1.)*(times as f64).powf(2.) *
                       beta_val.powf(2.));
                       X.push(Ms[1] / ((self.total_nodes as f64 -1.)*(times as f64).powf(2.)) - 
                       Ms[0].powf(2.) / (self.total_nodes as f64 -1.)*(times as f64).powf(2.) *
                       beta_val);
                       bar1.inc(1);
                    }
                energy_vec.push(energy.clone());
            }
        bar1.abandon();
        self.stop_threads();
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
        for item in M.into_iter() {
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
        for item in E.into_iter() {
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
        for item in X.into_iter() {
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
                println!("Created files sucessfully, now writing data, please wait...");
                file // forward file to outer scope
            },
        };
        for item in C.into_iter() {
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

    // fn wolff_iter(&mut self, beta: &f64, energy: f64) -> f64 {
    //     //
    //     // Purpose
    //     // -------
    //     // Mutates the self.LinkedLat lattice of spins by one Iteration of
    //     // the Wolff Algorithm.
    //     //
    //     // if the required threads are not alive already, launch them.
    //     let balcond = 1. - consts::E.powf(-2.*beta);
    //     let mut rngx = thread_rng();
    //     let mut rngy = thread_rng();
    //     let mut rng_flip = thread_rng();
    //     let mut x_y = 0;
    //     let mut target_spin = 0.;
    //     // pick random point
    //     'node_selection : loop {
    //         x_y = rngx.gen_range(0..self.x_size) + rngy.gen_range(0..self.y_size) * self.x_size;
    //         target_spin = self.get_spin_at(x_y);
    //         // println!("{}", &x_y);
    //         if target_spin != 0. { break 'node_selection; }
    //         else { continue 'node_selection; }
    //     }

    //     // push the random node to the work queue
    //     for nbr in rand_node:
    //         if nbr.get_spin() == rand_node.get_spin():
    //             self.work_queue_path.put(nbr.get_index())
    //             nbr.unmark_node()
    //     // wait for cluster to be generated
    //     self.GetCluster(balcond)
    //     try:
    //         // evaluate energy change path integral (discrete sum? lol)
    //         dE = float64(0)
    //         cur = self.LinkedLat[self.cluster_queue.get(block=False)]
    //         S_i = cur.get_spin()
    //         while True:
    //             if random.random() < balcond:
    //                 nbr_sum = float64(0)
    //                 _dE = 2*S_i*nbr_sum*self.J
    //                 for nbr in cur:
    //                     nbr_sum += nbr.get_spin()
    //                 cur.flip_spin()
    //                 dE += _dE
    //             cur = self.LinkedLat[self.cluster_queue.get(block=False)]
    //     except queue.Empty:
    //         // exit while loop when queue is empty and end the current
    //         // Iteration
    //         pass

    //     return energy + dE;
    // }

    fn metropolis_iter(&mut self, beta: &f64, energy: f64) -> f64 {
        let mut rngx = thread_rng();
        let mut rngy = thread_rng();
        let mut rng_flip = thread_rng();
        let mut x_y = 0;
        let mut target_spin = 0.;
        // Evolves the lattice by one iteration.
        'node_selection : loop {
            x_y = rngx.gen_range(0..self.x_size) + rngy.gen_range(0..self.y_size) * self.x_size;
            target_spin = self.get_spin_at(x_y);
            // println!("{}", &x_y);
            if target_spin != 0. { break 'node_selection; }
            else { continue 'node_selection; }
        }

        let mut dE: f64 = 0.;
        let mut nbr_E: f64 = 0.;
        match self.sym_type {
            SymmetryType::C4V => {
                // calculate the sum of energy of the neighbors
                for i in 0..4 {
                    nbr_E += self.get_spin_at(indexmod(x_y, i, self.x_size, self.y_size, self.sym_type).unwrap_or_default());
                }
                dE = 2.*target_spin*nbr_E;
                if dE <= 0. {
                    self.flip_node_at(x_y);
                    // dE = 2.*target_spin*nbr_E;
                }
                else if rng_flip.gen_range(0_f64..1_f64) < consts::E.powf(-beta*dE) {
                    self.flip_node_at(x_y);
                    // dE = 2.*target_spin*nbr_E;
                }
                else {
                    dE = 0.;
                }
            }, 
            SymmetryType::C3V => {

            },
            SymmetryType::C6V => {

            }
        }
        return energy + dE;
    }
}
    