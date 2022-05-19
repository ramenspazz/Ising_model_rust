extern crate rand;
use ising_model::atomistic_simulation::SymmetryType;
use ising_model::*;
use ndarray::prelude::*;
use rand::{thread_rng, Rng};


fn main() {
    println!("Enter the number of iterations to evolve the system for: ");

    let times = 20;
    let start = 1.;
    let end = 0.2;
    let num_points = 10;
    let step = (end - start) / num_points as f64;
    let mut beta_list = vec![];

    for i in 1..(num_points + 1) {
        beta_list.push((i as f64) * step + start);
    }
    
    let mut c4v_driver = ising_model::atomistic_simulation::Driver::new(
        
    );

    c4v_driver.spin_energy(beta_list, times, usrin);
    c4v_driver.stop_threads();
}
