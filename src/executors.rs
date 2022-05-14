

extern crate rand;

use std::io;
use std::sync::mpsc::channel;
use ising_model;

// import commonly used items from the prelude:
use executors::*;
use ising_model::atomistic_simulation::SymmetryType;
use ndarray::prelude::*;


fn get_input_as_usize() -> usize {
    let mut usrin = String::from("");
    loop {
        io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");
    
        let usrin: usize = match usrin.trim().parse() {
            Ok(num) => return num,
            Err(_) => continue,
        };
    }
}

fn main() {
    println!("enter the x, then the y size of the system in terms of two basis vectors: ");
    let mut test_driver = ising_model::atomistic_simulation::Driver::new(
        get_input_as_usize(),
        get_input_as_usize(),
        SymmetryType::C4V,
        array![[1., 0.], [0., 1.]],
    );

    let times = 100000;
    let start = 0.01;
    let end = 10.;
    println!("Enter the number of beta values to calculate in the range {} to {}: ", start, end);
    let num_points = get_input_as_usize();
    let step = (end - start) / num_points as f64;
    let mut beta_list = vec![];
    for i in 1..(num_points + 1) {
        beta_list.push((i as f64) * step + start);
    }
    test_driver.spin_energy(beta_list, times);
}
