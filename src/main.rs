extern crate rand;

use std::io;
use std::sync::mpsc::channel;
use ising_model;

// import commonly used items from the prelude:
use executors::*;
use ising_model::atomistic_simulation::SymmetryType;
use ndarray::prelude::*;

// https://play.rust-lang.org/?gist=9ca489c1de25feb2f09a1770e40b62b1&version=stable
fn is_only_numbers(input: &str) -> bool {
    let chars_are_numeric: Vec<bool> = input.chars().map(|c|c.is_numeric()).collect();
    !chars_are_numeric.contains(&false)
}

fn get_input_as_usize(msg: Option<&str>) -> usize {
    loop {
        if msg != None {
            println!("{}", msg.unwrap());
        }
        let mut usrin = String::from("");

        io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");

        if is_only_numbers(&usrin) { 
            println!("Invalid input!");
            continue
        }
        else if usrin.contains(".") {
            println!("Invalid input!");
            continue
        }

        match usrin.trim().parse::<usize>() {
            Ok(num) => return num,
            Err(_) => {
                println!("Invalid input!");
                continue
            },
        };
    }
}

fn get_input_as_f64(msg: Option<&str>) -> f64 {
    loop {
        if msg != None {
            println!("{}", msg.unwrap());
        }
        let mut usrin = String::from("");

        io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");

        if is_only_numbers(&usrin) { 
            println!("Invalid input!");
            continue
        }

        match usrin.trim().parse::<f64>() {
            Ok(num) => return num,
            Err(_) => {
                println!("Invalid input!");
                continue
            },
        };
    }
}

fn main() {
    println!("Enter the number of iterations to evolve the system for: ");

    let times = get_input_as_usize(None);
    let start = get_input_as_f64(Some("Enter the starting value for beta: "));
    let end = get_input_as_f64(Some("Enter the ending value for beta: "));
    let num_points = get_input_as_usize(Some(format!("Enter the number of beta values to calculate in the range {} to {}: ", start, end).as_str()));
    let step = (end - start) / num_points as f64;
    let mut beta_list = vec![];

    for i in 1..(num_points + 1) {
        beta_list.push((i as f64) * step + start);
    }

    let mut usrin: usize = 0;
    loop {
        usrin = get_input_as_usize(Some("Enter 0 to run a C3V lattice, or 1 to run C4V symmety lattice: "));
        println!("{usrin}");
        if usrin == 0 || usrin == 1 {
            break
        }
        else {
            println!("Invalid input!");
            continue
        }
    }

    println!("Enter the J spin coupling value, x, then the y size of the system in terms of two basis vectors: ");
    if usrin == 0 {
        let mut c3v_driver = ising_model::atomistic_simulation::Driver::new(
            get_input_as_f64(Some("Enter J coupling: ")),
            get_input_as_usize(Some("Enter x size: ")),
            get_input_as_usize(Some("Enter y size: ")),
            SymmetryType::C3V,
            array![[1., 0.], [0.5, (3 as f64).sqrt()/2.]],
        );
        loop {
            usrin = get_input_as_usize(Some("Enter 0 to run the Metropolis-Hastings algorithm, or 1 to run the Wolff algorithm: "));
            if usrin == 0 || usrin == 1 {
                break
            }
            else {
                println!("Invalid input!");
                continue
            }
        }
        c3v_driver.spin_energy(beta_list, times, usrin);
    } else if usrin == 1 {
        let mut c4v_driver = ising_model::atomistic_simulation::Driver::new(
            get_input_as_f64(Some("Enter J coupling: ")),
            get_input_as_usize(Some("Enter x size: ")),
            get_input_as_usize(Some("Enter y size: ")),
            SymmetryType::C4V,
            array![[1., 0.], [0., 1.]],
        );
        loop {
            usrin = get_input_as_usize(Some("Enter 0 to run the Metropolis-Hastings algorithm, or 1 to run the Wolff algorithm: "));
            if usrin == 0 || usrin == 1 {
                break
            }
            else {
                println!("Invalid input!");
                continue
            }
        }
        c4v_driver.spin_energy(beta_list, times, usrin);
    }
    else {
        panic!("idk what happened but it shouldnt have...");
    }


    // test_driver.stop_threads();
}
