extern crate rand;
use ising_model::atomistic_simulation::SymmetryType;
use ising_model::*;
use ndarray::prelude::*;
use rand::{thread_rng, Rng};

fn main() {
    let times = get_input_as_usize(Some("Enter the number of iterations to evolve the system for: "));
    let ignore_n_runs = get_input_as_usize(Some("Enter the number of iterations to ignore so the system can relax before collecting data: "));
    let start = get_input_as_f64(Some("Enter the starting value for beta: "));
    let end = get_input_as_f64(Some("Enter the ending value for beta: "));
    let num_points = get_input_as_usize(Some(
        format!(
            "Enter the number of beta values to calculate in the range {} to {}: ",
            start, end
        )
        .as_str(),
    ));
    let step = (end - start) / num_points as f64;
    let mut beta_list = vec![];

    for i in 1..(num_points + 1) {
        beta_list.push((i as f64) * step + start);
    }

    let mut usrin: usize;
    loop {
        usrin = get_input_as_usize(Some(
            "Enter 0 to run a C3V lattice, or 1 to run C4V symmety lattice: ",
        ));
        println!("{usrin}");
        if usrin == 0 || usrin == 1 {
            break;
        } else {
            println!("Invalid input!");
            continue;
        }
    }
    let mut rand_gen = thread_rng();
    println!("Enter the J spin coupling value, x, then the y size of the system in terms of two basis vectors: ");
    if usrin == 0 {
        let fname = "c3v.dat";
        let mut c3v_driver = ising_model::atomistic_simulation::Driver::new(
            get_input_as_f64(Some("Enter J coupling: ")),
            get_input_as_usize(Some("Enter x size: ")),
            get_input_as_usize(Some("Enter y size: ")),
            SymmetryType::C3V,
            array![[1., 0.], [0.5, 3_f64.sqrt() / 2.]],
            rand_gen.gen_range(0_f64..1_f64),
            0.5,
            fname.to_string(),
        );
        c3v_driver.save_state(fname);
        loop {
            usrin = get_input_as_usize(Some("Enter 0 to run the Metropolis-Hastings algorithm, or 1 to run the Wolff algorithm: "));
            if usrin == 0 || usrin == 1 {
                break;
            } else {
                println!("Invalid input!");
                continue;
            }
        }
        c3v_driver.spin_energy(beta_list, times, usrin, ignore_n_runs);
        c3v_driver.stop_threads();
    } else if usrin == 1 {
        let fname = "c4v.dat";
        let mut c4v_driver = ising_model::atomistic_simulation::Driver::new(
            get_input_as_f64(Some("Enter J coupling: ")),
            get_input_as_usize(Some("Enter x size: ")),
            get_input_as_usize(Some("Enter y size: ")),
            SymmetryType::C4V,
            array![[1., 0.], [0., 1.]],
            rand_gen.gen_range(0_f64..1_f64),
            0.5,
            fname.to_string(),
        );
        c4v_driver.save_state(fname);
        loop {
            usrin = get_input_as_usize(Some("Enter 0 to run the Metropolis-Hastings algorithm, or 1 to run the Wolff algorithm: "));
            if usrin == 0 || usrin == 1 {
                break;
            } else {
                println!("Invalid input!");
                continue;
            }
        }
        c4v_driver.spin_energy(beta_list, times, usrin, ignore_n_runs);
        c4v_driver.stop_threads();
    } else {
        panic!("idk what happened but it shouldnt have...");
    }
}