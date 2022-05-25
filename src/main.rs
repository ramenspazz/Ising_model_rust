extern crate rand;
use ising_model::atomistic_simulation::SymmetryType;
use ising_model::*;
use ndarray::prelude::*;
use rand::{thread_rng, Rng};

fn main() {
    println!("Enter a name to save lattice state files: ");
    let mut fname = String::from("");
    std::io::stdin()
        .read_line(&mut fname)
        .expect("Failed to read line");
    
    let mut rand_gen = thread_rng();

    let mut usrin: usize;
    loop {
        usrin = get_input_as_usize(Some(
            "Enter 0 to run a C3V lattice, or 1 to run C4V symmety lattice: ",
        ));
        if usrin == 0 || usrin == 1 {
            break;
        } else {
            println!("Invalid input!");
            continue;
        }
    }

    println!("Lattice parameters and external magnetic field value:\n");
    let mut parameters = sim_params::SimulationParameters::new(
        get_input_as_f64(Some("Enter J coupling: ")),
        get_input_as_usize(Some("Enter x size: ")),
        get_input_as_usize(Some("Enter y size: ")),
        if usrin == 0 {
            SymmetryType::C3V
        } else {
            SymmetryType::C4V
        },
        get_input_as_f64(Some("Enter magnetic field value: ")),
        array![[1., 0.], [0.5, 3_f64.sqrt() / 2.]],
        rand_gen.gen_range(0_f64..1_f64),
        0.5,
        fname.to_string(),
    );

    parameters.set_fname(&fname);
    let mut driver_obj = ising_model::atomistic_simulation::Driver::new(parameters.clone());
    driver_obj.save_state(&fname);
    println!("Note, you can plot at any time with python -O plot.py");
    loop {
        let user_run_stop = get_input_as_i64(Some("Enter 0 to run a simulation, or -1 to quit: "));
        match user_run_stop.cmp(&0) {
            std::cmp::Ordering::Less => break,
            std::cmp::Ordering::Equal => {
                let times = get_input_as_usize(Some(
                    "Enter the number of iterations to evolve the system for: ",
                ));
                let ignore_n_runs = get_input_as_usize(Some("Enter the number of iterations to ignore so the system can relax before collecting data: "));
                let start = get_input_as_f64(Some("Enter the starting value for T = 1/beta: "));
                let end = get_input_as_f64(Some("Enter the ending value for T = 1/beta: "));
                let num_points = get_input_as_usize(Some(
                    format!(
                        "Enter the number of beta values to calculate in the range {} to {}: ",
                        start, end
                    )
                    .as_str(),
                ));
                let step = (1. / end - 1. / start) / num_points as f64;
                let mut beta_list = vec![];

                for i in 1..(num_points + 1) {
                    beta_list.push((i as f64) * step + 1. / start);
                }
                loop {
                    usrin = get_input_as_usize(Some("Enter 0 to run the Metropolis-Hastings algorithm, or 1 to run the Wolff algorithm: "));
                    if usrin == 0 || usrin == 1 {
                        break;
                    } else {
                        println!("Invalid input!");
                        continue;
                    }
                }
                let iteration_scheme = usrin;
                loop {
                    usrin = get_input_as_usize(Some("Enter 0 to anneal, 1 to run without anneal: "));
                    if usrin == 0 || usrin == 1 {
                        break;
                    } else {
                        println!("Invalid input!");
                        continue;
                    }
                }
                let anneal = if usrin == 0 { true } else { false };
                if anneal == true {
                    let anneal_betas = vec![100., 75., 50., 25., 10., 5., 1., 0.01];
                    println!("\nAnnealing system into a minimum energy state using beta values {:?} for 1,000,000 iterations for each beta value in the Metropolis-Hastings scheme.\n", &anneal_betas);
                    driver_obj.spin_energy(
                        anneal_betas,
                        1_000_000,
                        0,
                        1_000_000,
                        true,
                    );
                    driver_obj.save_state(&fname);
                }
                else {
                    loop {
                        println!("Enter 0 to load state from user specified filename or 1 to load from {}: ", &parameters.get_fname());
                        usrin = get_input_as_usize(None);
                        if usrin == 0 || usrin == 1 {
                            break;
                        } else {
                            println!("Invalid input!");
                            continue;
                        }
                    }
                    let load_state = if usrin == 0 { true } else { false };
                    if load_state == true {
                        let mut load_fname = String::from("");
                        std::io::stdin()
                            .read_line(&mut load_fname)
                            .expect("Failed to read line");
                        driver_obj.load_state(&load_fname);
                    } else {
                        driver_obj.load_state(&fname);
                    }
                }
                driver_obj.spin_energy(beta_list, times, iteration_scheme, ignore_n_runs, false);
            }
            std::cmp::Ordering::Greater => {
                println!("invalid input!");
                continue;
            }
        }
    }
    // end
}
