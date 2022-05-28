extern crate rand;
use ising_model::atomistic_simulation::SymmetryType;
use ising_model::*;
use ndarray::prelude::*;
use std::f64::consts;

fn main() {
    println!("Enter a name to save lattice state files: ");
    let mut fname_in = String::from("");
    std::io::stdin()
        .read_line(&mut fname_in)
        .expect("Failed to read line");
    let temp = fname_in.to_string();
    let fname = temp.trim();

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
    let symmetry = if usrin == 0 {
        SymmetryType::C3V
    } else {
        SymmetryType::C4V
    };
    let mut parameters = sim_params::SimulationParameters::new(
        1.,
        get_input_as_usize(Some("Enter x size: ")),
        get_input_as_usize(Some("Enter y size: ")),
        symmetry,
        0.,
        match symmetry {
            SymmetryType::C3V => array![[1., 0.], [0.5, 3_f64.sqrt() / 2.]],
            SymmetryType::C4V => array![[1., 0.], [0., 1.]],
            SymmetryType::C6V => todo!(),
        },
        get_input_as_f64(Some("Enter decimal percentage of spins to start spin up: ")),
        0.5,
        fname.to_string(),
    );
    //c4v_50000iter_0_01-20T_40betavals_wolff_90percent_up

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
                let beta_list = logspace(start, end, num_points);
                driver_obj.update_j(get_input_as_f64(Some("Enter a value for the exchange constant")));
                driver_obj.update_b(get_input_as_f64(Some("Enter magnetic field value: ")));
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
                        println!("Enter a lattice state file name to load: ");
                        let mut fname_in = String::from("");
                        std::io::stdin()
                            .read_line(&mut fname_in)
                            .expect("Failed to read line");
                        let temp = fname_in.to_string();
                        let load_fname = temp.trim();
                        driver_obj.load_state(&load_fname);
                    } else {
                        driver_obj.load_state(&fname);
                    }
                }
                let start_time = std::time::Instant::now();
                driver_obj.spin_energy(beta_list, times, iteration_scheme, ignore_n_runs, false);
                let elapsed_time = start_time.elapsed();
                println!("Spin energy finished in {}.", elapsed_time.as_secs_f32());
            }
            std::cmp::Ordering::Greater => {
                println!("invalid input!");
                continue;
            }
        }
    }
    // end
}

fn logspace(min:f64,max: f64, logbins: usize) -> Vec<f64> {
    let logmin = min.log(consts::E);
    let logmax = max.log(consts::E);
    let delta = (logmax - logmin) / logbins as f64;
    let mut accdelta: f64 = 0.;
    let mut v = vec![];
    for _ in 0..(logbins + 1)
    {
        v.push(consts::E.powf(logmin + accdelta));
        accdelta += delta;
    }
    v
}
