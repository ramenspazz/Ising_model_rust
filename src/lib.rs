/// Purpose
/// -------
/// Returns a tuple with the quotient and remainder of `dividend` / `divisor`.
/// I designed this to run in exponential jumps of the power of 2 using the
/// left bitshift operator, so it functions faster than the standard
/// implimentation of the remainder algorithm that I have seen.
/// Runs in Omega(log(n)), Theta(n), O(nlog(n))
///
/// Returns
/// -------
/// (Quotient, Remainder) : `tuple`
/// - `tuple` containing the integer quotient and remainder of
/// division.
pub fn dividend_remainder(dividend: usize, divisor: usize) -> (usize, usize) {
    // Might be necessary but for now doesnt appear to be relavant for my use
    // case. Included just incase, just uncomment and define MAX_INT and
    // MIN_INT.
    // if !(MAX_INT > area > MIN_INT) {
    //     panic('Input n and m are too large!')
    // }
    if divisor > dividend {
        return (0, dividend);
    } else if divisor == dividend {
        return (1, 0);
    } else if dividend % divisor == 0 {
        return (dividend / divisor, 0);
    }

    let mut div_power: usize = 0;
    let mut current_quotient: usize = 0;
    let mut prev: usize = 0;

    loop {
        let test_quotient = (divisor << div_power) + (divisor * current_quotient);

        if test_quotient > dividend && prev < dividend {
            if prev + divisor > dividend {
                return (current_quotient + (2 << (div_power - 2)), dividend - prev);
            }
            current_quotient += 2 << (div_power - 2);
            div_power = 0;
            prev = current_quotient * divisor;
            continue;
        } else if test_quotient < dividend {
            prev = test_quotient;
            div_power += 1;
            continue;
        }
    }
}

pub fn is_only_numbers(input: &str) -> bool {
    for value in input.as_bytes().iter() {
        if (48 <= *value && *value <= 57) || *value == 46 || *value == 45 || *value == 10 {
            continue;
        } else {
            println!("failing character is {}", *value);
            return false
        }
    }
    true
}

pub fn get_input_as_i64(msg: Option<&str>) -> i64 {
    loop {
        if msg != None {
            println!("{}", msg.unwrap());
        }
        let mut usrin = String::from("");

        std::io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");

        if is_only_numbers(&usrin) != true {
            println!("Invalid input!");
            continue;
        }

        match usrin.trim().parse::<i64>() {
            Ok(num) => return num,
            Err(_) => {
                println!("Invalid input!");
                continue;
            }
        };
    }
}

pub fn get_input_as_usize(msg: Option<&str>) -> usize {
    loop {
        if msg != None {
            println!("{}", msg.unwrap());
        }
        let mut usrin = String::from("");

        std::io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");

        if is_only_numbers(&usrin) != true {
            println!("Invalid input!");
            continue;
        }

        match usrin.trim().parse::<usize>() {
            Ok(num) => return num,
            Err(_) => {
                println!("Invalid input!");
                continue;
            }
        };
    }
}

pub fn get_input_as_f64(msg: Option<&str>) -> f64 {
    loop {
        if msg != None {
            println!("{}", msg.unwrap());
        }
        let mut usrin = String::from("");

        std::io::stdin()
            .read_line(&mut usrin)
            .expect("Failed to read line");

        if is_only_numbers(&usrin) != true {
            println!("Invalid input!");
            continue;
        }

        match usrin.trim().parse::<f64>() {
            Ok(num) => return num,
            Err(_) => {
                println!("Invalid input!");
                continue;
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use crate::atomistic_simulation::Driver;
    use crate::atomistic_simulation::SymmetryType;
    use crate::sim_params;
    use ndarray::prelude::*;
    #[test]
    fn test() {
        // let times = 20;
        let start = 1.;
        let end = 0.2;
        let num_points = 10;
        let step = (end - start) / num_points as f64;
        let mut beta_list = vec![];

        for i in 1..(num_points + 1) {
            beta_list.push((i as f64) * step + start);
        }

        let _fname = "c4v_test.dat";

        let parameters = sim_params::SimulationParameters::new(
            1.,
            4,
            4,
            SymmetryType::C4V,
            0.,
            array![[1., 0.], [0., 1.]],
            0.,
            0.5,
            "c4v_test.dat".to_string(),
        );

        let mut c4v_driver = Driver::new(parameters.clone());

        c4v_driver.load_state(parameters.get_fname());
        let mut initial_mag: f64;
        for _ in 0..6 {
            initial_mag = c4v_driver.get_magnetization();
            assert_eq!(initial_mag, 2.);
        }

        for i in 0..16 {
            let delta_energy_calc = c4v_driver.calculate_energy_change_of_suggested_flip(i);
            let energy_i = c4v_driver.get_energy();
            c4v_driver.flip_node_at(i);
            let energy_f = c4v_driver.get_energy();
            c4v_driver.flip_node_at(i);
            let delta_energy_get_energy = energy_f - energy_i;

            assert_eq!(delta_energy_calc, delta_energy_get_energy);
        }
    }
}

pub mod thread_workers;
pub mod atomistic_simulation;
pub mod lat_node;
pub mod lattice_structure;
pub mod signal_container;
pub mod sim_params;
