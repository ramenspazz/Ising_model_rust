/// Purpose
/// -------
/// Returns a list with the quotient and remainder of `dividend` / `divisor`.
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

// https://play.rust-lang.org/?gist=9ca489c1de25feb2f09a1770e40b62b1&version=stable
pub fn is_only_numbers(input: &str) -> bool {
    let chars_are_numeric: Vec<bool> = input.chars().map(|c| c.is_numeric()).collect();
    !chars_are_numeric.contains(&false)
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

        if is_only_numbers(&usrin) || usrin.contains('.') {
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

        if is_only_numbers(&usrin) {
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
    use ndarray::prelude::*;
    #[test]
    fn test2() {
        // println!("Enter the number of iterations to evolve the system for: ");
    
        let times = 20;
        let start = 1.;
        let end = 0.2;
        let num_points = 10;
        let step = (end - start) / num_points as f64;
        let mut beta_list = vec![];
    
        for i in 1..(num_points + 1) {
            beta_list.push((i as f64) * step + start);
        }
    
        let fname = "c4v_test.dat";
        let mut c4v_driver = Driver::new(
            1.,
            4,
            4,
            SymmetryType::C4V,
            array![[1., 0.], [0., 1.]],
            0.,
            0.5,
            fname.to_string(),
        );

        c4v_driver.load_state(fname);

        for _ in 0..6 {
            println!(
                "initial magnitization = {}",
                &c4v_driver.get_magnitization()
            );
            println!("initial energy = {}", &c4v_driver.get_energy());
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
        c4v_driver.stop_threads();
    }
}

pub mod signal_container;
pub mod atomistic_simulation;
pub mod lat_node;
pub mod lattice_structure;
