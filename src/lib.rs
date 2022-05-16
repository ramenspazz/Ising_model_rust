pub fn dividend_remainder(dividend: usize,
    divisor: usize) -> (usize, usize) {
        // Purpose
        // -------
        // Returns a list with the quotient and remainder of `dividend` / `divisor`.
        // I designed this to run in exponential jumps of the power of 2 using the
        // left bitshift operator, so it functions faster than the standard
        // implimentation of the remainder algorithm that I have seen.
        // Runs in Omega(log(n)), Theta(n), O(nlog(n))
        // 
        // Returns
        // -------
        // [Quotient, Remainder] : `list`
        // - ``list` containing the integer quotient and remainder of
        // division.
    
        // Might be necessary but for now doesnt appear to be relavant for my use
        // case. Included just incase, just uncomment and define MAX_INT and
        // MIN_INT.
        // if !(MAX_INT > area > MIN_INT) {
        //     panic('Input n and m are too large!')
        // }
        if divisor > dividend {
            return (0, dividend);
        }
        else if divisor == dividend {
            return (1, 0);
        }
        else if dividend % divisor == 0 {
            return (dividend / divisor, 0);
        }
    
        let mut div_power: usize = 0;
        let mut current_quotient: usize = 0;
        let mut prev: usize = 0;

        loop {
            let test_quotient = (divisor << div_power) + (divisor * current_quotient);
        
            if test_quotient > dividend && prev < dividend {
                if prev + divisor > dividend {
                    return (current_quotient + (2 << div_power - 2), dividend - prev);
                }
                current_quotient += 2 << (div_power - 2);
                div_power = 0;
                prev = current_quotient * divisor;
                continue;
            }
            else if test_quotient < dividend {
                prev = test_quotient;
                div_power += 1;
                continue;
            }
        }
    }

// https://play.rust-lang.org/?gist=9ca489c1de25feb2f09a1770e40b62b1&version=stable
pub fn is_only_numbers(input: &str) -> bool {
    let chars_are_numeric: Vec<bool> = input.chars().map(|c|c.is_numeric()).collect();
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

pub mod atomistic_simulation;
pub mod lattice_structure;
pub mod lat_node;
