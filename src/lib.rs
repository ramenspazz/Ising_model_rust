pub fn DividendRemainder(dividend: usize,
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
    
        let mut test: usize = 0;
        let mut div_power: usize = 0;
        let mut quotient: usize = 0;
        let mut prev: usize = 0;

        loop {
            test = (divisor << div_power) + (divisor * quotient);
        
            if test > dividend && prev < dividend {
                if prev + divisor > dividend {
                    return (quotient + (2 << div_power - 2), dividend - prev);
                }
                quotient += 2 << (div_power - 2);
                div_power = 0;
                prev = quotient * divisor;
                continue;
            }
            else if test < dividend {
                prev = test;
                div_power += 1;
                continue;
            }
        }
    }

pub mod atomistic_simulation;
pub mod lattice_structure;
pub mod lat_node;