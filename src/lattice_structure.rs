use crate::is_only_numbers;
use crate::lat_node::SpinNode;
use std::io::BufRead;
use std::{
    io::Write,
    sync::{Arc, RwLock},
};

pub struct Lattice {
    pub internal_vector: Arc<RwLock<Vec<SpinNode>>>,
}

impl Lattice {
    pub fn new(internal_vec: Vec<SpinNode>) -> Self {
        Self {
            internal_vector: Arc::new(RwLock::new(internal_vec)),
        }
    }

    pub fn push(&mut self, item: SpinNode) {
        if let Ok(value) = self.internal_vector.write().as_mut() {
            value.push(item)
        }
    }

    fn read_line(&mut self, file_name: &str) -> Result<(), std::io::Error> {
        // open target file
        let file = std::fs::File::open(&file_name)?;

        // uses a reader buffer
        let mut reader = std::io::BufReader::new(file);
        let mut line = String::new();
        let mut cursor: usize = 0;
        let mut write_lock_internal_vector = self.internal_vector.write().unwrap();

        loop {
            match reader.read_line(&mut line) {
                Ok(bytes_read) => {
                    // EOF: save last file address to restart from this address for next run
                    if bytes_read == 0 {
                        break;
                    }

                    if is_only_numbers(&line) {
                        println!("Invalid input!");
                        panic!("Invalid entry in file! exiting.")
                    }

                    match line.trim().parse::<f64>() {
                        Ok(num) => {
                            if let Some(cur_node) = write_lock_internal_vector.get_mut(cursor) {
                                cur_node.set_spin(num);
                                // println!("node {} initilized to spin {}", &cursor, &num);
                            }
                        }
                        Err(_) => {
                            println!("Invalid input!");
                            continue;
                        }
                    };

                    // do not accumulate data
                    line.clear();
                    cursor += 1;
                }
                Err(err) => {
                    return Err(err);
                }
            };
        }
        Ok(())
    }

    pub fn load_state_from_file(&mut self, fname: String) {
        self.read_line(fname.as_str()).unwrap();
    }

    pub fn export_state_to_file(&self, fname: &str) {
        // setup saving the output file
        let path = std::path::Path::new(fname);
        let display = path.display();
        // Write output data to a file
        let mut file = match std::fs::File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => {
                file // forward file to outer scope
            }
        };
        // start writing data
        let read_lock_internal_vector = self.internal_vector.read().unwrap();
        for index in 0..read_lock_internal_vector.len() {
            let write_str = format!("{}\n", unsafe {
                read_lock_internal_vector.get_unchecked(index).get_spin()
            });
            match file.write_all(write_str.as_bytes()) {
                Err(why) => panic!("couldn't write to {}: {}", display, why),
                Ok(_) => {
                    continue;
                }
            }
        }
    }
}
