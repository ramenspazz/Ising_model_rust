use ndarray::prelude::*;
use std::sync::RwLock;

pub struct TrueFalse {
    marked: bool,
}

impl TrueFalse {
    pub fn new() -> Self { Self { marked: false } }

    pub fn flip(&mut self) {
        if self.marked == true {
            self.marked = false;
        } else {
            self.marked = true;
        }
    }

    pub fn unmark(&mut self) { self.marked = false; }

    pub fn get_marked(&self) -> bool {
        self.marked.clone()
    }

    pub fn mark(&mut self) {
        self.marked = true;
    }

}

pub struct SpinNode {
    spin: QuantumSpin,
    pos: Array1<f64>,
    coords: Array1<f64>, // included for plotting purposes
    pub neighbors: RwLock<Vec<usize>>,
    pub marked: RwLock<TrueFalse>,
}

impl SpinNode {

    pub fn cons_node(init_spin: f64, pos: Array1<f64>, coords: Array1<f64>, neighbors: RwLock<Vec<usize>>) -> SpinNode {
        SpinNode {
            spin: QuantumSpin::new(Some(init_spin)),
            pos,
            coords,
            neighbors,
            marked: RwLock::new(TrueFalse::new()),
        }
    }

    pub fn new(pos: Vec<i64>) -> Self {
        Self { 
            spin: QuantumSpin::new(None),
            pos: Array1::<f64>::zeros(2),
            coords: Array1::<f64>::zeros(2),
            neighbors: RwLock::new(vec![]),
            marked: RwLock::new(TrueFalse::new()),
            } 
    }

    /// Set the spin node's pos.
    pub fn set_pos(&mut self, new_pos: Array1<f64>) {
        self.pos = new_pos;
    }
    
    /// Get a reference to the spin node's pos.
    #[must_use]
    pub fn get_pos(&self) -> Array1<f64> {
        self.pos.clone()
    }

    pub fn get_spin(&self) -> f64 {
        match self.spin.get_spin() {
            Some(spin) => spin,
            None => 0.,
        }
    }

    pub fn set_spin(&mut self, spin: f64) {
        self.spin.set_spin(spin);
    }

    /// Flip the spin node's spin.
    pub fn flip_spin(&mut self) {
        self.spin.set_spin(self.get_spin() * -1.)
    }
}

pub struct QuantumSpin {
    spin: Option<f64>,
}

impl QuantumSpin {
    // construct a new spin
    pub fn new(spin: Option<f64>) -> Self { Self { spin } }
    
    /// Get a reference to the spin node's spin.
    #[must_use]
    pub fn get_spin(&self) -> Option<f64> {
        self.spin.clone()
    }
    
    /// Set the spin node's spin.
    pub fn set_spin(&mut self, spin: f64) {
        self.spin = Some(spin);
    }
}