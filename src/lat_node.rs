use ndarray::prelude::*;

pub struct SpinNode {
    spin: QuantumSpin,
    pos: Array1<f64>,
    coords: Array1<f64>,
}

impl SpinNode {

    pub fn cons_node(init_spin: f64, pos: Array1<f64>, coords: Array1<f64>) -> SpinNode {
        SpinNode {
            spin: QuantumSpin::new(Some(init_spin)),
            pos,
            coords,
        }
    }

    pub fn new(pos: Vec<i64>) -> Self {
        Self { 
            spin: QuantumSpin::new(None),
            pos: Array1::<f64>::zeros(2),
            coords: Array1::<f64>::zeros(2),
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