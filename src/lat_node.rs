use ndarray::prelude::*;
use std::sync::RwLock;

#[derive(Clone, PartialEq)]
pub enum StateValue {
    Unmarked,
    Processing,
    Marked,
    Pushed,
}
pub struct TrueFalse {
    marked: StateValue,
}

impl Default for TrueFalse {
    fn default() -> Self {
        Self::new()
    }
}

impl TrueFalse {
    pub fn new() -> Self {
        Self {
            marked: StateValue::Unmarked,
        }
    }

    pub fn set_marked(&mut self) {
        self.marked = StateValue::Marked;
    }

    pub fn set_processing(&mut self) {
        self.marked = StateValue::Processing;
    }

    pub fn set_unmark(&mut self) {
        self.marked = StateValue::Unmarked;
    }

    pub fn set_pushed(&mut self) {
        self.marked = StateValue::Pushed;
    }

    pub fn get_status(&self) -> StateValue {
        self.marked.clone()
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
    pub fn cons_node(
        init_spin: f64,
        pos: Array1<f64>,
        coords: Array1<f64>,
        neighbors: RwLock<Vec<usize>>,
    ) -> SpinNode {
        SpinNode {
            spin: QuantumSpin::new(Some(init_spin)),
            pos,
            coords,
            neighbors,
            marked: RwLock::new(TrueFalse::new()),
        }
    }

    pub fn new(_pos: Vec<i64>) -> Self {
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
        self.spin.get_spin().unwrap_or(0.)
    }

    pub fn set_spin(&mut self, spin: f64) {
        self.spin.set_spin(spin);
    }

    /// Flip the spin node's spin.
    pub fn flip_spin(&mut self) {
        self.spin.set_spin(self.get_spin() * -1.)
    }

    pub fn get_status(&self) -> StateValue {
        if let Ok(node_read_lock) = self.marked.read() {
            node_read_lock.get_status()
        } else {
            panic!("Unable to read node!")
        }
    }
}

pub struct QuantumSpin {
    spin: Option<f64>,
}

impl QuantumSpin {
    // construct a new spin
    pub fn new(spin: Option<f64>) -> Self {
        Self { spin }
    }

    /// Get a reference to the spin node's spin.
    #[must_use]
    pub fn get_spin(&self) -> Option<f64> {
        self.spin
    }

    /// Set the spin node's spin.
    pub fn set_spin(&mut self, spin: f64) {
        self.spin = Some(spin);
    }
}
