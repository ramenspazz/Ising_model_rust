use crate::atomistic_simulation::{self, SymmetryType};

#[derive(Clone)]
pub struct SimulationParameters {
    j: f64,
    xsize: usize,
    ysize: usize,
    symtype: atomistic_simulation::SymmetryType,
    b_field: f64,
    basis: ndarray::Array2<f64>,
    spin_up_percent: f64,
    spin_unit: f64,
    fname: String,
}

impl SimulationParameters {
    pub fn new(
        j: f64,
        xsize: usize,
        ysize: usize,
        symtype: atomistic_simulation::SymmetryType,
        b_field: f64,
        basis: ndarray::Array2<f64>,
        spin_up_percent: f64,
        spin_unit: f64,
        fname: String,
    ) -> Self {
        assert!(
            xsize > 0 && ysize > 0,
            "X and Y size must be greater than 0!"
        );
        Self {
            j,
            xsize,
            ysize,
            symtype,
            b_field,
            basis,
            spin_up_percent,
            spin_unit,
            fname,
        }
    }

    pub fn num_nodes(&self) -> usize {
        self.xsize * self.ysize
    }
    pub fn get_xsize(&self) -> usize {
        self.xsize
    }
    pub fn get_ysize(&self) -> usize {
        self.ysize
    }
    pub fn get_spin_up_chance(&self) -> f64 {
        self.spin_up_percent
    }
    pub fn get_spin_unit(&self) -> f64 {
        self.spin_unit
    }
    pub fn get_basis(&self) -> ndarray::Array2<f64> {
        self.basis.clone()
    }
    pub fn get_symtype(&self) -> SymmetryType {
        self.symtype
    }
    pub fn get_b_field(&self) -> f64 {
        self.b_field
    }
    pub fn get_j(&self) -> f64 {
        self.j
    }
    pub fn get_fname(&self) -> &str {
        &self.fname
    }
    pub fn set_fname(&mut self, fname: &str) {
        self.fname = fname.to_string();
    }
}
