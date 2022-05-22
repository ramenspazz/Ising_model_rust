

use crate::atomistic_simulation::{self, SymmetryType};


#[derive(Clone)]
pub struct SimulationParameters {
    J: f64,
    XSize: usize,
    YSize: usize,
    SymType: atomistic_simulation::SymmetryType,
    BField: f64,
    Basis: ndarray::Array2<f64>,
    SpinUpPercent: f64,
    SpinUnit: f64,
    Fname: String,
}

impl SimulationParameters {
    pub fn new(
        J: f64,
        XSize: usize,
        YSize: usize,
        SymType: atomistic_simulation::SymmetryType,
        BField: f64,
        Basis: ndarray::Array2<f64>,
        SpinUpPercent: f64,
        SpinUnit: f64,
        Fname: String) -> Self {
        
        assert!(XSize > 0 && YSize > 0, "X and Y size must be greater than 0!");
        Self {
            J,
            XSize,
            YSize,
            SymType,
            BField,
            Basis,
            SpinUpPercent,
            SpinUnit,
            Fname,
        }
    }

    pub fn num_nodes(&self) -> usize { self.XSize*self.YSize }
    pub fn get_xsize(&self) -> usize { self.XSize }
    pub fn get_ysize(&self) -> usize { self.YSize }
    pub fn get_spin_up_chance(&self) -> f64 { self.SpinUpPercent }
    pub fn get_spin_unit(&self) -> f64 { self.SpinUnit }
    pub fn get_basis(&self) -> ndarray::Array2<f64> { self.Basis.clone() }
    pub fn get_symtype(&self) -> SymmetryType { self.SymType }
    pub fn get_b_field(&self) -> f64 { self.BField }
    pub fn get_J(&self) -> f64 { self.J }
    pub fn get_fname(&self) -> &str { &self.Fname }
    pub fn set_fname(&mut self, fname: &str) { self.Fname = fname.to_string(); }

}