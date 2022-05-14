use std::sync::{Arc, RwLock};
use crate::lat_node::SpinNode;

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
        match self.internal_vector.write().as_mut() {
            Ok(value) => value.push(item),
            Err(_) => return,
        }
    }
}