use std::sync::{RwLock, Arc};
use ndarray::prelude::*;
use rand::{thread_rng, Rng};

use crate::{signal_container::SignalContainer, lat_node::{SpinNode, StateValue}};

pub fn thread_cluster_worker(
    cluster_queue_signaler: SignalContainer<usize>,
    shared_data: Arc<RwLock<Vec<SpinNode>>>,
    shared_touched_vec: Arc<RwLock<Vec<usize>>>,
    target_spin: f64,
    balance_condition: f64,
    ext_mag_field: f64,
    spin_unit: f64,
) -> Array1<f64> {
    let mut rng = thread_rng();
    let mut deltas = Array1::from(vec![0., 0.]);
    'cluster_loop : loop {
        if let Ok(node_index) = cluster_queue_signaler.try_recv() {
            let read_lock_internal_vec = shared_data.read().unwrap();
            if read_lock_internal_vec
                .get(node_index)
                .unwrap()
                .get_status()
                != StateValue::Unmarked
            {
                continue 'cluster_loop;
            } else {
                let mut write_lock_touched_index_vector =
                    shared_touched_vec.write().unwrap();
                write_lock_touched_index_vector.push(node_index);
                drop(write_lock_touched_index_vector);
                read_lock_internal_vec
                    .get(node_index)
                    .unwrap()
                    .marked
                    .write()
                    .unwrap()
                    .set_marked();
                if target_spin
                    != read_lock_internal_vec
                        .get(node_index)
                        .unwrap()
                        .get_spin()
                {
                    continue 'cluster_loop;
                }
            }
            drop(read_lock_internal_vec);
            if rng.gen_range(0_f64..1_f64) < balance_condition {
                let mut energy_i = 0.;
                let mut energy_f = 0.;
                let mut mag_energy_i = 0.;
                let mut mag_energy_f = 0.;
                // iterate through the neighbors of the suggested node to flip
                for nbrs_index_of_flip_node in shared_data
                    .read()
                    .unwrap()
                    .get(node_index)
                    .unwrap()
                    .neighbors
                    .read()
                    .unwrap()
                    .iter()
                {
                    let read_lock_internal_vec =
                        shared_data.read().unwrap();
                    let target_of_cur_target = read_lock_internal_vec
                        .get(*nbrs_index_of_flip_node)
                        .unwrap();
                    let target_of_cur_target_nbrs_lock =
                        target_of_cur_target.neighbors.read().unwrap();
                    // cycle through the neighbors of the neighbors of the suggested node to flip (wordy, yeah)
                    for nbrs_of_nbrs_index_of_flip_node in
                        target_of_cur_target_nbrs_lock.iter()
                    {
                        let read_lock_internal_vec =
                            shared_data.read().unwrap();
                        let nbr_of_nbrs_of_flip_node =
                            read_lock_internal_vec.get(
                                *nbrs_of_nbrs_index_of_flip_node,
                            ).unwrap();
                        // if nbrs_of_nbrs_index_of_flip_node is the same as the target_index, calculate the E_i and E_f
                        if *nbrs_of_nbrs_index_of_flip_node
                            == node_index
                        {
                            // get the energy using a flipped value of spin for nbr_of_nbrs_of_flip_node
                            energy_f -= -nbr_of_nbrs_of_flip_node
                                .get_spin()
                                * target_of_cur_target.get_spin();
                            mag_energy_f -= -ext_mag_field
                                * target_of_cur_target.get_spin();
                        } else {
                            energy_f -= nbr_of_nbrs_of_flip_node
                                .get_spin()
                                * target_of_cur_target.get_spin();
                            mag_energy_i -= ext_mag_field
                                * target_of_cur_target.get_spin();
                        }
                        // get the regular energy
                        energy_i -= nbr_of_nbrs_of_flip_node.get_spin()
                            * target_of_cur_target.get_spin();
                    }
                }
                deltas = deltas + array![(energy_f - energy_i) / spin_unit
                    + (mag_energy_f - mag_energy_i), 0.];
                deltas = deltas + array![0., -target_spin];
                let mut write_lock_internal_vector =
                    shared_data.write().unwrap();
                // flip the spin
                write_lock_internal_vector
                    .get_mut(node_index)
                    .unwrap()
                    .flip_spin();
                drop(write_lock_internal_vector);
                let read_lock_internal_vec =
                    shared_data.read().unwrap();
                let nbs_lock_of_target = read_lock_internal_vec
                    .get(node_index)
                    .unwrap()
                    .neighbors
                    .read()
                    .unwrap();
                for nbr in nbs_lock_of_target.iter() {
                    let read_lock_nbr =
                        read_lock_internal_vec.get(*nbr).unwrap();
                    if read_lock_nbr.get_spin() == target_spin {
                        read_lock_internal_vec
                            .get(*nbr)
                            .unwrap()
                            .marked
                            .write()
                            .unwrap()
                            .set_pushed();
                        cluster_queue_signaler.send(*nbr).unwrap();
                    }
                }
            }
        } else {
            break 'cluster_loop
        }
    }
    deltas
}

pub fn thread_energy_worker(
    n_jobs: usize,
    this_thread_num: usize,
    div: usize,
    rem:usize,
    shared_internal_vector: Arc<RwLock<Vec<SpinNode>>>,
    ext_mag_field: f64,
    energy_psum_signaler: SignalContainer<f64>
) {
    let mut energy_psum = 0.;
    let mut mag_field_psum = 0.;
    let range = if this_thread_num == n_jobs - 1 {
        (this_thread_num * div)..((this_thread_num + 1) * div + rem)
    } else {
        (this_thread_num * div)..((this_thread_num + 1) * div)
    };
    if let Ok(read_lock_internal_vector) = shared_internal_vector.read() {
        for cur_node_index in range {
            // iterate through all nodes in this threads range
            if let Some(cur_node) =
                read_lock_internal_vector.get(cur_node_index)
            {
                let read_lock_neighbors_of_cur_node = cur_node
                    .neighbors
                    .read()
                    .expect("Couldnt read from internal_vector!");
                    let mut sum_of_cur_node_nbrs_spins = 0.;
                    
                    if ext_mag_field != 0. {
                        mag_field_psum -= ext_mag_field * cur_node.get_spin();
                    }
                    for nbr_index in read_lock_neighbors_of_cur_node.iter() {
                        // iterate over the current nodes neighbors
                        if let Some(cur_nodes_neighbor_ref) =
                            read_lock_internal_vector.get(*nbr_index)
                        {
                            sum_of_cur_node_nbrs_spins += cur_nodes_neighbor_ref.get_spin();
                        }
                    }
                    // multiply the sum of spins by the current spin value
                    energy_psum -= sum_of_cur_node_nbrs_spins * cur_node.get_spin();
            }
        }
    }
    energy_psum_signaler
        .send(energy_psum + mag_field_psum)
        .unwrap();
}