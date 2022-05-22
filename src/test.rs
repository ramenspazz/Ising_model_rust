// push the random node neighbors indicies to the work queue
if let Ok(read_lock_internal_vector) = self.internal_lattice.internal_vector.read() {
    if let Some(target_node) = read_lock_internal_vector.get(target_index) {
        let read_lock_nbrs = target_node.neighbors.read().unwrap();
        for i in 0..read_lock_nbrs.len() {
            // SAFETY: bounds checked in Driver::new, garaunteed valid return.
            if let Some(nbr_index) = read_lock_nbrs.get(i) {
                let nbr_spin = self.get_spin_at(*nbr_index);
                if nbr_spin == target_spin {
                    // if the spin is the same as the randomly picked node, add it to the
                    // queue.
                    self.cluster_queue_signaler
                        .send(*nbr_index)
                        .expect("Couldnt send message through cluster_queue_signaler!");
                }
            }
        }
    }
}



// spin up threads for generating a cluster flip
for _ in 0..self.num_threads {
    self.cluster_go_stop_signaler
        .send((SignalType::SigStart, balance_condition, target_spin))
        .unwrap();
}
// waiting for threads to finish
while !self.cluster_done_signaler.is_full() {
    stdth::sleep(time::Duration::from_micros(1));
}
// Clear queue
for _ in 0..self.num_threads {
    // clear the queue
    _ = self.cluster_done_signaler.recv().unwrap();
}

// now we flip the resulting cluster and calculate the change in energy and magnitization
delta_mag += -1. * target_spin * self.flip_index_vector.read().unwrap().len() as f64;
// get a write lock to flip_index_vecor
unsafe {
    'flip: loop {
        let mut write_lock_flip_index_vector = self.flip_index_vector.write().unwrap();
        // for item in write_lock_flip_index_vector.iter() {
        //     println!("flip nodes spin is {} and target was {}",
        //         self.internal_lattice.internal_vector.read().unwrap().get_unchecked(*item).get_spin(),
        //         &target_spin);
        // }
        // get_input_as_usize(Some("Enter to continue"));
        let opt_cur_flip_node_index= write_lock_flip_index_vector.pop();
        if let Some(cur_flip_node_index) = opt_cur_flip_node_index {
            if self.internal_lattice.internal_vector.read().unwrap().get_unchecked(cur_flip_node_index).get_spin() != target_spin {
                continue 'flip
            }
            // assert_eq!(self.internal_lattice.internal_vector.read().unwrap().get_unchecked(cur_flip_node_index).get_spin(), target_spin);
            let mut write_lock_internal_vector =
            self.internal_lattice.internal_vector.write().unwrap();
            // preform the accepted spin flip
            write_lock_internal_vector
            .get_unchecked_mut(cur_flip_node_index)
            .flip_spin();
        } else { break 'flip; }
        drop(write_lock_flip_index_vector);
        // flip already happened, so take the negative of the change in energy
        delta_energy += self.calculate_energy_change_of_suggested_flip(opt_cur_flip_node_index.unwrap());
    }
    // unmark all touched nodes, I cant really think of a better way to do this at the moment.
    let mut write_lock_touched_index_vector = self.touched_index_vec.write().unwrap();
    if let Some(touched_index) = write_lock_touched_index_vector.pop() {
        let mut write_lock_internal_vector =
            self.internal_lattice.internal_vector.write().unwrap();
        write_lock_internal_vector
            .get_unchecked_mut(touched_index)
            .marked
                .write()
                .unwrap()
                .set_unmark();
    }
}