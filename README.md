# Ising_model_rust
The spin glass ising model programmed in rust

# Running
The easiest way to run this is to run ```cargo run --release && python -O plot.py```. My build target was UNIX, so if you are on windows you might need to build this yourself.

# Install RUST
Please follow the instructions here on the RUST-lang website https://www.rust-lang.org/tools/install.

# Whats working and whats not
The metropolis-hastings algorithm works, the wolff algorithm doesnt for reasons I am not entirely sure of, but it appears to be a lot of parameter tweaking, as that is what led to the metropolis working.
