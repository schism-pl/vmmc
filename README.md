# vmmc
This code implements a fast (and probably correct) implementation of [Virtual Move Monte Carlo](https://pubs.aip.org/aip/jcp/article/127/15/154101/915022) in Rust.
The goal of this code was to provide a clean, fast, and well-tested implementation of the code that will be more amenable to outside experimentation.



### Running the code

First, make sure you have https://rustup.rs/ installed.

Then, you can build and run the code using:
```bash
cargo build --release
cargo run -- --help // see info about default parameters  
cargo run --release // runs 1000 sweeps of 1000 step attempts per particle. 
```

### Visualizing the simulation
The code produces an XYZ trajectory `trajectory.xyz` and a VMD script `vmd.tcl`
To visualize the simulation, just run:
```Bash
vmd trajectory.xyz -e vmd.tcl
```
