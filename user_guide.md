# VMMC
This code implements a fast (and probably correct) implementation of [Virtual Move Monte Carlo](https://pubs.aip.org/aip/jcp/article/127/15/154101/915022) in Rust.
The goal of this code was to provide a clean, fast, and well-tested implementation of the code that will be more amenable to outside experimentation.

# Basic Instructions

### Building

First, make sure you have https://rustup.rs/ installed.

Then, you can build and run the code using:
```bash
cargo build --release
cargo run -- --help // see info about default parameters  
cargo run --release // runs 1000 sweeps of 1000 step attempts per particle. 
```

### Running


# Inputs
The tool reads parameters from a config file in the [Toml](https://toml.io/en/) format.
Below is a description of the Toml file, examples of valid config files, and an exhaustive list of the parameters that the tool uses.

### Config.toml

description of types

Short example
Link to defaults

how to specify which file to use

### Parameters
Below I describe the parameters that the tool uses. 
Each parameter has a name, type, a short description, its default value, as well as the expected value range for it.
Any value not specified will take its default value.
The tool should reject parameters outside of the expected value range, but if it doesn't file and issue or ping Evan (the same if you think the expected value range is wrong / unreasonable).

##### box_width
Type: `f64`  
Expected Value Range: `10.0 <= box_width <= 200.0`  
Description: Width of simulation box in units of particle diameter  

##### box_height
Type: `f64`  
Expected Value Range: `10.0 <= box_height <= 200.0`  
Description: Height of simulation box in units of particle diameter  

##### num_particles
Type: `usize`  
Expected Value Range: `0 <= num_particles <= 2000`  
Description: Number of particles. Placed uniformly at random within the specified box.  


##### prob_translate
Type: `f64`  
Expected Value Range:  `0.0 <= prob_translate <= 1.0`  
Description: Probability of a particle move being a translate (as opposed to a rotation)   

##### max_translate
Type: `f64`  
Expected Value Range:  `0.0 <= max_translate <= 1.0`  
Description: Maximum distance of a translation move, in units of particle diameter 

##### max_rotation
Type: `f64`  
Expected Value Range: `0.0 <= max_rotation <= 1.0`  
Description: Maximum rotation angle of a rotation move, in units of radians 


##### num_sweeps
Type: `usize`  
Expected Value Range: `0 >= num_sweeps`   
Description: Monte-Carlo steps for the simulation to take. 1 sweep = 10e6 moves.  

##### protocol
Type: [ProtocolStep]
Expected Value Range: `protocol.len() == num_sweeps. Foreach ProtocolStep(interaction_energy: f64, chemical_potential: f64). 0.01 <= interaction_energy <= 20.0 && -20.0 <= chemical_potential <=  `
Description: Protocol for simulation. For each sweep, the protocol contains 1 ProtocolStep. 
Each ProtocolStep contains an interaction_energy and a chemical potential.
TODO: finish


##### shapes
Type: [Shape]
Expected Value Range:
Description:






# Outputs


### Visualizing the simulation
The code produces an XYZ trajectory `trajectory.xyz` and a VMD script `vmd.tcl`
To visualize the simulation, just run:
```Bash
vmd trajectory.xyz -e vmd.tcl
```


