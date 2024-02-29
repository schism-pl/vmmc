# VMMC
This code implements a fast (and probably correct) implementation of [Virtual Move Monte Carlo](https://pubs.aip.org/aip/jcp/article/127/15/154101/915022) in Rust.
The goal of this code was to provide a clean, fast, and well-tested implementation of the code that will be more amenable to outside experimentation.

# Basic Instructions

### Building

First, make sure you have https://rustup.rs/ installed.

Then, you can build and run the code using:
```bash
cargo build --release
```

### Running
```
// runs 1000 sweeps of 1000 step attempts per particle.
cargo run --release
// run, but do so by reading `config.toml` for config info and writing files to the directory `out`
cargo run --release -- --input=config.toml output=out
```

# Inputs
The tool reads parameters from a config file in the [Toml](https://toml.io/en/) format.
Below is a description of the Toml file, examples of valid config files, and an exhaustive list of the parameters that the tool uses.

### Config.toml

Below is an example of a VMMC config file.
The first 4 parameters (seed, initial_particles, box_width, and box_height) are all self explanatory (although they are described in detail in the section below).
The protocol argument describes 

```TOML
seed = 2037276277152067866
initial_particles = 500
box_width = 50.0
box_height = 50.0

[protocol]
chemical_potential_eq = "(0.05 * t)"
interaction_energy_eq = "(5 + (0.05 * (t^2)))"
num_megasteps = 20

[[shapes]]

[[shapes.patches]]
radius = 0.05
theta = 0.0
chemtype = 0

[[shapes.patches]]
radius = 0.05
theta = 90.0
chemtype = 0

[[shapes.patches]]
radius = 0.05
theta = 180.0
chemtype = 0

[[shapes.patches]]
radius = 0.05
theta = 270.0
chemtype = 0

```



### Parameters
Below I describe the parameters that the tool uses. 
Each parameter has a name, type, a short description, its default value, as well as the expected value range for it.
Any value not specified will take its default value.
The tool should reject parameters outside of the expected value range, but if it doesn't file and issue or ping Evan (the same if you think the expected value range is wrong / unreasonable).

##### seed
Description: Random seed for the simulation. 
Type: `i64`  
Expected Value Range: `any`  
Default Value: `random`  

##### initial_particles
Description: Number of initial particles. Placed uniformly at random within the simulation box.  
Type: `usize`  
Expected Value Range: `0 <= initial_particles <= 2000`  
Default Value: `500`  

##### box_width
Description: Width of simulation box in units of particle diameter  
Type: `f64`  
Expected Value Range: `10.0 <= box_width <= 200.0`  
Default Value: `50.0`  

##### box_height
Description: Height of simulation box in units of particle diameter  
Type: `f64`  
Expected Value Range: `10.0 <= box_height <= 200.0`  
Default Value: `50.0`  




##### protocol
Description: Protocol for simulation. For each sweep, the protocol contains 1 ProtocolStep. 
Each ProtocolStep contains a chemical potential (mu) in KbT (?) and interaction_energy (epsilon) in KbT (?) for that timestep.

Type: `[ProtocolStep(mu: f64, epsilon: f64)]`  
Expected Value Range:  
```
protocol.len() == num_sweeps 
foreach ProtocolStep(mu, epsilon). -20.0 <= mu <= 20.0 && 0.01 <= epsilon <= 20.0
```



The default protocol is one that keeps mu fixed at 0.0 and epsilon at 10.0. 


##### shapes
Type: 
```
shapes:: [Shape]
Shape:: [Patch]
Patch:: (color: , theta: f64, radius: f64)
```

Expected Value Range:
```
foreach shape. in shapes:
  3 <= shape.patches.len() <= 6
  foreach patch. in shape.patches:
     color = any
     0 <= theta <= 360.0 
     0.01 <= radius <= 0.25 
```

Description:






# Outputs


### Visualizing the simulation
The code produces an XYZ trajectory `trajectory.xyz` and a VMD script `vmd.tcl`
To visualize the simulation, just run:
```Bash
vmd trajectory.xyz -e vmd.tcl
```


