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
//"&" designates a background task
cargo run --release -- --input=config.toml --output-dir=out &
```

# Inputs
The tool reads parameters from a config file in the [Toml](https://toml.io/en/) format.
Below is a description of the Toml file, examples of valid config files, and an exhaustive list of the parameters that the tool uses.

### Config.toml

Below is an example of a VMMC config file.
The first 4 parameters (`seed`, `initial_particles`, `box_width`, and `box_height`) are all self explanatory (although they are described in detail in the section below).

The `protocol` argument describes the synthesis protocol that we are simulating. 
It has 3 fields, `chemical_potential_eq` and `interaction_energy_eq` are two equations that describe the chemical potential and interaction energy over the length of the simulation.
Currently, these equations can have at max only 1 variable, and this variable is the time variable, counted in Megasteps (so t=1 is equivalent to 10e6 simulation steps).
The `num_megasteps` field describes the length of the simulation. For the below config, the simulation would end after 20\*10e6 steps.

The `shapes` argument describes the morphologies of the particles in the simulation. Each `shapes` entry describes one morphology (in the below example, there are two, a uniform 4 patch and a uniform 3 patch particle). If there is more than 1 morphology, the simulation uses a uniform mixture of the particles (non-uniform mixtures are not implemented, but can be added easily).

![Screenshot from 2024-03-11 11-14-06](https://github.com/schism-pl/vmmc/assets/82984409/c0228ed3-825c-4856-bdb1-c099f960f126)

Each patch consists of a `radius` and `theta`. Where `radius` is defined as half the length of a chord subtended by the disc, and `theta` is defined as the location for the center of a patch. `chemtype` describes the chemical selectivity of the patch, i.e., a `chemtype=0` patch can bond with a `chemtype=0` patch but not a `chemtype=1` patch.

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

[[shapes]]

[[shapes.patches]]
radius = 0.05
theta = 0.0
chemtype = 0

[[shapes.patches]]
radius = 0.05
theta = 120.0
chemtype = 0

[[shapes.patches]]
radius = 0.05
theta = 240.0
chemtype = 0

```



### Parameters
Below I describe the parameters that the tool uses. 
Each parameter has a name, type, a short description, its default value, as well as the expected value range for it.
Any value not specified will take its default value.
The tool should reject parameters outside of the expected value range. If it doesn't, file an issue or message Evan. Let us know if you think the expected value range is wrong / unreasonable for any of the parameters.

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
Default Value: `75.0`  

##### box_height
Description: Height of simulation box in units of particle diameter  
Type: `f64`  
Expected Value Range: `10.0 <= box_height <= 200.0`  
Default Value: `75.0`  




##### protocol
Description: Synthesis protocol for the simulation. See above for more in-depth description
Type: `Protocol` 

Expected Value Range:  
```
forall. 0 <= t <= num_megasteps: 0.01 <= interaction_energy_eq(t) <= 20.0
forall. 0 <= t <= num_megasteps: -20.0 <= chemical_potential_eq(t) <= 20.0
```

Default value: The default protocol is one that keeps chemical_potential fixed at 0.0 and interaction_energy fixed at 8.0 for a total of 1000 megasteps. 


##### shapes
Description:

Type: `Shapes`

Expected Value Range:
```
foreach shape. in shapes:
  3 <= shape.patches.len() <= 6
  foreach patch. in shape.patches:
     color = any
     0 <= theta <= 360.0 
     0.01 <= radius <= 0.25 
```

Default value: The default `shapes` is a a regular 4-patch particle with `radius=0.05` and all `chemtype=0`.




# Outputs


### Visualizing the simulation
The code produces 4 files: an XYZ trajectory `trajectory.xyz`, a VMD script `vmd.tcl`, a copy of the config `config.toml` (useful for reproducing runs), and an image `geometry.png` that shows the patch interactions for the end state of the simulation.
To visualize the simulation, just run:
```Bash
vmd trajectory.xyz -e vmd.tcl
```


