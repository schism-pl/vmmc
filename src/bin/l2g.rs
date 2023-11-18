use clap::Parser;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;
use vmmc::cli::VmmcConfig;
use vmmc::particle::IsParticle;
use vmmc::position::Position;
use vmmc::{
    io::{read_xyz_snapshot, write_tcl, XYZWriter},
    particle::Particle,
    patchy_discs::PatchyDiscParams,
    position::Orientation,
    simbox::SimBox,
    vmmc::{Vmmc, VmmcParams},
};

// grab first frame from xyz and load particles
fn particles_from_xyz(path: &str) -> Vec<Particle> {
    let mut particles = Vec::new();
    let (positions, orientations) = read_xyz_snapshot(path);

    for idx in 0..positions.len() {
        let pos = positions[idx];
        let or = orientations[idx];
        let particle = Particle::new(idx as u16, pos, or);
        particles.push(particle);
    }
    particles
}

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)


// TODO: builder pattern
fn main() {
    env_logger::init();
    // Get commandline arguments
    let config = VmmcConfig::parse();
    // Get default params

    // Seed the rng
    let seed = config.seed();
    println!("Using seed = {:?}", seed);
    let mut rng = SmallRng::seed_from_u64(seed);

    // TODO: implement fitness function class
    // TODO: implement pruning functionality
    // TODO: implement Metasim class
    // TODO: parallelize children
    // TODO: randomize starts?
    // TODO: modify something meaningful each run (chemical potential / interaction range?)
    // for idx in 0..ip.num_generations {
    //     println!("Starting generation {:?}", idx);
    //     for jdx in 0..ip.children_per_generation {
    //         println!("Simulating child {:?}-{:?}", idx, jdx);
    //         let mut vmmc = vmmc_from_config(&ip, &mut rng);
    //         let avg_energy = run_sim(&mut vmmc, &ip, &mut rng);
    //         println!("average energy of sim {:?}-{:?} = {:?} kBT\n", idx, jdx, avg_energy);
    //     }
    // }
}
