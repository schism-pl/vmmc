use std::f64::consts::PI;

use rand::rngs::SmallRng;

use crate::{vmmc::{Vmmc, VmmcParams}, patchy_discs::PatchyDiscParams, simbox::SimBox, position::{Orientation, Position}, particle::{Particle, IsParticle}};

// check if a particle overlaps any other particles
fn overlaps(simbox: &SimBox, pos: &Position, particles: &[Particle]) -> bool {
    for other in particles.iter() {
        // dist < 1.0 (hard sphere radius)
        if simbox.sep_in_box(*pos, other.pos()).norm() < 1.0 {
            return true;
        }
    }
    false
}

// TODO: dedup with other impl of randomized_particles and overlaps
// TODO: probably need to check for overlapping and particles and do some other checks
fn randomized_particles(simbox: &SimBox, n: usize, rng: &mut SmallRng) -> Vec<Particle> {
    let mut particles = Vec::new();

    for idx in 0..n {
        let mut pos = simbox.random_pos(rng);
        while overlaps(simbox, &pos, &particles) {
            pos = simbox.random_pos(rng);
        }
        let or = Orientation::unit_vector(rng);
        let particle = Particle::new(idx as u16, pos, or);
        particles.push(particle);
    }
    particles
}

fn vmmc_from_config(ip: &L2GInputParams, rng: &mut SmallRng) -> Vmmc {

    let base_length = ((ip.num_particles as f64 * PI) / (4.0 * ip.density)).sqrt();
    let box_dimensions = [base_length, base_length];
    // lower bound on cell length (max distance that cells can interact at)
    let cell_range = 1.0 + ip.interaction_range;

    let cells_per_axis = (base_length / cell_range).floor();
    let cell_length = base_length / cells_per_axis;
    let cells_per_axis = [cells_per_axis as usize, cells_per_axis as usize];

    let cell_dimensions = [cell_length, cell_length];

    // let particles = if !config.start_frame().is_empty(){
    //     particles_from_xyz(config.start_frame())
    // }
    // else {
    //     randomized_particles(simbox, n, &mut rng)
    // };
    // let particles = particles_from_xyz(config.start_frame());
    let particles = randomized_particles(simbox, ip.num_particles, &mut rng);

    let simbox = SimBox::new(box_dimensions, cells_per_axis, cell_dimensions, &particles);

    let pd_params =
        PatchyDiscParams::new(ip.num_patches, ip.interaction_energy, ip.interaction_range);

    let params = VmmcParams::new(
        ip.prob_translate,
        ip.max_trial_translation,
        ip.max_trial_rotation,
        ip.reference_radius,
    );

    // println!("Box dimensions: {:?}", box_dimensions);
    // println!("Cell dimensions: {:?}", cell_dimensions);
    Vmmc::new(particles, simbox, params, pd_params)
}


fn run_sim(vmmc: &mut Vmmc, ip: &L2GInputParams, rng: &mut SmallRng) -> f64 {
    for idx in 0..ip.num_sweeps {
        vmmc.step_n(ip.steps_per_sweep * ip.num_particles, rng);
        // println!(
        //     "Step {:?}: average particle energy = {:?}",
        //     (idx + 1) * 1000 * ip.num_particles,
        //     vmmc.get_average_energy()
        // );
    };
    vmmc.get_average_energy()
}


struct L2GInputParams {
    num_particles: usize,
    interaction_energy: f64, // kBT
    interaction_range: f64,  // diameter of patch (in units of particle diameter)
    density: f64,
    num_patches: usize,

    prob_translate: f64,
    max_trial_translation: f64,
    max_trial_rotation: f64,
    reference_radius: f64,
    num_sweeps: usize,
    steps_per_sweep: usize,

    // L2G parameters
    num_generations: usize,
    children_per_generation: usize,
}

impl Default for L2GInputParams {
    fn default() -> Self {
        let num_particles = 500; // starts at 500 and changes over time
        let interaction_energy = 8.0; // The GA will change this. For now it can stay at 8.0
        let interaction_range = 0.05; // checked
        let density = 0.2; // doesn't use density? we'll keep it anyways lol
        let num_patches = 3;

        let prob_translate = 0.5; // checked
        let max_trial_translation = 0.15; // checked
        let max_trial_rotation = 0.06; // checked
        let reference_radius = 0.5; // checked
        let num_sweeps = 5;
        let steps_per_sweep = 1000;

        let num_generations = 3;
        let children_per_generation = 3;

        Self {
            num_particles,
            interaction_energy,
            interaction_range,
            density,
            num_patches,

            prob_translate,
            max_trial_translation,
            max_trial_rotation,
            reference_radius,
            num_sweeps,
            steps_per_sweep,

            num_generations,
            children_per_generation,
        }
    }
}


pub enum FitnessFunc {
    Random,
    AvgEnergy,
}

pub struct EvoVmmcParams {
    children_per_generation: usize,
}

pub struct EvoVmmc {
    active_sims: Vec<Vmmc>,
    fitness_func: FitnessFunc,
    params: EvoVmmcParams,
}


impl EvoVmmc {
    pub fn new() -> Self {
        panic!("TODO")
    }

    /// returns a number between 0.0 and 1.0
    pub fn fitness(&self, vmmc: &Vmmc) -> f64 {
        panic!("TODO")
    }

    fn step_generation(&mut self, rng: &mut SmallRng){
        // println!("Starting generation {:?}", idx);
        for jdx in 0..self.params.children_per_generation {
            // println!("Simulating child {:?}-{:?}", idx, child_id);
            let ip = L2GInputParams::default();
            let mut vmmc = vmmc_from_config(&ip, rng);
            let avg_energy = run_sim(&mut vmmc, &ip, rng);
            // println!("average energy of sim {:?}-{:?} = {:?} kBT\n", jdx, avg_energy);
        }
    }


    pub fn step_generation_n(&mut self, n: usize){
        panic!("TODO")
    }


    pub fn prune(){}
}
