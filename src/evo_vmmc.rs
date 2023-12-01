use std::f64::consts::PI;

use rand::rngs::SmallRng;

use crate::{
    patchy_discs::PatchyDiscParams,
    simbox::SimBox,
    tiling::TilingGraph,
    vmmc::{Vmmc, VmmcParams},
};

// fn overlapping_polygon_sum(vmmc: &Vmmc) -> usize {
//     // 1. generate edge list
//     let tiling_graph = TilingGraph::from_vmmc(vmmc);
//     // 2. get polgyon count
//     let polygons = tiling_graph.get_polygon_counts();
//     polygons.sum()
// }

pub struct L2GInputParams {
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
    survivors_per_generation: usize, // # of children post-pruning
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
        let survivors_per_generation = 1;

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
            survivors_per_generation,
        }
    }
}

pub enum FitnessFunc {
    Random,
    AvgEnergy,
    OverlappingPolygonSum,
}

// pub enum MutationFunc {
//     NumParticles,
//     Null,
// }

pub struct EvoVmmc {
    // active_sims: Vec<Vmmc>,
    fitness_func: FitnessFunc,
    params: L2GInputParams,
}

impl EvoVmmc {
    pub fn new(fitness_func: FitnessFunc) -> Self {
        let params = L2GInputParams::default();
        Self {
            fitness_func,
            params,
        }
    }

    // TODO: remove config? need seperate vmmc_from_config type functions
    fn fresh_vmmc(&self, rng: &mut SmallRng) -> Vmmc {
        let ip = &self.params;

        let base_length = ((ip.num_particles as f64 * PI) / (4.0 * ip.density)).sqrt();
        let box_dimensions = [base_length, base_length];
        // lower bound on cell length (max distance that cells can interact at)
        let cell_range = 1.0 + self.params.interaction_range;

        let cells_per_axis = (base_length / cell_range).floor();
        let cell_length = base_length / cells_per_axis;
        let cells_per_axis = [cells_per_axis as usize, cells_per_axis as usize];

        let cell_dimensions = [cell_length, cell_length];

        // No initial position, so we will use a randomized start position
        let simbox = SimBox::new_with_randomized_particles(
            box_dimensions,
            cells_per_axis,
            cell_dimensions,
            ip.num_particles,
            rng,
        );

        let pd_params =
            PatchyDiscParams::new(ip.num_patches, ip.interaction_energy, ip.interaction_range);

        let params = VmmcParams::new(
            ip.prob_translate,
            ip.max_trial_translation,
            ip.max_trial_rotation,
            ip.reference_radius,
        );

        Vmmc::new(simbox, params, pd_params)
    }

    fn run_sim(&self, vmmc: &mut Vmmc, rng: &mut SmallRng) -> f64 {
        let ip = &self.params;
        for _idx in 0..ip.num_sweeps {
            vmmc.step_n(ip.steps_per_sweep * ip.num_particles, rng);
        }
        vmmc.get_average_energy()
    }

    /// returns a number between 0.0 and 1.0
    pub fn fitness(&self, vmmc: &Vmmc) -> f64 {
        match self.fitness_func {
            FitnessFunc::AvgEnergy => 1.0 - 1.0 / vmmc.get_average_energy().abs(),
            FitnessFunc::Random => 0.5, // TODO: fix (need to pass through rng here?)
            FitnessFunc::OverlappingPolygonSum => panic!("Unimplemented"), //1.0 - 1.0 / overlapping_polygon_sum(vmmc),
        }
    }

    // TODO: return vec of vmmcs?
    fn step_generation(&mut self, rng: &mut SmallRng) -> Vec<Vmmc> {
        let mut children = Vec::new();
        for jdx in 0..self.params.children_per_generation {
            let mut vmmc = self.fresh_vmmc(rng);
            let avg_energy = self.run_sim(&mut vmmc, rng);
            println!(
                "average energy of child-sim {:?} = {:?} kBT\n",
                jdx, avg_energy
            );
            children.push(vmmc);
        }
        children
    }

    // TODO: how to derive new generation from old?
    pub fn step_generation_n(&mut self, n: usize, rng: &mut SmallRng) {
        for idx in 0..n {
            println!("Starting generation {:?}", idx);
            let children = self.step_generation(rng);
            let survivors = self.prune(children);
            let survivor_ids: Vec<usize> = survivors.iter().map(|v| v.0).collect();
            println!("{:?} survive", survivor_ids);
        }
    }

    pub fn step_all(&mut self, rng: &mut SmallRng) {
        self.step_generation_n(self.params.num_generations, rng)
    }

    fn get_index_of_least_fit(&self, values: &[(usize, Vmmc)]) -> usize {
        let mut lowest = f64::MIN;
        let mut index = 0;
        for i in 0..values.len() {
            if self.fitness(&values[i].1) < lowest {
                lowest = self.fitness(&values[i].1);
                index = i;
            }
        }
        index
    }

    // generate children from survivors of previously generations
    pub fn spawn_children(&self, vmmcs: Vec<Vmmc>) -> Vec<Vmmc> {
        panic!("TODO")
    }

    // Note: this function is implemented assuming that computing fitness is cheap
    // pretty easy to optimize if that isn't the case
    pub fn prune(&self, vmmcs: Vec<Vmmc>) -> Vec<(usize, Vmmc)> {
        let n = self.params.survivors_per_generation;
        let mut survivors = Vec::new();

        for (idx, vmmc) in vmmcs.into_iter().enumerate() {
            // the first n children get to start as survivors
            if idx < n {
                survivors.push((idx, vmmc));
                continue;
            }
            // attempt to replace one of the survivors
            let fitness = self.fitness(&vmmc);
            let index = self.get_index_of_least_fit(&survivors);
            if fitness > self.fitness(&survivors[index].1) {
                survivors[index] = (idx, vmmc);
            }
        }
        survivors
    }
}
