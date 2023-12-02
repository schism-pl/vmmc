use quickcheck::{Arbitrary, Gen};

pub mod cli;
mod consts;
pub mod io;
pub mod morphology;
pub mod particle;
pub mod patchy_discs;
pub mod position;
pub mod simbox;
pub mod stats;
pub mod vmmc;

#[derive(Clone)]
pub struct InputParams {
    pub num_particles: usize,
    pub interaction_energy: f64, // kBT
    pub patch_radius: f64,       // diameter of patch (in units of particle diameter)
    // density: f64,
    pub box_width: f64,
    pub box_height: f64,

    pub prob_translate: f64,
    pub max_translation: f64, // TODO: related to temperature somehow?
    pub max_rotation: f64,    // TODO: related to temperature somehow?
    pub num_sweeps: usize,
    pub steps_per_sweep: usize,
}

// TODO: need to eliminate patch radius from here
// TODO: eliminate density and replace with box length

impl Default for InputParams {
    fn default() -> Self {
        let num_particles = 1000;
        let interaction_energy = 8.0; // kBT
        let patch_radius = 0.1; // diameter of patch (in units of particle diameter)
        let box_width = 60.0;
        let box_height = 60.0;

        let prob_translate = 0.5;
        let max_translation = 0.15;
        let max_rotation = 0.2;
        let num_sweeps = 20;
        let steps_per_sweep = 1000;

        Self {
            num_particles,
            interaction_energy,
            patch_radius,
            box_width,
            box_height,

            prob_translate,
            max_translation,
            max_rotation,
            num_sweeps,
            steps_per_sweep,
        }
    }
}

// for arbitrary trait
fn usize_in_range(g: &mut Gen, min: usize, max: usize) -> usize {
    usize::arbitrary(g).min(min).max(max)
}

fn f64_in_range(g: &mut Gen, min: f64, max: f64) -> f64 {
    f64::arbitrary(g).min(min).max(max)
}

// for testing
impl Arbitrary for InputParams {
    fn arbitrary(g: &mut Gen) -> Self {
        let num_particles = usize_in_range(g, 0, 2500);
        let interaction_energy = f64_in_range(g, 0.1, 10.0); // kBT
        let patch_radius = f64_in_range(g, 0.01, 0.2); // diameter of patch (in units of particle diameter)
        let box_width = f64_in_range(g, 10.0, 100.0); // TODO: need to make sure its big enough to deal with all the particles
        let box_height = f64_in_range(g, 10.0, 100.0); // TODO: need to make sure its big enough to deal with all the particles
        let prob_translate = f64_in_range(g, 0.0, 1.0);
        let max_translation = f64_in_range(g, 0.0, 1.0);
        let max_rotation = f64_in_range(g, 0.0, 1.0);
        let num_sweeps = 10;
        let steps_per_sweep = 1000;

        Self {
            num_particles,
            interaction_energy,
            patch_radius,
            box_width,
            box_height,

            prob_translate,
            max_translation,
            max_rotation,
            num_sweeps,
            steps_per_sweep,
        }
    }
}
