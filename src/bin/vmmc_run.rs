use clap::Parser;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::f64::consts::PI;
use vmmc::cli::VmmcConfig;
use vmmc::morphology::Morphology;
use vmmc::{
    io::{read_xyz_snapshot, write_tcl, XYZWriter},
    particle::Particle,
    simbox::SimBox,
    vmmc::{Vmmc, VmmcParams},
};

// grab first frame from xyz and load particles
// TODO: dedup with other particles_from_xyz
fn particles_from_xyz_nomix(path: &str) -> Vec<Particle> {
    let mut particles = Vec::new();
    let (positions, orientations) = read_xyz_snapshot(path);

    for idx in 0..positions.len() {
        let pos = positions[idx];
        let or = orientations[idx];
        let particle = Particle::new(idx as u16, pos, or, 0);
        particles.push(particle);
    }
    particles
}

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)

struct InputParams {
    num_particles: usize,
    interaction_energy: f64, // kBT
    patch_radius: f64,       // diameter of patch (in units of particle diameter)
    density: f64,

    prob_translate: f64,
    max_translation: f64,
    max_rotation: f64,
    reference_radius: f64,
    num_sweeps: usize,
    steps_per_sweep: usize,
}

impl Default for InputParams {
    fn default() -> Self {
        let num_particles = 1000;
        let interaction_energy = 8.0; // kBT
        let patch_radius = 0.1; // diameter of patch (in units of particle diameter)
        let density = 0.2;

        let prob_translate = 0.5;
        let max_translation = 0.15;
        let max_rotation = 0.2;
        let reference_radius = 0.5;
        let num_sweeps = 1000;
        let steps_per_sweep = 1000;

        Self {
            num_particles,
            interaction_energy,
            patch_radius,
            density,

            prob_translate,
            max_translation,
            max_rotation,
            reference_radius,
            num_sweeps,
            steps_per_sweep,
        }
    }
}

fn vmmc_from_config(config: &VmmcConfig, ip: &InputParams, rng: &mut SmallRng) -> Vmmc {
    let base_length = ((ip.num_particles as f64 * PI) / (4.0 * ip.density)).sqrt();
    let box_dimensions = [base_length, base_length];
    // lower bound on cell length (max distance that cells can interact at)
    let cell_width = 1.0 + ip.patch_radius;

    let cells_per_axis = (base_length / cell_width).floor();
    let cell_length = base_length / cells_per_axis;
    let cells_per_axis = [cells_per_axis as usize, cells_per_axis as usize];

    let cell_dimensions = [cell_length, cell_length];

    // only using 1 shape for now
    let shape = Morphology::regular_3patch(ip.patch_radius);

    let simbox = if !config.start_frame().is_empty() {
        // We have an initial position, so just use that
        let particles = particles_from_xyz_nomix(config.start_frame());
        SimBox::new(
            box_dimensions,
            cells_per_axis,
            cell_dimensions,
            particles,
            vec![shape],
        )
    } else {
        // No initial position, so we will use a randomized start position
        SimBox::new_with_randomized_particles_nomix(
            box_dimensions,
            cells_per_axis,
            cell_dimensions,
            ip.num_particles,
            shape,
            rng,
        )
    };

    // let pd_params = PatchyDiscParams::new(ip.num_patches, ip.interaction_energy, ip.patch_radius);

    let params = VmmcParams::new(
        ip.prob_translate,
        ip.max_translation,
        ip.max_rotation,
        ip.reference_radius,
    );

    println!("Box dimensions: {:?}", box_dimensions);
    println!("Cell dimensions: {:?}", cell_dimensions);
    Vmmc::new(simbox, params, ip.interaction_energy)
}

// TODO: builder pattern
fn main() {
    env_logger::init();
    // Get commandline arguments
    let config = VmmcConfig::parse();
    // Get default params
    let ip = InputParams::default();

    // Seed the rng
    let seed = config.seed();
    println!("Using seed = {:?}", seed);
    let mut rng = SmallRng::seed_from_u64(seed);

    // Generate the simulator
    let mut vmmc = vmmc_from_config(&config, &ip, &mut rng);

    // Init I/O
    let mut writer = XYZWriter::new(config.xyz_output());

    println!("Checking initial configuration");
    debug_assert!(vmmc.well_formed());
    println!("Initial configuration ok");

    println!("# of particles: {:?}", vmmc.particles().len());
    println!("Initial average energy: {:?}", vmmc.get_average_energy());
    println!("------------------------------");

    for idx in 0..ip.num_sweeps {
        writer.write_xyz_frame(&vmmc);
        vmmc.step_n(ip.steps_per_sweep * ip.num_particles, &mut rng);
        println!(
            "Step {:?}: average particle energy = {:?}",
            (idx + 1) * 1000 * ip.num_particles,
            vmmc.get_average_energy()
        );
    }
    // write the final frame
    writer.write_xyz_frame(&vmmc);

    // Write script to generate animation
    write_tcl(&vmmc, config.vmd_output());
}
