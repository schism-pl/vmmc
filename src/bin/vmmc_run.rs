use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;
use vmmc::{
    io::{read_xyz_snapshot, write_tcl, XYZReader, XYZWriter},
    particle::Particle,
    patchy_discs::PatchyDiscParams,
    position::Orientation,
    simbox::SimBox,
    vmmc::{Vmmc, VmmcParams},
};

// Try to formalize correctness conditions of box / particles?

// TODO: probably need to check for overlapping and particles and do some other checks
fn randomized_particles(simbox: &SimBox, n: usize, rng: &mut SmallRng) -> Vec<Particle> {
    let mut particles = Vec::new();

    for idx in 0..n {
        let pos = simbox.random_pos(rng);
        let or = Orientation::unit_vector(rng);
        let particle = Particle::new(idx as u16, pos, or);
        particles.push(particle);
    }
    particles
}

// grab first frame from xyz and load particles
fn particles_from_xyz(path: &str) -> Vec<Particle> {
    let mut particles = Vec::new();
    // let mut rng = rand::thread_rng();
    // let mut rdr = XYZReader::new(path);
    let (positions, orientations) = read_xyz_snapshot(path);

    // let mut or_rdr = XYZReader::new(orientations);
    // let orientations = or_rdr.read_xyz_frame();

    for idx in 0..positions.len() {
        // let or = Orientation::unit_vector(&mut rng);
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

    let num_particles = 1000;
    let interaction_energy = 8.0; // kBT
    let interaction_range = 0.1; // diameter of patch (in units of particle diameter)
    let density = 0.2;
    let max_interactions = 3;

    let base_length = ((num_particles as f64 * PI) / (4.0 * density)).sqrt();
    let box_dimensions = [base_length, base_length];
    // approximate cell length
    let cell_range = 1.0 + 0.5 * interaction_range;

    let cells_per_axis = (base_length / cell_range).floor();
    let cell_length = base_length / cells_per_axis;
    let cells_per_axis = [cells_per_axis as usize, cells_per_axis as usize];

    // let cells_per_axis = (box_dimensions / cell_length).floor()

    let cell_dimensions = [cell_length, cell_length];

    let sqd_cutoff_distance = interaction_range * interaction_range;
    // let particles = randomized_particles(&simbox, num_particles);
    let particles = particles_from_xyz("snapshot.xyz");
    let simbox = SimBox::new(box_dimensions, cells_per_axis, cell_dimensions, &particles);

    let pd_params =
        PatchyDiscParams::new(max_interactions, interaction_energy, sqd_cutoff_distance);

    let prob_translate = 0.5;
    let max_trial_translation = 0.15;
    let max_trial_rotation = 0.2;
    let reference_radius = 0.5;

    let params = VmmcParams::new(
        prob_translate,
        max_trial_translation,
        max_trial_rotation,
        reference_radius,
    );

    let mut writer = XYZWriter::new("trajectory.xyz");

    let mut vmmc = Vmmc::new(particles, simbox, params, pd_params);

    println!("Checking initial configuration");
    debug_assert!(vmmc.well_formed());
    println!("Initial configuration ok");

    println!("# of particles: {:?}", vmmc.particles().len());
    println!("Box dimensions: {:?}", box_dimensions);
    println!("Cell dimensions: {:?}", cell_dimensions);
    println!("Initial average energy: {:?}", vmmc.get_average_energy());

    //let mut rng = rand::thread_rng();
    let mut rng = SmallRng::seed_from_u64(1337);

    for idx in 0..10 {
        vmmc.step_n(1000 * num_particles, &mut rng);
        // writer.write_xyz_frame(&vmmc);
        println!(
            "Step {:?}: average particle energy = {:?}",
            (idx + 1) * 1000 * num_particles,
            vmmc.get_average_energy()
        );
    }
    // writer.write_xyz_frame(&vmmc);

    // vmmc.step_n(1000 * 1000);

    // write_tcl(&vmmc, "vmd.tcl");

    // write_svg(&vmmc, "box0.svg");
    // println!("stats: {:?}", vmmc.st);
    // write_svg(&vmmc, "box1.svg");
}
