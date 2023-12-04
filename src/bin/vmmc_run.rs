use clap::Parser;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use vmmc::cli::VmmcConfig;
use vmmc::io::write_geometry_png;
use vmmc::morphology::Morphology;
use vmmc::position::DimVec;
use vmmc::InputParams;
use vmmc::{
    io::{read_xyz_snapshot, write_tcl, XYZWriter},
    particle::Particle,
    simbox::SimBox,
    vmmc::{Vmmc, VmmcParams},
};

// grab first frame from xyz and load particles
// TODO: dedup with other particles_from_xyz
// TODO: need to specify shape id, and specify morphologies somehow
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

fn vmmc_from_config(config: &VmmcConfig, ip: &InputParams, rng: &mut SmallRng) -> Vmmc {
    let box_dimensions = DimVec::new([ip.box_width, ip.box_height]);

    let max_interaction_range = 1.0 + ip.patch_radius;

    let shapes = vec![
        // Morphology::regular_3patch(ip.patch_radius),
        Morphology::regular_3patch(ip.patch_radius),
    ];

    let simbox = if !config.start_frame().is_empty() {
        // We have an initial position, so just use that
        let particles = particles_from_xyz_nomix(config.start_frame());
        let cells_x_axis = (box_dimensions.x() / max_interaction_range).floor();
        let cells_y_axis = (box_dimensions.y() / max_interaction_range).floor();

        let cell_width = box_dimensions.x() / cells_x_axis;
        assert!(cell_width >= max_interaction_range);

        let cells_per_axis = [cells_x_axis as usize, cells_y_axis as usize];

        // cells are always square
        let cell_dimensions = DimVec::new([cell_width, cell_width]);

        SimBox::new(
            box_dimensions,
            cells_per_axis,
            cell_dimensions,
            particles,
            shapes,
        )
    } else {
        // No initial position, so we will use a randomized start position
        SimBox::new_with_randomized_particles(
            box_dimensions,
            max_interaction_range,
            ip.num_particles,
            shapes,
            rng,
        )
    };

    let params = VmmcParams::new(ip.prob_translate, ip.max_translation, ip.max_rotation);

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
    // println!("Box dimensions: {:?}", box_dimensions);
    // println!("Cell per axis: {:?}", cells_per_axis);
    // println!("Cell dimensions: {:?}", cell_dimensions);
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
        vmmc.step_n(ip.num_particles * ip.num_particles, &mut rng);
        println!(
            "Step {:?}: bonds per particle = {:?}",
            (idx + 1) * 1000 * ip.num_particles,
            -vmmc.get_average_energy() / ip.interaction_energy
        );
    }
    // write the final frame
    writer.write_xyz_frame(&vmmc);

    // Write script to generate animation
    write_tcl(&vmmc, config.vmd_output());
    write_geometry_png(&vmmc, "geometry.png"); // TODO: configurable path
}
