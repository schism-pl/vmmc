use clap::Parser;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use vmmc::cli::VmmcConfig;
use vmmc::io::write_geometry_png;
use vmmc::morphology::Morphology;
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
    let box_dimensions = [ip.box_width, ip.box_height];
    // lower bound on cell length (max distance that cells can interact at)
    let cell_width_min = 1.0 + ip.patch_radius;

    let cells_x_axis = (ip.box_width / cell_width_min).floor();
    let cells_y_axis = (ip.box_height / cell_width_min).floor();

    let cell_width = ip.box_width / cells_x_axis;
    assert!(cell_width >= cell_width_min);

    let cells_per_axis = [cells_x_axis as usize, cells_y_axis as usize];

    // cells are always square
    let cell_dimensions = [cell_width, cell_width];

    println!("Box dimensions: {:?}", box_dimensions);
    println!("Cell per axis: {:?}", cells_per_axis);
    println!("Cell dimensions: {:?}", cell_dimensions);

    let shapes = vec![
        Morphology::regular_3patch(ip.patch_radius),
        Morphology::regular_4patch(ip.patch_radius),
    ];

    let simbox = if !config.start_frame().is_empty() {
        // We have an initial position, so just use that
        let particles = particles_from_xyz_nomix(config.start_frame());
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
            cells_per_axis,
            cell_dimensions,
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
    write_geometry_png(&vmmc, "geometry.png"); // TODO: configurable path
}
