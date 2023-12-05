use clap::Parser;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use vmmc::cli::VmmcConfig;
use vmmc::io::write_geometry_png;
use vmmc::morphology::Morphology;
use vmmc::position::DimVec;
use vmmc::protocol::{FixedProtocol, ProtocolStep};
use vmmc::{
    io::{read_xyz_snapshot, write_tcl, XYZWriter},
    particle::Particle,
    simbox::SimBox,
    vmmc::{Vmmc, VmmcParams},
};
use vmmc::{protocol, InputParams};

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

    let max_interaction_range = 1.0 + ip.max_patch_radius();

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
            ip.shapes.clone(),
        )
    } else {
        // No initial position, so we will use a randomized start position
        SimBox::new_with_randomized_particles(
            box_dimensions,
            max_interaction_range,
            ip.num_particles,
            ip.shapes.clone(),
            rng,
        )
    };

    let params = VmmcParams::new(ip.prob_translate, ip.max_translation, ip.max_rotation);
    // let interaction_energy = ip.protocol.next().unwrap().interaction_energy();
    let interaction_energy = ip.protocol.initial_interaction_energy();
    Vmmc::new(simbox, params, interaction_energy)
}

fn maybe_remove_particle(vmmc: &mut Vmmc, rng: &mut SmallRng) {
    panic!("TODO")

    // if (vmmc.particles().len() > 0) {
    //     remove_particle();
    // }
}

fn maybe_insert_particle(vmmc: &mut Vmmc, rng: &mut SmallRng) {
    panic!("TODO")
    // if (vmmc.particles().len() < MAX_PARTICLES) {
    //     insert_particle();
    // }
}

fn particle_exchange(vmmc: &mut Vmmc, rng: &mut SmallRng) {
    if (rng.gen::<f64>() < 0.5) {
        maybe_remove_particle(vmmc, rng);
    } else {
        maybe_insert_particle(vmmc, rng);
    }
}

// TODO: better name?
// equations G1-G3 in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.011044
// TODO: this equation satisfies detailed balance, but I'm not quite sure how/why
fn maybe_particle_exchange(vmmc: &mut Vmmc, rng: &mut SmallRng) {
    let p_exchange = 1.0 / (1.0 + vmmc.particles().len() as f64);
    if rng.gen::<f64>() < p_exchange {
        particle_exchange(vmmc, rng);
    }
}

fn run_vmmc(
    vmmc: &mut Vmmc,
    mut protocol: FixedProtocol,
    writer: &mut XYZWriter,
    num_sweeps: usize,
    rng: &mut SmallRng,
) {
    // TODO: num sweeps is not the right name for this
    for idx in 0..num_sweeps {
        writer.write_xyz_frame(&vmmc);
        let protocol_update = protocol.next().unwrap();
        vmmc.set_interaction_energy(protocol_update.interaction_energy());
        let stats = vmmc.step_n(1000 * 1000, rng);
        println!(
            "Step {:?}: bonds per particle = {:?}",
            (idx + 1) * 1000 * 1000,
            -vmmc.get_average_energy() / vmmc.potential().interaction_energy()
        );
        println!(
            "Acceptance ratio: {:?}",
            stats.num_accepts() as f64 / (1000.0 * 1000.0)
        );
    }
    // write the final frame
    writer.write_xyz_frame(&vmmc);
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

    // let protocol = FixedProtocol::flat_protocol(0.0, ip.interaction_energy, ip.num_sweeps);
    run_vmmc(&mut vmmc, ip.protocol, &mut writer, ip.num_sweeps, &mut rng);
    // Write script to generate animation
    write_tcl(&vmmc, config.vmd_output());
    write_geometry_png(&vmmc, "geometry.png"); // TODO: configurable path
}
