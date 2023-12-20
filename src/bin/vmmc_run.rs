use std::fs::{self, create_dir_all};

use clap::Parser;
use rand::rngs::SmallRng;
use rand::{Rng, SeedableRng};
use vmmc::chemical_potential::maybe_particle_exchange;
use vmmc::cli::VmmcConfig;
use vmmc::consts::MAX_PARTICLES;
use vmmc::io::{clear_out_files, write_geometry_png};
use vmmc::polygons::{calc_polygon_count, Polygon};
use vmmc::protocol::{FixedProtocol, ProtocolStep};
use vmmc::stats::RunStats;
use vmmc::{
    io::{write_tcl, XYZWriter},
    vmmc::Vmmc,
};
use vmmc::{vmmc_from_config, InputParams};

// maps shapes to bond distribution
fn calc_bond_distribution(vmmc: &Vmmc) -> Vec<Vec<usize>> {
    let mut bond_counts_per_shape = Vec::new();
    for shape in vmmc.simbox().shapes() {
        let mut bond_counts = vec![0; shape.patches().len() + 1];
        for p in vmmc.particles().iter() {
            let bond_count = vmmc.determine_interactions(p).iter().count();
            bond_counts[bond_count] += 1;
        }
        bond_counts_per_shape.push(bond_counts);
    }
    bond_counts_per_shape
}

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)

pub trait VmmcCallback {
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats);
}

struct StdCallback {
    writer: Box<XYZWriter>,
}
impl VmmcCallback for StdCallback {
    // runs after every million steps
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats) {
        self.writer.write_xyz_frame(vmmc);

        println!(
            "-----------------------------------------\nStep {:?}",
            (idx + 1) * 1000 * 1000,
        );

        println!("# of particles: {:?}", vmmc.particles().num_particles());
        println!("# of polygons: {:?}", calc_polygon_count(vmmc, 6));
        println!(
            "Interaction Energy (epsilon): {:.4}",
            step.interaction_energy()
        );
        println!("Chemical potential (mu): {:.4}", step.chemical_potential());
        println!(
            "Acceptance ratio: {:.4}",
            run_stats.num_accepts() as f64 / run_stats.num_attempts() as f64
        );
        for (idx, shape_stats) in calc_bond_distribution(vmmc).iter().enumerate() {
            let weighted_sum: usize = shape_stats.iter().enumerate().map(|(i, c)| i * c).sum();
            println!(
                "shape_{:?} bond distribution: {:?} average degree = {:.4}",
                idx,
                shape_stats,
                weighted_sum as f64 / vmmc.particles().num_particles() as f64
            );
        }
    }
}

pub fn run_vmmc_w_callback(
    vmmc: &mut Vmmc,
    protocol: FixedProtocol,
    mut callback: Option<Box<dyn VmmcCallback>>,
    rng: &mut SmallRng,
) {
    for (idx, protocol_step) in protocol.enumerate() {
        let mut run_stats = RunStats::new();
        vmmc.set_interaction_energy(protocol_step.interaction_energy());
        let chemical_potential = protocol_step.chemical_potential();
        for _ in 0..1000 {
            let stats = vmmc.step_n(1000, rng);
            run_stats = stats + run_stats;
            maybe_particle_exchange(vmmc, chemical_potential, rng);
        }
        if let Some(ref mut cb) = callback {
            cb.run(&vmmc, &protocol_step, idx, &run_stats);
        }
    }
}

fn main() -> anyhow::Result<()> {
    env_logger::init();

    // Get commandline arguments
    let config = VmmcConfig::parse();

    let ip = if config.input() != "" {
        let contents = fs::read_to_string(config.input())?;
        toml::from_str(&contents)?
    } else {
        InputParams::default()
    };

    // Seed the rng
    let seed = config.seed();
    println!("Using seed = {:?}", seed);
    let mut rng = SmallRng::seed_from_u64(seed);

    // Generate the simulator
    let mut vmmc = vmmc_from_config(&ip, &mut rng);

    // Init I/O
    println!("Writing output to {}", config.output_dir());
    let out_path = std::path::Path::new(config.output_dir());
    create_dir_all(out_path).unwrap();
    clear_out_files(&config)?;

    // dump full config toml to output directory
    let toml = toml::to_string(&ip).unwrap();
    fs::write(config.toml(), toml).expect("Unable to write file");

    // let mut writer = XYZWriter::new(&config.trajectory());
    // let wr
    let mut writer = Box::new(XYZWriter::new(&config.trajectory()));

    // Check that initial conditions are reasonable
    debug_assert!(vmmc.well_formed());
    println!("Initial configuration ok");

    // Run the simulation
    writer.write_xyz_frame(&vmmc);
    run_vmmc_w_callback(
        &mut vmmc,
        ip.protocol,
        Some(Box::new(StdCallback { writer })),
        &mut rng,
    );

    // Write visualizations to disc
    write_tcl(&vmmc, &config.vmd());
    write_geometry_png(&vmmc, &config.geometry());

    Ok(())
}
