use std::f64::consts::PI;
use std::fs::{self, create_dir_all};
use std::process::exit;

use clap::Parser;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use vmmc::cli::VmmcConfig;
use vmmc::consts::{MAX_INITIAL_PACKING_FRACTION, PARTICLE_RADIUS};
use vmmc::io::{clear_out_files, write_geometry_png};
use vmmc::polygons::{calc_bond_distribution, calc_polygon_count};
use vmmc::protocol::ProtocolStep;
use vmmc::stats::RunStats;
use vmmc::{
    io::{write_tcl, XYZWriter},
    vmmc::Vmmc,
};
use vmmc::{run_vmmc, vmmc_from_config, InputParams, VmmcCallback};

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)

fn packing_fraction(num_particles: usize, volume: f64) -> f64 {
    let particle_volume = num_particles as f64 * (PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    particle_volume / volume
}

struct StdCallback {
    writer: Box<XYZWriter>,
}
impl VmmcCallback for StdCallback {
    type CbResult = ();
    // runs after every million steps
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats) {
        self.writer.write_xyz_frame(vmmc);

        println!(
            "------------------------------------\nStep {:?} x 10e6",
            (idx + 1),
        );

        assert!(vmmc.well_formed());

        println!("# of particles: {:?}", vmmc.particles().num_particles());
        println!(
            "Packing fraction: {:?}",
            packing_fraction(vmmc.particles().num_particles(), vmmc.simbox().volume())
        );
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

    fn state(&self) {}
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
    let seed = ip.seed;
    println!("Using seed = {:x?}", seed);
    let mut rng = SmallRng::seed_from_u64(seed as u64);

    let packing_fraction = packing_fraction(ip.initial_particles, ip.box_height * ip.box_width);
    println!("Requested initial packing fraction: {packing_fraction}");
    if packing_fraction >= MAX_INITIAL_PACKING_FRACTION {
        println!("Initial packing fraction > {MAX_INITIAL_PACKING_FRACTION}");
        println!("Exiting, as it is not likely that we can pack particles this densely");
        exit(1);
    }

    // Generate the simulator
    let mut vmmc = vmmc_from_config(&ip, &mut rng);

    println!("Using {} Mb", vmmc.needed_mem() as f64 / (1024.0 * 1024.0));
    println!(
        "Could use up to {} Mb",
        vmmc.max_needed_mem() as f64 / (1024.0 * 1024.0)
    );

    // Init I/O
    println!("Writing output to {}", config.output_dir());
    let out_path = std::path::Path::new(config.output_dir());
    create_dir_all(out_path).unwrap();
    clear_out_files(&config)?;

    // dump full config toml to output directory
    let toml = toml::to_string(&ip).unwrap();
    fs::write(config.toml(), toml).expect("Unable to write file");

    let mut writer = Box::new(XYZWriter::new(&config.trajectory()));

    // Check that initial conditions are reasonable
    debug_assert!(vmmc.well_formed());
    println!("Initial configuration ok");

    // Write initial frame
    writer.write_xyz_frame(&vmmc);

    // Run the simulation
    let cb = Box::new(StdCallback { writer });
    run_vmmc(&mut vmmc, ip.protocol, cb, &mut rng)?;

    // Write visualizations to disc
    write_tcl(&vmmc, &config.vmd());
    write_geometry_png(&vmmc, &config.geometry());

    println!("Done!");

    Ok(())
}
