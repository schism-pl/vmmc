use std::f64::consts::PI;
use std::fs::{self, create_dir_all};
use std::process::exit;
use std::time::Instant;

use vmmc::protocol::ProtocolStep;
use vmmc::stats::RunStats;
use std::env;
use vmmc::chemical_potential::maybe_particle_exchange;
use serde::{Serialize, Deserialize};
use vmmc::vmmc::Vmmc;
use rand::rngs::SmallRng;

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)

// fn packing_fraction(num_particles: usize, volume: f64) -> f64 {
//     let particle_volume = num_particles as f64 * (PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
//     particle_volume / volume
// }

// struct StdCallback {
//     writer: Box<XYZWriter>,
//     start_time: Instant,
//     timestamp: Instant,
// }
// impl VmmcCallback for StdCallback {
//     type CbResult = ();
//     // runs after every million steps
//     fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats) {
//         self.writer.write_xyz_frame(vmmc);

//         println!(
//             "------------------------------------\nStep {:?} x 10e6",
//             (idx + 1),
//         );

//         assert!(vmmc.well_formed());

//         println!("# of particles: {:?}", vmmc.particles().num_particles());
//         println!(
//             "Packing fraction: {:?}",
//             packing_fraction(vmmc.particles().num_particles(), vmmc.simbox().volume())
//         );
//         let polygon_dist = calc_polygon_distribution(vmmc, 12);
//         println!("Polygon distribution: {:?}", polygon_dist);
//         println!("Total polygons: {:?}", polygon_dist.iter().sum::<usize>());
//         println!(
//             "Interaction Energy (epsilon): {:.4}",
//             step.interaction_energy()
//         );
//         println!("Chemical potential (mu): {:.4}", step.chemical_potential());
//         println!(
//             "Acceptance ratio: {:.4}",
//             run_stats.num_accepts() as f64 / run_stats.num_attempts() as f64
//         );
//         for (idx, shape_stats) in calc_bond_distribution(vmmc).iter().enumerate() {
//             let weighted_sum: usize = shape_stats.iter().enumerate().map(|(i, c)| i * c).sum();
//             println!(
//                 "shape_{:?} bond distribution: {:?} average degree = {:.4}",
//                 idx,
//                 shape_stats,
//                 weighted_sum as f64 / vmmc.particles().num_particles() as f64
//             );
//         }
//         let t = Instant::now();
//         println!("Execution time: {:?}", t - self.timestamp);
//         println!("Total execution time: {:?}", t - self.start_time);
//         self.timestamp = t;
//     }

//     fn state(&self) {}
// }

#[derive(Serialize, Deserialize)]
struct Snapshot {
    vmmc: Vmmc,
    step: ProtocolStep, 
    rng: SmallRng,
}

fn main() -> anyhow::Result<()> {
    env_logger::init();

    // // Get commandline arguments
    // let config = VmmcConfig::parse();

    // let ip = if config.input() != "" {
    //     println!("Using requested config: = {}", config.input());
    //     let contents = fs::read_to_string(config.input())?;
    //     toml::from_str(&contents)?
    // } else {
    //     println!("No config provided. Using default config");
    //     InputParams::default()
    // };

    // // Seed the rng
    // let seed = ip.seed;
    // println!("Using seed = {:x?}", seed);
    // let mut rng = SmallRng::seed_from_u64(seed as u64);

    // let packing_fraction = packing_fraction(
    //     ip.sim_params.initial_particles,
    //     ip.sim_params.box_height * ip.sim_params.box_width,
    // );
    // println!("Requested initial packing fraction: {packing_fraction}");
    // if packing_fraction >= MAX_INITIAL_PACKING_FRACTION {
    //     println!("Initial packing fraction > {MAX_INITIAL_PACKING_FRACTION}");
    //     println!("Exiting, as it is not likely that we can pack particles this densely");
    //     exit(1);
    // }

    // // Generate the simulator
    // let mut vmmc = vmmc_from_inputparams(&ip, &mut rng);

    // println!("Using {} Mb", vmmc.needed_mem() as f64 / (1024.0 * 1024.0));
    // println!(
    //     "Could use up to {} Mb",
    //     vmmc.max_needed_mem() as f64 / (1024.0 * 1024.0)
    // );

    // // Init I/O
    // println!("Writing output to {}", config.output_dir());
    // let out_path = std::path::Path::new(config.output_dir());
    // create_dir_all(out_path).unwrap();
    // clear_out_files(&config)?;

    // // dump full config toml to output directory
    // let toml = toml::to_string(&ip).unwrap();
    // fs::write(config.toml(), toml).expect("Unable to write file");

    // let mut writer = Box::new(XYZWriter::new(&config.trajectory()));

    // // Check that initial conditions are reasonable
    // debug_assert!(vmmc.well_formed());
    // println!("Initial configuration ok");

    // // Write initial frame
    // writer.write_xyz_frame(&vmmc);
    let args: Vec<String> = env::args().collect();
    let filename = args[1];

    println!("Using snapshot {}", filename);

    let contents = fs::read_to_string(filename)?;
    let state: Snapshot = toml::from_str(&contents)?;

    let vmmc = state.vmmc;
    let protocol_step = state.step;

    // Run the simulation
    let mut run_stats = RunStats::new();
    vmmc.set_interaction_energy(protocol_step.interaction_energy());
    let chemical_potential = protocol_step.chemical_potential();
    for _ in 0..(1000 * 1000) {
        let _ = vmmc.step(rng, &mut run_stats);
        if vmmc.dynamic_particle_count() {
            maybe_particle_exchange(&mut vmmc, chemical_potential, rng);
        }
    }

    
    // Write visualizations to disc
    // write_tcl(&vmmc, &config.vmd());
    // write_geometry_png(&vmmc, &config.geometry());
    // write_colored_geometry_png(&vmmc, &config.colored_geometry());
    // write_protocols_png(protocol, &config.protocols());

    // write_stats(&vmmc, &config.stats());

    println!("Done!");

    Ok(())
}
