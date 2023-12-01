use clap::Parser;
use rand::rngs::SmallRng;
use rand::SeedableRng;
use vmmc::cli::VmmcConfig;
use vmmc::evo_vmmc::{EvoVmmc, FitnessFunc};

// correctness criteria:
// 1. average energy monotonically increases (decreases?)
// 2. particles visibly stick together in visualization
// 3. values match other impls (approximately)

// TODO: builder pattern
fn main() {
    env_logger::init();
    // Get commandline arguments
    // TODO: use seperate commandline arguments
    let config = VmmcConfig::parse();
    // Get default params

    // Seed the rng
    let seed = config.seed();
    println!("Using seed = {:?}", seed);
    let mut rng = SmallRng::seed_from_u64(seed);

    // let ip = L2GInputParams::default();
    let mut evo_vmmc = EvoVmmc::new(FitnessFunc::AvgEnergy);
    evo_vmmc.step_all(&mut rng);

    // evo_vmmc.step_generation_n(ip.num_generations, &mut rng);

    // TODO: implement fitness function class
    // TODO: implement pruning functionality
    // TODO: implement Metasim class
    // TODO: parallelize children
    // TODO: randomize starts?
    // TODO: modify something meaningful each run (chemical potential / interaction range?)
    // for idx in 0..ip.num_generations {
    //     println!("Starting generation {:?}", idx);
    //     for jdx in 0..ip.children_per_generation {
    //         println!("Simulating child {:?}-{:?}", idx, jdx);
    //         let mut vmmc = vmmc_from_config(&ip, &mut rng);
    //         let avg_energy = run_sim(&mut vmmc, &ip, &mut rng);
    //         println!("average energy of sim {:?}-{:?} = {:?} kBT\n", idx, jdx, avg_energy);
    //     }
    // }
}
