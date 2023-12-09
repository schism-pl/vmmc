use clap::Parser;
use rand::{rngs::SmallRng, Rng, SeedableRng};

#[derive(Parser, Debug)]
pub struct VmmcConfig {
    #[arg(long, default_value_t = SmallRng::from_entropy().gen::<u64>())]
    seed: u64,
    #[arg(long, default_value = "")]
    start_frame: String,
    // #[arg(long, default_value = "out/")]
    // output_dir: String,
    #[arg(long, default_value = "vmd.tcl")]
    vmd_output: String,
    #[arg(long, default_value = "trajectory.xyz")]
    xyz_output: String,
}

impl VmmcConfig {
    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn start_frame(&self) -> &str {
        &self.start_frame
    }

    // pub fn output_dir(&self) -> &str {
    //     &self.output_dir
    // }

    pub fn vmd_output(&self) -> &str {
        &self.vmd_output
    }

    pub fn xyz_output(&self) -> &str {
        &self.xyz_output
    }
}
