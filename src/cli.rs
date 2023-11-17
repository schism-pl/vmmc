use clap::Parser;

#[derive(Parser, Debug)]
pub struct VmmcConfig {
    #[arg(long, default_value = None)]
    seed: Option<u64>,
    #[arg(long, default_value = "snapshot.xyz")]
    start_frame: String,
    #[arg(long, default_value = "vmd.tcl")]
    vmd_output: String,
    #[arg(long, default_value = "trajectory.xyz")]
    xyz_output: String,
}

impl VmmcConfig {
    pub fn seed(&self) -> Option<u64> {
        self.seed
    }

    pub fn start_frame(&self) -> &str {
        &self.start_frame
    }

    pub fn vmd_output(&self) -> &str {
        &self.vmd_output
    }

    pub fn xyz_output(&self) -> &str {
        &self.xyz_output
    }
}
