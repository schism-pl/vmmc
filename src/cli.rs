use clap::Parser;
use rand::{rngs::SmallRng, Rng, SeedableRng};

#[derive(Parser, Debug)]
pub struct VmmcConfig {
    #[arg(long, default_value_t = SmallRng::from_entropy().gen::<u64>())]
    seed: u64,
    #[arg(long, default_value = "")]
    input: String,
    #[arg(long, default_value = "./out")]
    output_dir: String,
}

impl VmmcConfig {
    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn input(&self) -> &str {
        &self.input
    }

    pub fn output_dir(&self) -> &str {
        &self.output_dir
    }

    pub fn set_output_dir(&mut self, s: &str) {
        self.output_dir = s.to_string();
    }

    pub fn toml(&self) -> String {
        format!("{}/config.toml", self.output_dir())
    }

    pub fn vmd(&self) -> String {
        format!("{}/vmd.tcl", self.output_dir())
    }

    pub fn geometry(&self) -> String {
        format!("{}/geometry.png", self.output_dir())
    }

    pub fn trajectory(&self) -> String {
        format!("{}/trajectory.xyz", self.output_dir())
    }
}
