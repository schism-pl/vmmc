use clap::Parser;

#[derive(Parser, Debug, Clone)]
pub struct VmmcConfig {
    #[arg(long, default_value = "")]
    input: String,
    #[arg(long, default_value = "./out")]
    output_dir: String,
}

impl VmmcConfig {
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

    pub fn protocols(&self) -> String {
        format!("{}/protocols.png", self.output_dir())
    }

    pub fn trajectory(&self) -> String {
        format!("{}/trajectory.xyz", self.output_dir())
    }

    pub fn stats(&self) -> String {
        format!("{}/stats.txt", self.output_dir())
    }
}
