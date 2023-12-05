#[derive(Debug, Clone)]
pub struct ProtocolStep {
    chemical_potential: f64, // mu / KbT
    interaction_energy: f64, // epsilon / KbT
}

impl ProtocolStep {
    pub fn new(chemical_potential: f64, interaction_energy: f64) -> Self {
        Self {
            chemical_potential,
            interaction_energy,
        }
    }

    pub fn chemical_potential(&self) -> f64 {
        self.chemical_potential
    }

    pub fn interaction_energy(&self) -> f64 {
        self.interaction_energy
    }
}

trait Protocol: Iterator<Item = ProtocolStep> {}
impl<T: Iterator<Item = ProtocolStep>> Protocol for T {}

#[derive(Debug, Clone)]
pub struct FixedProtocol {
    inner: Vec<ProtocolStep>,
    t: usize,
}

impl Iterator for FixedProtocol {
    type Item = ProtocolStep;

    fn next(&mut self) -> Option<Self::Item> {
        if self.t >= self.inner.len() {
            return None;
        }

        let step = self.inner[self.t].clone();
        self.t += 1;
        Some(step)
    }
}

impl FixedProtocol {
    pub fn new(inner: Vec<ProtocolStep>) -> Self {
        Self { inner, t: 0 }
    }

    fn initial_config(&self) -> ProtocolStep {
        self.inner[0].clone()
    }

    pub fn initial_interaction_energy(&self) -> f64 {
        self.initial_config().interaction_energy()
    }

    pub fn initial_chemical_potential(&self) -> f64 {
        self.initial_config().chemical_potential()
    }

    pub fn flat_protocol(chemical_potential: f64, interaction_energy: f64, len: usize) -> Self {
        let inner = vec![ProtocolStep::new(chemical_potential, interaction_energy); len];
        Self::new(inner)
    }
}
