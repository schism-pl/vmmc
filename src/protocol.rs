use serde::{Deserialize, Serialize};
use equations::Equation;
use std::str::FromStr;

#[derive(Debug, Clone, Serialize, Deserialize)]
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

// trait Protocol: Iterator<Item = ProtocolStep> {}
// impl<T: Iterator<Item = ProtocolStep>> Protocol for T {}

// TODO: [(Equation, Duration)]
#[derive(Clone, Serialize, Deserialize)]
pub struct SynthesisProtocol {
    chemical_potential_eq: Equation,
    interaction_energy_eq: Equation,
    t: f64,
    end: usize, // in megasteps
}

// kilostep iterator
impl Iterator for SynthesisProtocol {
    type Item = ProtocolStep;

    fn next(&mut self) -> Option<Self::Item> {
        if self.t >= self.end as f64 {
            return None;
        }

        let chemical_potential = self.chemical_potential_eq.eval(self.t as f64);
        let interaction_energy = self.interaction_energy_eq.eval(self.t as f64);
        let step = ProtocolStep::new(chemical_potential, interaction_energy);
        self.t += 0.001; // t is counted in megasteps, but we iter every 1000 steps
        Some(step)
    }
}

impl SynthesisProtocol {
    pub fn new(chemical_potential_s: &str, interaction_energy_s: &str, end: usize) -> Self {
        let chemical_potential_eq = Equation::from_str(chemical_potential_s).unwrap();
        let interaction_energy_eq = Equation::from_str(interaction_energy_s).unwrap();
        Self { chemical_potential_eq, interaction_energy_eq, t: 0.0, end }
    }

    pub fn initial_interaction_energy(&self) -> f64 {
        self.interaction_energy_eq.eval(0.0)
    }

    pub fn initial_chemical_potential(&self) -> f64 {
        self.chemical_potential_eq.eval(0.0)
    }

    pub fn flat_protocol(chemical_potential: f64, interaction_energy: f64, end: usize) -> Self {
        let chemical_potential_s = format!("mu = {}", chemical_potential);
        let interaction_energy_s = format!("tau = {}", interaction_energy);
        Self::new(&chemical_potential_s, &interaction_energy_s, end)
    }
}
