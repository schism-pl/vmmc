use equationx::Expr;
use serde::{Deserialize, Serialize};
use std::fmt;
use std::str::FromStr;

use crate::vmmc::Vmmc;

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

#[derive(Clone, Serialize, Deserialize)]
pub struct SynthesisProtocol {
    chemical_potential_eq: Expr<f64>,
    interaction_energy_eq: Expr<f64>,
    num_megasteps: usize,
}

impl SynthesisProtocol {
    pub fn new(
        chemical_potential_s: &str,
        interaction_energy_s: &str,
        num_megasteps: usize,
    ) -> Self {
        let chemical_potential_eq = Expr::from_str(chemical_potential_s).unwrap();
        let interaction_energy_eq = Expr::from_str(interaction_energy_s).unwrap();
        Self {
            chemical_potential_eq,
            interaction_energy_eq,
            num_megasteps,
        }
    }

    pub fn len(&self) -> usize {
        self.num_megasteps
    }

    pub fn interaction_energy(&self, t: usize) -> f64 {
        self.interaction_energy_eq.eval(t as f64)
    }

    pub fn chemical_potential(&self, t: usize) -> f64 {
        self.chemical_potential_eq.eval(t as f64)
    }

    pub fn initial_interaction_energy(&self) -> f64 {
        self.interaction_energy_eq.eval(0.0)
    }

    pub fn initial_chemical_potential(&self) -> f64 {
        self.chemical_potential_eq.eval(0.0)
    }

    pub fn num_megasteps(&self) -> usize {
        self.num_megasteps
    }

    pub fn megastep_iter(&self) -> ProtocolMegastepIter {
        ProtocolMegastepIter::new(self)
    }

    pub fn flat_protocol(chemical_potential: f64, interaction_energy: f64, end: usize) -> Self {
        let chemical_potential_s = format!("{}", chemical_potential);
        let interaction_energy_s = format!("{}", interaction_energy);
        Self::new(&chemical_potential_s, &interaction_energy_s, end)
    }
}

impl fmt::Display for SynthesisProtocol {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "chemical_potential: {} interaction_energy: {} t={}",
            self.chemical_potential_eq, self.interaction_energy_eq, self.num_megasteps
        )
    }
}

pub trait ProtocolIter {
    fn start(&self) -> ProtocolStep;
    fn peek(&self, vmmc: &Vmmc) -> ProtocolStep;
    fn next(&mut self, vmmc: &Vmmc) -> Option<ProtocolStep>;
    fn len(&self) -> usize;
}

// gets a new protocol step every 1000 steps of the simulation
pub struct ProtocolMegastepIter {
    protocol: SynthesisProtocol,
    t: f64,
}

impl ProtocolMegastepIter {
    fn new(protocol: &SynthesisProtocol) -> Self {
        Self {
            protocol: protocol.clone(),
            t: 0.0,
        }
    }
}

// next, len, peek

// kilostep iterator
impl ProtocolIter for ProtocolMegastepIter {
    // type Item = ProtocolStep;

    fn next(&mut self, _vmmc: &Vmmc) -> Option<ProtocolStep> {
        if self.t >= self.protocol.num_megasteps as f64 {
            return None;
        }

        let chemical_potential = self.protocol.chemical_potential_eq.eval(self.t);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(self.t);
        let step = ProtocolStep::new(chemical_potential, interaction_energy);
        self.t += 1.0; // t is counted in megasteps
        Some(step)
    }

    fn start(&self) -> ProtocolStep {
        let chemical_potential = self.protocol.chemical_potential_eq.eval(0.0);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(0.0);
        ProtocolStep::new(chemical_potential, interaction_energy)
    }

    fn peek(&self, _vmmc: &Vmmc) -> ProtocolStep {
        let chemical_potential = self.protocol.chemical_potential_eq.eval(self.t);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(self.t);
        ProtocolStep::new(chemical_potential, interaction_energy)
    }

    fn len(&self) -> usize {
        self.protocol.num_megasteps() - self.t as usize
    }
}
