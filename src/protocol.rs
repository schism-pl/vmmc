use equationx::Expr;
use serde::{Deserialize, Serialize};
use std::fmt;
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
    chemical_potential_eq: Expr,
    interaction_energy_eq: Expr,
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

// gets a new protocol step every 1000 steps of the simulation
pub struct ProtocolMegastepIter<'a> {
    protocol: &'a SynthesisProtocol,
    t: f64,
}

impl<'a> ProtocolMegastepIter<'a> {
    fn new(protocol: &'a SynthesisProtocol) -> Self {
        Self { protocol, t: 0.0 }
    }
}

// kilostep iterator
impl<'a> Iterator for ProtocolMegastepIter<'a> {
    type Item = ProtocolStep;

    fn next(&mut self) -> Option<Self::Item> {
        if self.t >= self.protocol.num_megasteps as f64 {
            return None;
        }

        let chemical_potential = self.protocol.chemical_potential_eq.eval(self.t);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(self.t);
        let step = ProtocolStep::new(chemical_potential, interaction_energy);
        self.t += 1.0; // t is counted in megasteps
        Some(step)
    }
}

impl<'a> ExactSizeIterator for ProtocolMegastepIter<'a> {
    // We can easily calculate the remaining number of iterations.
    fn len(&self) -> usize {
        self.protocol.num_megasteps() - self.t as usize
    }
}

pub trait Peekable {
    type Output;
    fn peek(&self) -> Self::Output;
}

impl<'a> Peekable for ProtocolMegastepIter<'a> {
    type Output = ProtocolStep;
    fn peek(&self) -> Self::Output {
        let chemical_potential = self.protocol.chemical_potential_eq.eval(self.t);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(self.t);
        ProtocolStep::new(chemical_potential, interaction_energy)
    }
}

pub trait ProtocolIter:
    Iterator<Item = ProtocolStep> + ExactSizeIterator + Peekable<Output = ProtocolStep>
{
}

impl<T: Iterator<Item = ProtocolStep> + ExactSizeIterator + Peekable<Output = ProtocolStep>>
    ProtocolIter for T
{
}
