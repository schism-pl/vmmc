use equationx::Expr;
use serde::{Deserialize, Serialize};
use std::str::FromStr;
use crate::types::Num;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProtocolStep {
    chemical_potential: Num, // mu / KbT
    interaction_energy: Num, // epsilon / KbT
}

impl ProtocolStep {
    pub fn new(chemical_potential: Num, interaction_energy: Num) -> Self {
        Self {
            chemical_potential,
            interaction_energy,
        }
    }

    pub fn chemical_potential(&self) -> Num {
        self.chemical_potential
    }

    pub fn interaction_energy(&self) -> Num {
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

    pub fn initial_interaction_energy(&self) -> Num {
        self.interaction_energy_eq.eval(0.0)
    }

    pub fn initial_chemical_potential(&self) -> Num {
        self.chemical_potential_eq.eval(0.0)
    }

    pub fn megastep_iter(&self) -> ProtocolMegastepIter {
        ProtocolMegastepIter::new(self)
    }

    pub fn flat_protocol(chemical_potential: Num, interaction_energy: Num, end: usize) -> Self {
        let chemical_potential_s = format!("{}", chemical_potential);
        let interaction_energy_s = format!("{}", interaction_energy);
        Self::new(&chemical_potential_s, &interaction_energy_s, end)
    }
}

// gets a new protocol step every 1000 steps of the simulation
pub struct ProtocolMegastepIter<'a> {
    protocol: &'a SynthesisProtocol,
    t: Num,
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
        if self.t >= self.protocol.num_megasteps as Num {
            return None;
        }

        let chemical_potential = self.protocol.chemical_potential_eq.eval(self.t);
        let interaction_energy = self.protocol.interaction_energy_eq.eval(self.t);
        let step = ProtocolStep::new(chemical_potential, interaction_energy);
        self.t += 1.0; // t is counted in megasteps
        Some(step)
    }
}
