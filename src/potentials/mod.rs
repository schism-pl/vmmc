use crate::{particle::IsParticle, simbox::SimBox};

pub mod gradated_patches;
pub mod patchy_discs;

pub trait Potential {
    fn interaction_energy(&mut self) -> &mut f64;

    fn compute_pair_energy<P1: IsParticle, P2: IsParticle>(
        &self,
        simbox: &SimBox,
        particle0: &P1,
        particle1: &P2,
    ) -> f64;
}
