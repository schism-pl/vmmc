#![allow(clippy::needless_range_loop)]

use anyhow::Result;
use chemical_potential::maybe_particle_exchange;
use consts::PARTICLE_DIAMETER;
use morphology::{Morphology, Patch};
use position::DimVec;
use protocol::{ProtocolIter, ProtocolStep, SynthesisProtocol};
use quickcheck::{Arbitrary, Gen};
use rand::{rngs::SmallRng, Rng, SeedableRng};
use rand_distr::num_traits::Zero;
use serde::{Deserialize, Serialize};
use simbox::SimBox;
use stats::RunStats;
use vmmc::Vmmc;

pub mod chemical_potential;
pub mod cli;
pub mod consts;
pub mod io;
pub mod morphology;
pub mod particle;
pub mod patchy_discs;
pub mod polygons;
pub mod position;
pub mod protocol;
pub mod simbox;
pub mod stats;
pub mod vmmc;

// TODO: use approx crate for floating point equality

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct InputParams {
    // TODO: change this to u64
    pub seed: i64, // toml crashes when I try to store as u64?

    pub initial_particles: usize,
    pub protocol: SynthesisProtocol,
    pub shapes: Vec<Morphology>,

    pub box_width: f64,
    pub box_height: f64,
}

impl Default for InputParams {
    fn default() -> Self {
        let seed = SmallRng::from_entropy().gen::<i64>();

        let initial_particles = 400;
        let box_width = 30.0;
        let box_height = 30.0;

        let protocol = SynthesisProtocol::flat_protocol(0.0, 8.0, 1000);

        let shapes = vec![Morphology::regular_4patch(0.05)];

        Self {
            seed,
            initial_particles,
            protocol,
            shapes,
            box_width,
            box_height,
        }
    }
}

impl InputParams {
    pub fn max_patch_radius(&self) -> f64 {
        self.shapes
            .iter()
            .map(|m| m.max_patch_radius())
            .fold(f64::MIN, |a, b| a.max(b))
    }

    // assert well-formedness predicate for InputParams
    // Note: can stash all assumptions about layout here
    // Note: this needs to be synced up with docs
    pub fn check(&self) {
        // Basic range checks
        assert!(self.initial_particles <= 2500);
        assert!(self.box_width >= 10.0 && self.box_height >= 10.0);
        assert!(self.box_width <= 200.0 && self.box_height <= 200.0);

        // check well-formedness of protocol
        for step in self.protocol.megastep_iter() {
            let mu = step.chemical_potential();
            let epsilon = step.interaction_energy();
            assert!((-10.0..=10.0).contains(&mu));
            assert!((0.01..=20.0).contains(&epsilon));
        }

        // make sure our patches are reasonable
        for shape in &self.shapes {
            assert!(shape.patches().len() >= 3);
            assert!(shape.patches().len() <= 6);
            for patch in shape.patches() {
                assert!(patch.theta() >= 0.0 && patch.theta() > 360.0);
                assert!(patch.radius() >= 0.01 && patch.radius() <= 0.25);
            }
        }
    }
}

// for arbitrary trait
fn usize_in_range(g: &mut Gen, min: usize, max: usize) -> usize {
    if min == max {
        return max;
    }
    let x = usize::arbitrary(g);
    let r = x % (max - min) + min;
    assert!(r >= min && r < max);
    r
}

fn f64_in_range(g: &mut Gen, min: f64, max: f64) -> f64 {
    let mut r = f64::INFINITY;
    while !(r.is_normal() || r.is_zero()) {
        let x = f64::arbitrary(g).abs();
        r = x % (max - min) + min;
    }

    assert!(r.is_normal() || r.is_zero());
    assert!(r >= min && r < max);
    r
}

// num_cells = box_x * box_y / (1 + patch_radius)^2
// we need num_cells >= 2*initial_particles

// for testing
impl Arbitrary for InputParams {
    fn arbitrary(g: &mut Gen) -> Self {
        let seed = SmallRng::from_entropy().gen::<i64>();

        let initial_particles = usize_in_range(g, 0, 2500);
        let interaction_energy = f64_in_range(g, 0.01, 20.0); // kBT
        let patch_radius = f64_in_range(g, 0.01, 0.1); // radius of patch (in units of particle diameter)

        let mut box_width = 0.0;
        let mut box_height = 0.0;

        // while less than 1.5 cell per particle, pick a new box radius
        while box_width * box_height / (PARTICLE_DIAMETER + patch_radius + patch_radius).powi(2)
            <= 1.5 * initial_particles as f64
        {
            box_width = f64_in_range(g, 10.0, 200.0);
            box_height = f64_in_range(g, 10.0, 200.0);
        }

        let num_megasteps = 10;

        // TODO: make this more varied
        let protocol = SynthesisProtocol::flat_protocol(0.0, interaction_energy, num_megasteps);

        let mut shapes = Vec::new();
        let num_shapes = usize_in_range(g, 1, 3);
        for _ in 0..num_shapes {
            let mut patches = Vec::new();
            let num_patches = usize_in_range(g, 3, 6);
            let num_colors = usize_in_range(g, 1, num_patches);
            for _ in 0..num_patches {
                let color = usize_in_range(g, 0, num_colors - 1);
                let theta = f64_in_range(g, 0.0, 360.0);
                let radius = f64_in_range(g, 0.01, 0.25);
                let patch = Patch::new(radius, theta, color as u8);
                // TODO: check that the morphology makes sense, i.e., patches don't overlap
                patches.push(patch);
            }
            shapes.push(Morphology::new(morphology::CoreShape::Circle, patches));
        }

        Self {
            seed,
            initial_particles,
            protocol,
            shapes,
            box_width,
            box_height,
        }
    }
}

// angle coverage of patch = +- acos 1-(patch_radius^2/2)
impl Arbitrary for SimBox {
    fn arbitrary(g: &mut Gen) -> Self {
        let ip = InputParams::arbitrary(g);
        let box_dimensions = DimVec::new([ip.box_width, ip.box_height]);
        let patch_radius = 0.05; // TODO: remove this
        let max_interaction_range = PARTICLE_DIAMETER + patch_radius + patch_radius;

        let mut rng = SmallRng::from_entropy();
        SimBox::new_with_randomized_particles(
            box_dimensions,
            max_interaction_range,
            ip.initial_particles,
            ip.shapes,
            &mut rng,
        )
    }
}

pub fn vmmc_from_config(ip: &InputParams, rng: &mut SmallRng) -> vmmc::Vmmc {
    let box_dimensions = DimVec::new([ip.box_width, ip.box_height]);

    let max_interaction_range = PARTICLE_DIAMETER + ip.max_patch_radius();

    // No initial position, so we will use a randomized start position
    let simbox = SimBox::new_with_randomized_particles(
        box_dimensions,
        max_interaction_range,
        ip.initial_particles,
        ip.shapes.clone(),
        rng,
    );

    let interaction_energy = ip.protocol.initial_interaction_energy();
    vmmc::Vmmc::new(simbox, interaction_energy)
}

pub trait VmmcCallback {
    type CbResult;
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats);
    fn state(&self) -> Self::CbResult;
}

pub struct NoCallback {}
impl VmmcCallback for NoCallback {
    type CbResult = ();
    fn run(&mut self, _: &Vmmc, _: &ProtocolStep, _: usize, _: &RunStats) {}
    fn state(&self) {}
}

/// helper for when we use no callback
pub fn no_callback() -> Box<dyn VmmcCallback<CbResult = ()>> {
    Box::new(NoCallback {})
}

// updates chemical potential and interaction energy every 10e6 steps
// attempts particle exchange every step
pub fn run_vmmc<Cbr>(
    vmmc: &mut Vmmc,
    protocol_iter: impl ProtocolIter,
    mut cb: Box<dyn VmmcCallback<CbResult = Cbr>>,
    rng: &mut SmallRng,
) -> Result<Cbr> {
    for (idx, protocol_step) in protocol_iter.enumerate() {
        let mut run_stats = RunStats::new();
        vmmc.set_interaction_energy(protocol_step.interaction_energy());
        let chemical_potential = protocol_step.chemical_potential();
        for _ in 0..(1000 * 1000) {
            let _ = vmmc.step(rng, &mut run_stats);
            maybe_particle_exchange(vmmc, chemical_potential, rng);
        }
        cb.run(vmmc, &protocol_step, idx, &run_stats);
    }
    Ok(cb.state())
}
