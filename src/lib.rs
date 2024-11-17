#![allow(clippy::needless_range_loop)]

use std::{f64::consts::PI, time::Instant};

use anyhow::Result;
use chemical_potential::maybe_particle_exchange;
use consts::{PARTICLE_DIAMETER, PARTICLE_RADIUS};
use morphology::{Morphology, Patch};
use polygons::{calc_bond_distribution, calc_polygon_distribution};
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
pub mod polygons;
pub mod position;
pub mod potentials;
pub mod protocol;
pub mod simbox;
pub mod stats;
pub mod vmmc;

// TODO: use approx crate for floating point equality

#[derive(Clone, Serialize, Deserialize)]
pub struct SimParams {
    pub initial_particles: usize,
    pub shapes: Vec<Morphology>,
    pub box_width: f64,
    pub box_height: f64,
    pub dynamic_particle_count: bool,
}

impl Default for SimParams {
    fn default() -> Self {
        let initial_particles = 550;
        let shapes = vec![Morphology::regular_4patch(0.05)];

        SimParams {
            initial_particles,
            shapes,
            box_width: 30.0,
            box_height: 30.0,
            dynamic_particle_count: true,
        }
    }
}

impl Arbitrary for SimParams {
    fn arbitrary(g: &mut Gen) -> Self {
        let initial_particles = usize_in_range(g, 0, 2500);

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
            shapes.push(Morphology::new(patches));
        }

        Self {
            initial_particles,
            shapes,
            box_width,
            box_height,
            dynamic_particle_count: false,
        }
    }
}

impl SimParams {
    pub fn max_interaction_range(&self) -> f64 {
        PARTICLE_DIAMETER
            + self
                .shapes
                .iter()
                .map(|m| m.max_patch_radius())
                .fold(f64::MIN, |a, b| a.max(b))
    }

    pub fn box_dimensions(&self) -> DimVec {
        DimVec::new([self.box_width, self.box_height])
    }

    pub fn check(&self) {
        // Basic range checks
        assert!(self.initial_particles <= 2500);
        assert!(self.box_width >= 10.0 && self.box_height >= 10.0);
        assert!(self.box_width <= 200.0 && self.box_height <= 200.0);

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

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct InputParams {
    pub seed: u32,
    pub protocol: SynthesisProtocol,
    pub sim_params: SimParams,
}

impl Default for InputParams {
    fn default() -> Self {
        let seed = SmallRng::from_entropy().gen::<u32>();

        // let initial_particles = 400;
        // let box_width = 30.0;
        // let box_height = 30.0;

        let protocol = SynthesisProtocol::flat_protocol(0.0, 8.0, 100);

        // let shapes = vec![Morphology::regular_4patch(0.05)];

        Self {
            seed,
            sim_params: Default::default(),
            // initial_particles,
            protocol,
            // shapes,
            // box_width,
            // box_height,
        }
    }
}

impl InputParams {
    // assert well-formedness predicate for InputParams
    // Note: can stash all assumptions about layout here
    // Note: this needs to be synced up with docs
    pub fn check(&self) {
        self.sim_params.check();

        // check well-formedness of protocol
        for step in self.protocol.megastep_iter() {
            let mu = step.chemical_potential();
            let epsilon = step.interaction_energy();
            assert!((-10.0..=10.0).contains(&mu));
            assert!((0.01..=20.0).contains(&epsilon));
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
        // TODO: seed should be generated from arbitrary?
        let seed = SmallRng::from_entropy().gen::<u32>();

        let num_megasteps = 10;

        let interaction_energy = f64_in_range(g, 0.01, 20.0); // kBT
        let chemical_potential = f64_in_range(g, -10.0, 10.0); // kBT
        let protocol =
            SynthesisProtocol::flat_protocol(chemical_potential, interaction_energy, num_megasteps);

        Self {
            seed,
            sim_params: SimParams::arbitrary(g),
            protocol,
        }
    }
}

// angle coverage of patch = +- acos 1-(patch_radius^2/2)
impl Arbitrary for SimBox {
    fn arbitrary(g: &mut Gen) -> Self {
        let sim_params = SimParams::arbitrary(g);
        let mut rng = SmallRng::from_entropy();
        SimBox::new_with_randomized_particles(&sim_params, &mut rng)
    }
}

pub fn vmmc_from_simparams(
    sim_params: &SimParams,
    initial_interaction_energy: f64,
    rng: &mut SmallRng,
) -> vmmc::Vmmc {
    let simbox = SimBox::new_with_randomized_particles(sim_params, rng);

    vmmc::Vmmc::new(
        simbox,
        initial_interaction_energy,
        sim_params.dynamic_particle_count,
    )
}

pub fn vmmc_from_inputparams(ip: &InputParams, rng: &mut SmallRng) -> vmmc::Vmmc {
    let simbox = SimBox::new_with_randomized_particles(&ip.sim_params, rng);

    vmmc::Vmmc::new(
        simbox,
        ip.protocol.initial_interaction_energy(),
        ip.sim_params.dynamic_particle_count,
    )
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
        println!("Running VMMC megastep {idx}!");
        let mut run_stats = RunStats::new();
        vmmc.set_interaction_energy(protocol_step.interaction_energy());
        let chemical_potential = protocol_step.chemical_potential();
        for jdx in 0..(1000 * 1000) {
            let _ = vmmc.step(rng, &mut run_stats);
            if vmmc.dynamic_particle_count() {
                maybe_particle_exchange(vmmc, chemical_potential, rng);
            }
        }
        cb.run(vmmc, &protocol_step, idx, &run_stats);
    }
    Ok(cb.state())
}


pub fn packing_fraction(num_particles: usize, volume: f64) -> f64 {
    let particle_volume = num_particles as f64 * (PI * PARTICLE_RADIUS * PARTICLE_RADIUS);
    particle_volume / volume
}

pub struct StdCallback {
    start_time: Instant,
    timestamp: Instant,
}

impl StdCallback {
    pub fn new() -> Self {
        StdCallback {
            start_time: Instant::now(),
            timestamp: Instant::now(),
        }
    }
}

impl VmmcCallback for StdCallback {
    type CbResult = ();
    // runs after every million steps
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats) {
        println!(
            "------------------------------------\nStep {:?} x 10e6",
            (idx + 1),
        );

        assert!(vmmc.well_formed());

        println!("# of particles: {:?}", vmmc.particles().num_particles());
        println!(
            "Packing fraction: {:?}",
            packing_fraction(vmmc.particles().num_particles(), vmmc.simbox().volume())
        );
        let polygon_dist = calc_polygon_distribution(vmmc, 12);
        println!("Polygon distribution: {:?}", polygon_dist);
        println!("Total polygons: {:?}", polygon_dist.iter().sum::<usize>());
        println!(
            "Interaction Energy (epsilon): {:.4}",
            step.interaction_energy()
        );
        println!("Chemical potential (mu): {:.4}", step.chemical_potential());
        println!(
            "Acceptance ratio: {:.4}",
            run_stats.num_accepts() as f64 / run_stats.num_attempts() as f64
        );
        for (idx, shape_stats) in calc_bond_distribution(vmmc).iter().enumerate() {
            let weighted_sum: usize = shape_stats.iter().enumerate().map(|(i, c)| i * c).sum();
            println!(
                "shape_{:?} bond distribution: {:?} average degree = {:.4}",
                idx,
                shape_stats,
                weighted_sum as f64 / vmmc.particles().num_particles() as f64
            );
        }
        let t = Instant::now();
        println!("Execution time: {:?}", t - self.timestamp);
        println!("Total execution time: {:?}", t - self.start_time);
        self.timestamp = t;
    }

    fn state(&self) {}
}
