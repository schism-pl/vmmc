use chemical_potential::maybe_particle_exchange;
use morphology::{Morphology, Patch};
use position::DimVec;
use protocol::{FixedProtocol, ProtocolStep};
use quickcheck::{Arbitrary, Gen};
use rand::{rngs::SmallRng, SeedableRng};
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

#[derive(Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct InputParams {
    pub num_particles: usize,
    pub protocol: FixedProtocol, // TODO: make more flexible
    pub shapes: Vec<Morphology>,

    pub box_width: f64,
    pub box_height: f64,

    pub prob_translate: f64,
    pub max_translation: f64,
    pub max_rotation: f64,
}

impl Default for InputParams {
    fn default() -> Self {
        let num_particles = 500;
        let box_width = 75.0;
        let box_height = 75.0;

        let prob_translate = 0.5;
        let max_translation = 0.15;
        let max_rotation = 0.2;

        let protocol = FixedProtocol::flat_protocol(0.0, 10.0, 1000);

        let shapes = vec![Morphology::regular_3patch(0.1)];

        Self {
            num_particles,
            protocol,
            shapes,

            box_width,
            box_height,

            prob_translate,
            max_translation,
            max_rotation,
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
    // TODO: check smoothness of protocol?
    pub fn check(&self) {
        // Basic range checks
        assert!(self.num_particles <= 2500);
        assert!(self.box_width >= 10.0 && self.box_height >= 10.0);
        assert!(self.box_width <= 200.0 && self.box_height <= 200.0);
        assert!(
            self.prob_translate >= 0.0
                && self.prob_translate <= 1.0
                && (self.prob_translate.is_normal() || self.prob_translate.is_zero())
        );
        assert!(
            self.max_translation > 0.0
                && self.max_translation <= 1.0
                && self.max_translation.is_normal()
        );
        assert!(
            self.max_rotation > 0.0 && self.max_rotation <= 1.0 && self.max_rotation.is_normal()
        );

        // check well-formedness of protocol
        for step in self.protocol.clone() {
            let mu = step.chemical_potential();
            let epsilon = step.interaction_energy();
            assert!(-20.0 <= mu && mu <= 20.0);
            assert!(0.01 <= epsilon && epsilon <= 20.0);
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
// we need num_cells >= 2*num_particles

// for testing
// TODO: the patch_radius we use is wrong
impl Arbitrary for InputParams {
    fn arbitrary(g: &mut Gen) -> Self {
        let num_particles = usize_in_range(g, 0, 2500);
        let interaction_energy = f64_in_range(g, 0.01, 20.0); // kBT
        let patch_radius = f64_in_range(g, 0.01, 0.2); // diameter of patch (in units of particle diameter)

        let mut box_width = 0.0;
        let mut box_height = 0.0;

        // while less than 1.5 cell per particle, pick a new box raidus
        while box_width * box_height / (1.0 + patch_radius).powi(2) <= 1.5 * num_particles as f64 {
            box_width = f64_in_range(g, 10.0, 200.0);
            box_height = f64_in_range(g, 10.0, 200.0);
        }

        let prob_translate = f64_in_range(g, 0.0, 1.0);
        let max_translation = f64_in_range(g, 0.0, 1.0);
        let max_rotation = f64_in_range(g, 0.0, 1.0);
        let num_sweeps = 10;

        // TODO: make this more varied
        let protocol = FixedProtocol::flat_protocol(0.0, interaction_energy, num_sweeps);

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
            num_particles,
            protocol,
            shapes,
            box_width,
            box_height,

            prob_translate,
            max_translation,
            max_rotation,
        }
    }
}

// angle coverage of patch = +- acos 1-(patch_radius^2/2)
impl Arbitrary for SimBox {
    fn arbitrary(g: &mut Gen) -> Self {
        let ip = InputParams::arbitrary(g);
        let box_dimensions = DimVec::new([ip.box_width, ip.box_height]);
        let patch_radius = 0.05; // TODO: remove this
        let max_interaction_range = 1.0 + patch_radius;

        let mut rng = SmallRng::from_entropy();
        SimBox::new_with_randomized_particles(
            box_dimensions,
            max_interaction_range,
            ip.num_particles,
            ip.shapes,
            &mut rng,
        )
    }
}

pub fn vmmc_from_config(ip: &InputParams, rng: &mut SmallRng) -> vmmc::Vmmc {
    let box_dimensions = DimVec::new([ip.box_width, ip.box_height]);

    let max_interaction_range = 1.0 + ip.max_patch_radius();

    // No initial position, so we will use a randomized start position
    let simbox = SimBox::new_with_randomized_particles(
        box_dimensions,
        max_interaction_range,
        ip.num_particles,
        ip.shapes.clone(),
        rng,
    );

    let params = vmmc::VmmcParams::new(ip.prob_translate, ip.max_translation, ip.max_rotation);
    let interaction_energy = ip.protocol.initial_interaction_energy();
    vmmc::Vmmc::new(simbox, params, interaction_energy)
}

pub trait VmmcCallback {
    type CbResult;
    fn run(&mut self, vmmc: &Vmmc, step: &ProtocolStep, idx: usize, run_stats: &RunStats);
    fn state(&self) -> Self::CbResult;
}

pub fn run_vmmc<Cbr>(
    vmmc: &mut Vmmc,
    protocol: FixedProtocol,
    mut callback: Option<Box<dyn VmmcCallback<CbResult = Cbr>>>,
    rng: &mut SmallRng,
) -> Option<Cbr> {
    for (idx, protocol_step) in protocol.enumerate() {
        let mut run_stats = RunStats::new();
        vmmc.set_interaction_energy(protocol_step.interaction_energy());
        let chemical_potential = protocol_step.chemical_potential();
        for _ in 0..1000 {
            let stats = vmmc.step_n(1000, rng);
            run_stats = stats + run_stats;
            maybe_particle_exchange(vmmc, chemical_potential, rng);
        }
        if let Some(ref mut cb) = callback {
            cb.run(vmmc, &protocol_step, idx, &run_stats);
        }
    }
    callback.map(|cb| cb.state())
}
