use rand::{rngs::SmallRng, Rng};

use crate::{consts::MAX_PARTICLES, vmmc::Vmmc};

fn maybe_remove_particle(vmmc: &mut Vmmc, chemical_potential: f64, rng: &mut SmallRng) {
    let num_particles = vmmc.particles().num_particles();
    if num_particles == 0 {
        return;
    }
    // P(remove) = (N+1)/V * e^(E - mu)
    // Note: E will be negative and mu can be positive, can be negative
    // more positive means more crowded substrate
    let p_id = vmmc.choose_random_p_id(rng);
    let p = vmmc.particle(p_id);
    let energy = vmmc.get_particle_energy(p);
    let energy_factor = energy - chemical_potential;
    let p_accept = (num_particles as f64 + 1.0) / vmmc.simbox().volume() * f64::exp(energy_factor);
    if rng.gen::<f64>() < p_accept {
        vmmc.simbox_mut().remove_particle(p_id);
    }
}

fn maybe_insert_particle(vmmc: &mut Vmmc, chemical_potential: f64, rng: &mut SmallRng) {
    let num_particles = vmmc.particles().num_particles();

    if num_particles >= MAX_PARTICLES {
        return;
    }
    // P(insert) = V/(N+2) * e^(mu - E)
    // Note: E will be negative and mu can be positive, can be negative
    // more positive means more crowded substrate

    let p = vmmc.simbox_mut().new_random_particle(rng);
    let energy = vmmc.get_particle_energy(&p);
    let energy_factor = energy + chemical_potential;
    let p_accept = vmmc.simbox().volume() / (num_particles as f64 + 2.0) * f64::exp(energy_factor);
    if rng.gen::<f64>() < p_accept {
        vmmc.simbox_mut().insert_particle(p);
    } else {
        // we aren't using this particle, so its id is up for grabs
        vmmc.simbox_mut().particles_mut().push_unused_p_id(p.id());
    }
}

fn particle_exchange(vmmc: &mut Vmmc, chemical_potential: f64, rng: &mut SmallRng) {
    if rng.gen::<f64>() < 0.5 {
        maybe_remove_particle(vmmc, chemical_potential, rng);
    } else {
        maybe_insert_particle(vmmc, chemical_potential, rng);
    }
}

// TODO: better name?
// equations G1-G3 in https://journals.aps.org/prx/pdf/10.1103/PhysRevX.4.011044
// TODO: this equation satisfies detailed balance, but I'm not quite sure how/why
pub fn maybe_particle_exchange(vmmc: &mut Vmmc, chemical_potential: f64, rng: &mut SmallRng) {
    let p_exchange = 1.0 / (1.0 + vmmc.particles().num_particles() as f64);
    if rng.gen::<f64>() < p_exchange {
        particle_exchange(vmmc, chemical_potential, rng);
    }
}
