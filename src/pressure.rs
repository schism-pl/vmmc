use rand::Rng;

use crate::{vmmc::Vmmc, Prng};

// Hyperparameter for log-volume moves: maximum magnitude of log-volume change per move.
const EPS_LOGV: f64 = 0.005;
// const GAMMA = 0.02; // step size for axis moves toward target

pub fn maybe_volume_change(
    vmmc: &mut Vmmc,
    target_x: Option<f64>,
    target_y: Option<f64>,
    rng: &mut Prng,
) {
    // Tune how often we try to volume change
    if rng.random_range(0_u32..100_000_u32) != 0 {
        return;
    }
    // Unpack variables
    let changing_x = target_x.is_some();
    let changing_y = target_y.is_some();
    let p_x = target_x.unwrap_or(0.0);
    let p_y = target_y.unwrap_or(0.0);

    // let target_volume = target_x * target_y;
    // let v_old = vmmc.simbox().volume();
    let x_old = vmmc.simbox().x();
    let y_old = vmmc.simbox().y();

    // 1) propose a volume change
    let x_new = if changing_x {
        propose_axis_move(rng, x_old).unwrap()
    } else {
        vmmc.simbox().x()
    };
    let y_new = if changing_y {
        propose_axis_move(rng, y_old).unwrap()
    } else {
        vmmc.simbox().y()
    };
    // let v_new = x_new + y_new;

    log::info!(
        "proposed volume change: x {:.4} -> {:.4}, y {:.4} -> {:.4}",
        x_old,
        x_new,
        y_old,
        y_new
    );
    // 2) scale particles + box
    let proposed_vmmc = match vmmc.rescaled_simbox(x_new, y_new) {
        Ok(sim) => sim,
        Err(_) => {println!("Failed to rescale simbox"); return},
    };

    // 3) recompute energy + accept/reject
    // 3a) compute total potential energies for old and proposed states
    let old_energy = vmmc.get_total_energy();
    let new_energy = proposed_vmmc.get_total_energy();
    // println!("old energy: {:.4}, new energy: {:.4}", old_energy, new_energy);

    // 3b) compute external work increment for anisotropic update convention
    let d_u = -(new_energy - old_energy);
    let d_w = p_x * 20.0 * y_old * (x_new - x_old) + p_y * x_new * (y_new - y_old);



    // 3c) compute Jacobian term from affine scaling (area change)
    let a_old = x_old * y_old;
    let a_new = x_new * y_new;
    if a_old <= 0.0 || a_new <= 0.0 {
        {println!("Invalid area: a_old = {:.4}, a_new = {:.4}", a_old, a_new); return};
    }
    let log_j = (vmmc.particles().num_particles() as f64) * (a_new / a_old).ln();

    // 3d) Metropolis-Hastings accept/reject in log space (beta = 1)
    // energies are already in units of kBT, so beta = 1.0
    let log_alpha = -(d_u + d_w) + log_j;
        log::info!(
        "energy change: ΔU = {:.4}, ΔW = {:.4} log_alpha = {:.4})",
        d_u, d_w, log_alpha
    );
    let u = rng.random_range(0.0..1.0);
    let u = if u <= 0.0 { f64::MIN_POSITIVE } else { u };
    log::info!(
        "volume change energy: ΔU = {:.4}, ΔW = {:.4}, log_J = {:.4}, log_alpha = {:.4}, u = {:.4}",
        d_u, d_w, log_j, log_alpha, u
    );
    if u.ln() < log_alpha.min(0.0) {
        // 3e) accept: commit proposed box and scaled positions
        *vmmc.simbox_mut() = proposed_vmmc.simbox().clone();
    }
}

/// Propose a new axis length by a random log-step:
///   delta ~ Uniform(-eps_logv, +eps_logv)
///   L_new = L_old * exp(delta)
pub fn propose_axis_move(rng: &mut Prng, l_old: f64) -> Result<f64, String> {
    if !l_old.is_finite() || l_old <= 0.0 {
        return Err(format!("l_old must be finite and > 0, got {}", l_old));
    }
    debug_assert!(EPS_LOGV.is_finite() && EPS_LOGV > 0.0);

    let delta = rng.random_range(-EPS_LOGV..=EPS_LOGV);
    let l_new = l_old * delta.exp();

    if !l_new.is_finite() || l_new <= 0.0 {
        return Err(format!(
            "proposed l_new became non-finite or non-positive (l_old={}, delta={}, l_new={})",
            l_old, delta, l_new
        ));
    }

    Ok(l_new)
}
