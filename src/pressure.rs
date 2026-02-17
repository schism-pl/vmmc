use rand::Rng;

use crate::{vmmc::Vmmc, Prng};

const EPS_LOGV: f64 = 0.02;

pub fn maybe_volume_change(
    vmmc: &mut Vmmc,
    target_volume_x: f64,
    target_volume_y: f64,
    rng: &mut Prng,
) {
    // 1) propose a volume change
    let proposed_volume_x = propose_log_volume_move(rng, target_volume_x).unwrap();
    let proposed_volume_y = propose_log_volume_move(rng, target_volume_y).unwrap();
    log::info!(
        "Proposed volume change: x: {} -> {}, y: {} -> {}",
        target_volume_x,
        proposed_volume_x,
        target_volume_y,
        proposed_volume_y
    );
    // 2) scale particles + box
    // 3) recompute energy + accept/reject
    let _ = vmmc;
    let _ = target_volume_x;
    let _ = target_volume_y;
    let _ = rng;
    todo!()
}

/// Propose a log-volume move:
///   delta ~ Uniform(-eps_logv, +eps_logv)
///   V_new = V_old * exp(delta)
///
/// Returns V_new.
pub fn propose_log_volume_move(rng: &mut Prng, v_old: f64) -> Result<f64, String> {
    if !v_old.is_finite() || v_old <= 0.0 {
        return Err(format!("v_old must be finite and > 0, got {}", v_old));
    }
    debug_assert!(EPS_LOGV.is_finite() && EPS_LOGV > 0.0);

    // delta is in the range [-eps_logv, +eps_logv]
    let delta = rng.gen_range(-EPS_LOGV..=EPS_LOGV);

    // V_new = V_old * exp(delta)
    let v_new = v_old * delta.exp();

    // Numerically defensive: exp(delta) is always > 0, but check finiteness.
    if !v_new.is_finite() || v_new <= 0.0 {
        return Err(format!(
            "proposed v_new became non-finite or non-positive (v_old={}, delta={}, v_new={})",
            v_old, delta, v_new
        ));
    }

    Ok(v_new)
}
