use crate::num;
use fixed::consts::{PI, TAU};
use fixed::types::I32F32;
use std::fmt;

use crate::types::Num;
use serde::{de, Deserialize, Serialize};
// a = upper_bound
// b = lower_bound
// c = theta

// TODO: extensively test
// Note: assumes theta, target_theta, and tolerance are all in [0, 2pi]
fn in_modular_range(theta: Num, target_theta: Num, tolerance: Num) -> bool {
    let lower_bound = ((target_theta - tolerance) + I32F32::TAU) % I32F32::TAU;
    let upper_bound = (target_theta + tolerance) % I32F32::TAU;
    if lower_bound <= upper_bound {
        theta >= lower_bound && theta <= upper_bound
    } else {
        !(theta < lower_bound && theta > upper_bound)
    }
}

// #[test]
// fn in_modular_range_test(){

// }

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Patch {
    radius: Num,  // radius of patch (in units of particle diameter)
    theta: Num,   // angle in degrees
    chemtype: u8, // patches must have same chemtype to be compatible
}

impl Patch {
    pub fn new(radius: Num, theta: Num, chemtype: u8) -> Self {
        Self {
            radius,
            theta,
            chemtype,
        }
    }

    pub fn radius(&self) -> Num {
        self.radius
    }

    pub fn theta(&self) -> Num {
        self.theta
    }

    pub fn chemtype(&self) -> u8 {
        self.chemtype
    }

    /// TODO: add derivation here
    /// +- radians where it can still interact
    ///  = +- acos 1-(patch_radius^2/2)
    pub fn angle_tolerance(&self) -> Num {
        let pr = self.radius;
        cordic::acos(num!(1.0) - num!(0.5) * pr * pr)
        //(1.0 - (r * r) / (2.0 * PARTICLE_RADIUS * PARTICLE_RADIUS)).acos()
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Morphology {
    patches: Vec<Patch>,
    #[serde(skip_serializing)]
    max_radius: Num,
    // maps theta -> sin(theta) for each patch angle
    #[serde(skip_serializing)]
    sin_theta: Vec<Num>,
    // maps theta -> cos(theta) for each patch angle
    #[serde(skip_serializing)]
    cos_theta: Vec<Num>,
    // TODO: better docs
    // maps patch -> angle tolerance (+- radians where it can still interact)
    #[serde(skip_serializing)]
    angle_tolerances: Vec<Num>,
}

impl<'de> de::Deserialize<'de> for Morphology {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: de::Deserializer<'de>,
    {
        struct MorphologyVisitor;

        impl<'de> de::Visitor<'de> for MorphologyVisitor {
            type Value = Morphology;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("Morphology")
            }

            fn visit_map<V>(self, mut map: V) -> Result<Self::Value, V::Error>
            where
                V: de::MapAccess<'de>,
            {
                if map.next_key::<String>()?.is_some() {
                    let v1: Vec<Patch> = map.next_value()?;
                    Ok(Morphology::new(v1))
                } else {
                    Err(de::Error::missing_field("shapes"))
                }
            }
        }

        deserializer.deserialize_map(MorphologyVisitor {})
    }
}

impl Morphology {
    pub fn new(patches: Vec<Patch>) -> Self {
        let max_radius = patches
            .iter()
            .map(|p| p.radius)
            .fold(Num::MIN, |a, b| a.max(b));

        let mut sin_theta: Vec<_> = Vec::new();
        let mut cos_theta: Vec<_> = Vec::new();
        let mut angle_tolerances: Vec<_> = Vec::new();
        for patch in &patches {
            let theta = patch.theta() * I32F32::PI / num!(180.0); // convert degrees to radians
            sin_theta.push(cordic::sin(theta));
            cos_theta.push(cordic::cos(theta));
            angle_tolerances.push(patch.angle_tolerance());
        }

        // println!("tolerances = {:?}", angle_tolerances);
        Self {
            patches,
            max_radius,
            sin_theta,
            cos_theta,
            angle_tolerances,
        }
    }

    pub fn patches(&self) -> &[Patch] {
        &self.patches
    }

    pub fn max_patch_radius(&self) -> Num {
        self.max_radius
    }

    pub fn sin_theta(&self, patch_idx: usize) -> Num {
        self.sin_theta[patch_idx]
    }

    pub fn cos_theta(&self, patch_idx: usize) -> Num {
        self.cos_theta[patch_idx]
    }

    // TODO: can probably be optimized
    /// given an angle, the patch that contains the angle
    /// if there is no such patch, returns None
    pub fn closest_patch(&self, theta: Num) -> Option<&Patch> {
        // println!("Finding closest patch to {theta}");
        for (p_idx, patch) in self.patches.iter().enumerate() {
            let patch_theta = patch.theta() * I32F32::PI / num!(180.0); // convert degrees to radians
            let angle_tolerance = self.angle_tolerances[p_idx];

            if in_modular_range(theta, patch_theta, angle_tolerance) {
                return Some(patch);
            }
        }
        None
    }

    pub fn regular_3patch(radius: Num) -> Self {
        let p0 = Patch::new(radius, num!(0.0), 0);
        let p1 = Patch::new(radius, num!(120.0), 0);
        let p2 = Patch::new(radius, num!(240.0), 0);
        Self::new(vec![p0, p1, p2])
    }

    pub fn regular_4patch(radius: Num) -> Self {
        let p0 = Patch::new(radius, num!(0.0), 0);
        let p1 = Patch::new(radius, num!(90.0), 0);
        let p2 = Patch::new(radius, num!(180.0), 0);
        let p3 = Patch::new(radius, num!(270.0), 0);
        Self::new(vec![p0, p1, p2, p3])
    }

    pub fn regular_6patch(radius: Num) -> Self {
        let p0 = Patch::new(radius, num!(0.0), 0);
        let p1 = Patch::new(radius, num!(60.0), 0);
        let p2 = Patch::new(radius, num!(120.0), 0);
        let p3 = Patch::new(radius, num!(180.0), 0);
        let p4 = Patch::new(radius, num!(240.0), 0);
        let p5 = Patch::new(radius, num!(300.0), 0);
        Self::new(vec![p0, p1, p2, p3, p4, p5])
    }
}
