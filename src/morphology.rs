use std::{f64::consts::PI, fmt};

use serde::{de, Deserialize, Serialize};

use crate::consts::PARTICLE_RADIUS;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Patch {
    radius: f64,  // radius of patch (in units of particle diameter)
    theta: f64,   // angle in degrees
    chemtype: u8, // patches must have same chemtype to be compatible
}

impl Patch {
    pub fn new(radius: f64, theta: f64, chemtype: u8) -> Self {
        Self {
            radius,
            theta,
            chemtype,
        }
    }

    pub fn radius(&self) -> f64 {
        self.radius
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    pub fn chemtype(&self) -> u8 {
        self.chemtype
    }

    /// TODO: add derivation here
    /// +- radians where it can still interact
    ///  = +- acos 1-(patch_radius^2/2)
    pub fn angle_tolerance(&self) -> f64 {
        let r = self.radius;
        1.0 - (r * r) / (2.0 * PARTICLE_RADIUS * PARTICLE_RADIUS)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct Morphology {
    patches: Vec<Patch>,
    #[serde(skip_serializing)]
    max_radius: f64,
    // maps theta -> sin(theta) for each patch angle
    #[serde(skip_serializing)]
    sin_theta: Vec<f64>,
    // maps theta -> cos(theta) for each patch angle
    #[serde(skip_serializing)]
    cos_theta: Vec<f64>,
    // TODO: better docs
    // maps patch -> angle tolerance (+- radians where it can still interact)
    #[serde(skip_serializing)]
    angle_tolerances: Vec<f64>,
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
                if let Some(_) = map.next_key::<String>()? {
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
            .fold(f64::MIN, |a, b| a.max(b));

        let mut sin_theta: Vec<_> = Vec::new();
        let mut cos_theta: Vec<_> = Vec::new();
        let mut angle_tolerances: Vec<_> = Vec::new();
        for patch in &patches {
            let theta = patch.theta() * PI / 180.0; // convert degrees to radians
            sin_theta.push(theta.sin());
            cos_theta.push(theta.cos());
            angle_tolerances.push(patch.angle_tolerance());
        }

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

    pub fn max_patch_radius(&self) -> f64 {
        self.max_radius
    }

    pub fn sin_theta(&self, patch_idx: usize) -> f64 {
        self.sin_theta[patch_idx]
    }

    pub fn cos_theta(&self, patch_idx: usize) -> f64 {
        self.cos_theta[patch_idx]
    }

    // TODO: can probably be optimized
    /// given an angle, the patch that contains the angle
    /// if there is no such patch, returns None
    pub fn closest_patch(&self, theta: f64) -> Option<&Patch> {
        for (p_idx, patch) in self.patches.iter().enumerate() {
            let patch_theta = patch.theta() * PI / 180.0; // convert degrees to radians
                                                          // TODO: use in_range?
            let angle_tolerance = self.angle_tolerances[p_idx];
            if theta >= patch_theta - angle_tolerance && theta < patch_theta + angle_tolerance {
                // TODO: wrap-around?
                return Some(patch);
            }
        }
        None
    }

    pub fn regular_3patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 120.0, 0);
        let p2 = Patch::new(radius, 240.0, 0);
        Self::new(vec![p0, p1, p2])
    }

    pub fn regular_4patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 90.0, 0);
        let p2 = Patch::new(radius, 180.0, 0);
        let p3 = Patch::new(radius, 270.0, 0);
        Self::new(vec![p0, p1, p2, p3])
    }

    pub fn regular_6patch(radius: f64) -> Self {
        let p0 = Patch::new(radius, 0.0, 0);
        let p1 = Patch::new(radius, 60.0, 0);
        let p2 = Patch::new(radius, 120.0, 0);
        let p3 = Patch::new(radius, 180.0, 0);
        let p4 = Patch::new(radius, 240.0, 0);
        let p5 = Patch::new(radius, 300.0, 0);
        Self::new(vec![p0, p1, p2, p3, p4, p5])
    }
}
