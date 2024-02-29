use std::{f64::consts::PI, fmt};

use serde::{de, Deserialize, Deserializer, Serialize};

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
}

// impl Serialize for Equation {
//     fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
//     where
//         S: Serializer,
//     {
//         serializer.serialize_str(&self.to_string())
//     }
// }

// impl<'de> Deserialize<'de> for Equation {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {
//         let s = String::deserialize(deserializer)?;
//         FromStr::from_str(&s).map_err(de::Error::custom)
//     }
// }

#[derive(Debug, Deserialize)]
struct Inner {
    patches: Vec<Patch>,
}

// impl<'de> Deserialize<'de> for Morphology {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {

//         let patches = Vec::<Patch>::deserialize(deserializer)?;
//         Ok(Morphology::new(patches))
//     }
// }

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
                // let patches: Inner = Vec::<Patch>::deserialize(deserializer)?;
                // Ok(Morphology::new(patches))

                if let Some(k1) = map.next_key::<Inner>()? {
                    println!("k1 = {:?}", k1);
                    let v1: Inner = map.next_value()?;
                    println!("v1 = {:?}", v1);
                    if let Some(k2) = map.next_key::<Inner>()? {
                        println!("k2 = {:?}", k2);
                        let value: Inner = map.next_value()?;
                        Ok(Morphology::new(value.patches))
                    } else {
                        Err(de::Error::missing_field("shapes"))
                    }
                } else {
                    unimplemented!()
                }

                // if let Some(key) = map.next_key()? {
                //     let value: Inner = map.next_value()?;
                //     if let Some(_) = map.next_key::<&str>()? {
                //         Err(de::Error::duplicate_field("name"))
                //     } else {
                //         Ok(Self::Value {
                //             name: key,
                //             topics: value.topics,
                //         })
                //     }
                // } else {
                //     Err(de::Error::missing_field("name"))
                // }
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
        for patch in &patches {
            let theta = patch.theta() * PI / 180.0; // convert degrees to radians
            sin_theta.push(theta.sin());
            cos_theta.push(theta.cos());
        }

        Self {
            patches,
            max_radius,
            sin_theta,
            cos_theta,
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
