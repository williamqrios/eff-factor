use std::{fmt::Display}; 
enum ReactionType {
    TypeI, // A + B <-> C + D
    TypeII, // 2A <-> C + D
    TypeIII, // A + B <-> 2 C 
    TypeIV, // A <-> C + D 
    TypeV, // A + B <-> C 
    TypeVI, // A <-> C
    TypeVII // 2 A + B <-> C + D 
}

enum GeometryType {
    Slab { length: f64 },
    Sphere { radius: f64 },
    Cylinder { radius: f64 }, 
    Other {surf: f64, vol: f64 },
}

struct CatalystData {
    geometry: GeometryType,
    epsilon: f64, // Porosity 
    tau: f64, // Tortuosity 
    rho_p: f64, // Density 
}


impl CatalystData {
    fn new(geometry: GeometryType, epsilon: f64, tau: Option<f64>, rho_p: f64) -> Self {
        // Estimate tortuosity from porosity as the inverse of the porosity
        let tau = tau.unwrap_or(epsilon.powf(-1.0)); 
        Self { geometry, epsilon, tau, rho_p }
    }
}

impl Display for ReactionType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReactionType::TypeI => write!(f, "A + B <-> C + D\nr = k (Ca Cb - 1/Kc Cc Cd)"),
            ReactionType::TypeII => write!(f, "2 A <-> C + D\nr = k (Ca^2 - 1/Kc Cc Cd)"),
            ReactionType::TypeIII => write!(f, "A + B <-> 2 C\nr = k (Ca Cb - 1/Kc Cc^2)"),
            ReactionType::TypeIV => write!(f, "A <-> C + D\nr = k (Ca - 1/Kc Cc Cd)"),
            ReactionType::TypeV => write!(f, "A + B <-> C\nr = k (Ca Cb - 1/Kc Cc)"), 
            ReactionType::TypeVI => write!(f, "A <-> C\nr = k (Ca - 1/Kc Cc)"),
            ReactionType::TypeVII => write!(f, "2 A + B <-> C + D\nr = k (Ca Cb - 1/Kc Cc Cd/Ca)"), 
        }
    }
}

/// I_0 function
fn bessel_i0(z: f64) -> f64 {
    let tol = 1e-12; 
    let mut series_sum = 0.0; 
    let mut k = 0; 
    loop {
        let factorial: u128 = if k >= 1 {
            (1..=k).product()
        } else { 1 };
        let current_term = f64::powf(z/2.0, 2.0 * k as f64) / (factorial as f64 * factorial as f64);
        k += 1; 
        series_sum += current_term; 
        if current_term < tol {
            break; 
        }
    }
    series_sum 
}

fn bessel_i1(z: f64) -> f64 {
    let tol = 1e-12; 
    let mut series_sum = 0.0; 
    let mut k = 0; 
    loop {
        let factorial_k: u128 = if k >= 1 {
            (1..=k).product()
        } else { 1 };
        let factorial_k_plus_1 = factorial_k * (k + 1); 
        let current_term = f64::powf(z/2.0, 2.0 * k as f64) / (factorial_k as f64 * factorial_k_plus_1 as f64);
        k += 1; 
        series_sum += current_term; 
        if current_term < tol {
            break; 
        }
    }
    series_sum * 0.5 * z 
}

/// Regular thiele modulus
fn thiele_modulus(catalyst: &CatalystData, k: f64, deff_mix_a: f64) -> f64 {
    // phi = l * sqrt ( k * rho_p/Deff,A ) 
    let deff_cat_a = deff_mix_a * catalyst.epsilon/catalyst.tau;
    // characteristic length  
    let ell = match catalyst.geometry {
        GeometryType::Other { surf, vol } => vol / surf, 
        GeometryType::Slab { length } => length, 
        GeometryType::Sphere { radius } => radius, 
        GeometryType::Cylinder { radius } => radius, 
    }; 
    let phi = ell * f64::powf( k * catalyst.rho_p/deff_cat_a, 0.5 ); 
    phi 
}

fn factor_0(deff_mix: &[f64; 4], conc: &[f64; 4], reaction: &ReactionType) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => 0.0,
        TypeII => 0.0, 
        TypeIII => 0.0, 
        TypeIV => 0.0, 
        TypeV => 0.0, 
        TypeVI => 0.0, 
        TypeVII => {
            conc[0] * conc[0] * (conc[2]/conc[0] + deff_mix[0]/(2.0 * deff_mix[2]) ) * (conc[3]/conc[0] + deff_mix[0]/(2.0 * deff_mix[3]))
        }
    }
}

fn factor_1(deff_mix: &[f64; 4], conc: &[f64; 4], reaction: &ReactionType) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            conc[0] * conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        },
        TypeII => {
            conc[0] * conc[0] * (conc[2]/conc[0] + 0.5 * deff_mix[0]/deff_mix[2]) * (conc[3]/conc[0] + 0.5 * deff_mix[0]/deff_mix[3])
        }, 
        TypeIII => {
            conc[0] * conc[0] * f64::powf(conc[2]/conc[0] + 2.0 * deff_mix[0]/deff_mix[2], 2.0)
        }, 
        TypeIV => {
            conc[0] * conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        }, 
        TypeV => {
            conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2])
        }, 
        TypeVI => {
            conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2])
        }, 
        TypeVII => {
            - 0.5 * deff_mix[0]/deff_mix[3] * conc[0] * (conc[2]/conc[0] + 0.5 * deff_mix[0]/deff_mix[2]) - conc[0] * 0.5 * deff_mix[0]/deff_mix[2] * (conc[3]/conc[0] + 0.5 * deff_mix[0]/deff_mix[3])
        }
    }
}

fn factor_2(deff_mix: &[f64; 4], conc: &[f64; 4], reaction: &ReactionType, Kc: f64) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + conc[0]/Kc * deff_mix[0]/deff_mix[3] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) + conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        },
        TypeII => {
            conc[0]/Kc * deff_mix[0]/(2.0 * deff_mix[3]) * (conc[2]/conc[0] + deff_mix[0]/(2.0 * deff_mix[2])) + conc[0]/Kc * deff_mix[0]/(2.0 * deff_mix[2]) * (conc[3]/conc[0] + deff_mix[0]/(2.0 * deff_mix[3]))
        }, 
        TypeIII => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + 4.0 * conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[2]/conc[0] + 2.0 * deff_mix[0]/deff_mix[2])
        }, 
        TypeIV => {
            1.0 + conc[0]/Kc * deff_mix[0]/deff_mix[3] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) + conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        }, 
        TypeV => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + 1.0/Kc * deff_mix[0]/deff_mix[2]
        }, 
        TypeVI => {
            1.0 + 1.0/Kc * deff_mix[0]/deff_mix[2]
        }, 
        TypeVII => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/(2.0 * deff_mix[1])) - 1.0/Kc * deff_mix[0]/(2.0 * deff_mix[2]) * deff_mix[0]/(2.0 * deff_mix[3]) 
        }
    }
}

fn factor_3(deff_mix: &[f64; 4], reaction: &ReactionType, Kc: f64) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            deff_mix[0]/deff_mix[1] - 1.0/Kc * deff_mix[0]/deff_mix[2] * deff_mix[0]/deff_mix[3]
        },
        TypeII => {
            1.0 - 1.0/Kc * deff_mix[0]/(2.0 * deff_mix[2]) * deff_mix[0]/(2.0 * deff_mix[3])
        }, 
        TypeIII => {
            deff_mix[0]/deff_mix[1] - 1.0/Kc * f64::powf(2.0 * deff_mix[0]/deff_mix[2], 2.0)
        }, 
        TypeIV => {
            -1.0/Kc * deff_mix[0]/deff_mix[2] * deff_mix[0]/deff_mix[3]
        }, 
        TypeV => deff_mix[0]/deff_mix[1], 
        TypeVI => 0.0, 
        TypeVII => {
            deff_mix[0]/(2.0 * deff_mix[1])
        }
    }
}


/// Analytical Equations
fn generalized_thiele_modulus(conc: &[f64; 4], deff_mix: &[f64; 4], k: f64, Kc: f64, Ca_eq: f64, catalyst: &CatalystData, reaction: &ReactionType) -> f64 {
    // Regular Thiele modulus
    let phi = thiele_modulus(catalyst, k, deff_mix[0]); 
    use ReactionType::*; 
    // Reaction rate at catalyst surface 
    let rs = match reaction {
        TypeI => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[3] )
        },
        TypeII => {
            k * ( conc[0] * conc[0] - 1.0/Kc * conc[2] * conc[3] )
        },
        TypeIII => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[2] )
        },
        TypeIV => {
            k * ( conc[0] - 1.0/Kc * conc[2] * conc[3] )
        }, 
        TypeV => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] )
        },
        TypeVI => {
            k * ( conc[0] - 1.0/Kc * conc[2] )
        },
        TypeVII => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[3]/conc[0] )
        },
    }; 
    let constant: f64 = factor_0(deff_mix, conc, reaction) * (conc[0].ln() - Ca_eq.ln())/(-Kc) + factor_1(deff_mix, conc, reaction) * (conc[0] - Ca_eq)/(-Kc) + factor_2(deff_mix, conc, reaction, Kc) * ( conc[0].powf(2.0) - Ca_eq.powf(2.0) )/2.0 + factor_3(deff_mix, reaction, Kc) * (conc[0].powf(3.0) - Ca_eq.powf(3.0) )/3.0; 
    let phi_g = phi * rs /( f64::powf(2.0 * constant, 0.5) * k );
    phi_g
}


/// Computes the effectiveness factor given a generalized Thiele modulus and geometry  
fn effectiveness_factor(phi_g: f64, geometry: &GeometryType) -> f64 {
    match geometry {
        GeometryType::Slab { length: _ } => {
            // eta = tanh phi_g / phi_g
            phi_g.tanh() / phi_g
        },
        GeometryType::Sphere { radius: _ } => {
            // eta = 3 / phi_g * (1 / tanh phi_g - 1 / phi_g )
            3.0 / phi_g * (1.0 / phi_g.tanh() - 1.0 / phi_g )
        }, 
        GeometryType::Cylinder { radius: _ } => { 
            // eta = 2 / phi_g * I_1 (phi_g) / I_0 (phi_g)
            2.0 / phi_g * bessel_i1(phi_g) / bessel_i0(phi_g)
        }, 
        GeometryType::Other { surf: _, vol: _ } => { 
            // use the generic equation
            // eta = tanh phi_g / phi_g
            phi_g.tanh() / phi_g
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bessel_i0_small_values() {
        // Test for small values of z
        assert!((bessel_i0(0.0) - 1.0).abs() < 1e-12);
        assert!((bessel_i0(0.1) - 1.0025).abs() < 1e-4);
    }

    #[test]
    fn test_bessel_i0_large_values() {
        // Test for larger values of z
        assert!((bessel_i0(1.0) - 1.2660658777520082).abs() < 1e-12);
        assert!((bessel_i0(2.0) - 2.2795853023360673).abs() < 1e-12);
    }

    #[test]
    fn test_bessel_i0_edge_case() {
        // Test for edge case (only works with u128)
        assert!((bessel_i0(10.0) - 2815.716628466254).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i1_small_values() {
        // Test for small values of z
        assert!((bessel_i1(0.0) - 0.0).abs() < 1e-12);
        assert!((bessel_i1(0.1) - 0.0500625).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i1_large_values() {
        // Test for larger values of z
        assert!((bessel_i1(1.0) - 0.565159103992485).abs() < 1e-12);
        assert!((bessel_i1(2.0) - 1.590636854637329).abs() < 1e-12);
    }

    #[test]
    fn test_bessel_i1_edge_case() {
        // Test for edge case (only works with u128)
        assert!((bessel_i1(10.0) - 2670.988303701255).abs() < 1e-6);
    }

    #[test]
    fn test_phi() {
        let catalyst = CatalystData::new(GeometryType::Slab { length: 0.05 }, 0.5, Some(0.5), 1.0);
        let k = 4.0; 
        let deff_mix_a = 5.0e-4; 
        let phi = thiele_modulus(&catalyst, k, deff_mix_a); 
        dbg!(phi);
        // assert_eq!(1, 0)
    }

    #[test]
    fn test_phi_g() {
        let catalyst = CatalystData::new(GeometryType::Sphere { radius: 0.744/2.0 * 1.0e-2 }, 0.4886, Some(1.3), 600.0);
        let reaction = ReactionType::TypeI;
        let k = 4.5747e-5; // dm6 mol-1 min-1 gcat-1 
        let Kc = 2.6700;
        let deff_mix = [0.3172e-4, 0.2203e-4, 0.2009e-4, 0.3684e-4]; // dm2/min
        let conc = [8.5326, 8.5326, 0.0, 0.0];
        let Ca_eq = 3.2769;
        // Note: Ca_eq should be smaller than Ca0
        let phi_g = generalized_thiele_modulus(&conc, &deff_mix, k, Kc, Ca_eq, &catalyst, &reaction);
        let eta = effectiveness_factor(phi_g, &catalyst.geometry);
        assert!((eta - 0.9626).abs() < 1e-3)
    }
}