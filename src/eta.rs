use crate::models::*; 
use crate::utils::{bessel_i0, bessel_i1};

/// Regular thiele modulus
fn thiele_modulus(catalyst: &CatalystData, k: f64, deff_mix_a: f64) -> f64 {
    // phi = l * sqrt ( k * rho_p/Deff,A ) 
    let deff_cat_a = deff_mix_a * catalyst.epsilon/catalyst.tau;
    
    let phi = catalyst.geometry.get_value() * f64::powf( k * catalyst.rho_p/deff_cat_a, 0.5 ); 
    phi 
}

fn factor_0(deff_mix: &[f64], conc: &[f64], reaction: &ReactionType) -> f64 {
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

fn factor_1(deff_mix: &[f64], conc: &[f64], reaction: &ReactionType) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            conc[0] * conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        },
        TypeII => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd 
            conc[0] * conc[0] * (conc[1]/conc[0] + 0.5 * deff_mix[0]/deff_mix[1]) * (conc[2]/conc[0] + 0.5 * deff_mix[0]/deff_mix[2])
        }, 
        TypeIII => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc 
            conc[0] * conc[0] * f64::powf(conc[2]/conc[0] + 2.0 * deff_mix[0]/deff_mix[2], 2.0)
        }, 
        TypeIV => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd 
            conc[0] * conc[0] * (conc[1]/conc[0] + deff_mix[0]/deff_mix[1]) * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2])
        }, 
        TypeV => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc
            conc[0] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2])
        }, 
        TypeVI => {
            // Legend: conc[0] = Ca, conc[1] = Cb
            conc[0] * (conc[1]/conc[0] + deff_mix[0]/deff_mix[1])
        }, 
        TypeVII => {
            - 0.5 * deff_mix[0]/deff_mix[3] * conc[0] * (conc[2]/conc[0] + 0.5 * deff_mix[0]/deff_mix[2]) - conc[0] * 0.5 * deff_mix[0]/deff_mix[2] * (conc[3]/conc[0] + 0.5 * deff_mix[0]/deff_mix[3])
        }
    }
}

fn factor_2(deff_mix: &[f64], conc: &[f64], reaction: &ReactionType, Kc: f64) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + conc[0]/Kc * deff_mix[0]/deff_mix[3] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2]) + conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[3]/conc[0] + deff_mix[0]/deff_mix[3])
        },
        TypeII => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd 
            conc[0]/Kc * deff_mix[0]/(2.0 * deff_mix[2]) * (conc[1]/conc[0] + deff_mix[0]/(2.0 * deff_mix[1])) + conc[0]/Kc * deff_mix[0]/(2.0 * deff_mix[1]) * (conc[2]/conc[0] + deff_mix[0]/(2.0 * deff_mix[2]))
        }, 
        TypeIII => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc 
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + 4.0 * conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[2]/conc[0] + 2.0 * deff_mix[0]/deff_mix[2])
        }, 
        TypeIV => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd
            1.0 + conc[0]/Kc * deff_mix[0]/deff_mix[2] * (conc[1]/conc[0] + deff_mix[0]/deff_mix[1]) + conc[0]/Kc * deff_mix[0]/deff_mix[1] * (conc[2]/conc[0] + deff_mix[0]/deff_mix[2])
        }, 
        TypeV => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/deff_mix[1]) + 1.0/Kc * deff_mix[0]/deff_mix[2]
        }, 
        TypeVI => {
            // Legend: conc[0] = Ca, conc[1] = Cb
            1.0 + 1.0/Kc * deff_mix[0]/deff_mix[1]
        }, 
        TypeVII => {
            conc[0] * (conc[1]/conc[0] - deff_mix[0]/(2.0 * deff_mix[1])) - 1.0/Kc * deff_mix[0]/(2.0 * deff_mix[2]) * deff_mix[0]/(2.0 * deff_mix[3]) 
        }
    }
}

fn factor_3(deff_mix: &[f64], reaction: &ReactionType, Kc: f64) -> f64 {
    use ReactionType::*;
    match reaction {
        TypeI => {
            deff_mix[0]/deff_mix[1] - 1.0/Kc * deff_mix[0]/deff_mix[2] * deff_mix[0]/deff_mix[3]
        },
        TypeII => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd 
            1.0 - 1.0/Kc * deff_mix[0]/(2.0 * deff_mix[1]) * deff_mix[0]/(2.0 * deff_mix[2])
        }, 
        TypeIII => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc 
            deff_mix[0]/deff_mix[1] - 1.0/Kc * f64::powf(2.0 * deff_mix[0]/deff_mix[2], 2.0)
        }, 
        TypeIV => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd
            -1.0/Kc * deff_mix[0]/deff_mix[1] * deff_mix[0]/deff_mix[2]
        }, 
        TypeV => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc
            deff_mix[0]/deff_mix[1]
        }, 
        TypeVI =>{
            // Legend: conc[0] = Ca, conc[1] = Cb
            0.0
        }, 
        TypeVII => {
            deff_mix[0]/(2.0 * deff_mix[1])
        }
    }
}


/// Analytical Equations
pub fn generalized_thiele_modulus(rx: &ReactionData, catalyst: &CatalystData) -> f64 {
    // Regular Thiele modulus
    let deff_mix = &rx.deff; 
    let phi = thiele_modulus(catalyst, rx.k, deff_mix[0]); 
    let conc = &rx.conc; 
    let k = rx.k;
    let Kc = rx.Kc;
    let Ca_eq = rx.Ca_eq;
    let reaction = &rx.reaction_type;
    use ReactionType::*; 
    // Reaction rate at catalyst surface 
    let rs = match rx.reaction_type {
        TypeI => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[3] )
        },
        TypeII => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd 
            k * ( conc[0] * conc[0] - 1.0/Kc * conc[1] * conc[2] )
        },
        TypeIII => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc 
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[2] )
        },
        TypeIV => {
            // Legend: conc[0] = Ca, conc[1] = Cc, conc[2] = Cd
            k * ( conc[0] - 1.0/Kc * conc[1] * conc[2] )
        }, 
        TypeV => {
            // Legend: conc[0] = Ca, conc[1] = Cb, conc[2] = Cc
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] )
        },
        TypeVI => {
            // Legend: conc[0] = Ca, conc[1] = Cb
            k * ( conc[0] - 1.0/Kc * conc[1] )
        },
        TypeVII => {
            k * ( conc[0] * conc[1] - 1.0/Kc * conc[2] * conc[3]/conc[0] )
        },
    }; 

    let f0 = factor_0(deff_mix, conc, reaction);
    let f1 = factor_1(deff_mix, conc, reaction);
    let f2 = factor_2(deff_mix, conc, reaction, Kc);
    let f3 = factor_3(deff_mix, reaction, Kc);

    let constant: f64 = f0 * ((conc[0]/Ca_eq).ln())/(-Kc) + f1 * (conc[0] - Ca_eq)/(-Kc) + f2 * ( conc[0].powf(2.0) - Ca_eq.powf(2.0) )/2.0 + f3 * (conc[0].powf(3.0) - Ca_eq.powf(3.0) )/3.0; 
    
    let phi_g = phi * rs /( f64::powf(2.0 * constant, 0.5) * k );
    phi_g
}


/// Computes the effectiveness factor given a generalized Thiele modulus and geometry  
pub fn eff_factor(phi_g: f64, geometry: &GeometryType) -> f64 {
    match geometry {
        GeometryType::Slab(_) => {
            phi_g.tanh() / phi_g
        },
        GeometryType::Sphere(_) => {
            3.0 / phi_g * (1.0 / phi_g.tanh() - 1.0 / phi_g )
        }, 
        GeometryType::Cylinder(_) => { 
            2.0 / phi_g * bessel_i1(phi_g) / bessel_i0(phi_g)
        }, 
        GeometryType::Other(_) => { 
            phi_g.tanh() / phi_g
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[ignore]
    #[test]
    fn test_bessel_i0_small_values() {
        // Test for small values of z
        assert!((bessel_i0(0.0) - 1.0).abs() < 1e-12);
        assert!((bessel_i0(0.1) - 1.0025).abs() < 1e-4);
    }
    #[ignore]
    #[test]
    fn test_bessel_i0_large_values() {
        // Test for larger values of z
        assert!((bessel_i0(1.0) - 1.2660658777520082).abs() < 1e-12);
        assert!((bessel_i0(2.0) - 2.2795853023360673).abs() < 1e-12);
    }
    #[ignore]
    #[test]
    fn test_bessel_i0_edge_case() {
        // Test for edge case (only works with u128)
        assert!((bessel_i0(10.0) - 2815.716628466254).abs() < 1e-6);
    }
    #[ignore]
    #[test]
    fn test_bessel_i1_small_values() {
        // Test for small values of z
        assert!((bessel_i1(0.0) - 0.0).abs() < 1e-12);
        assert!((bessel_i1(0.1) - 0.0500625).abs() < 1e-6);
    }
    #[ignore]
    #[test]
    fn test_bessel_i1_large_values() {
        // Test for larger values of z
        assert!((bessel_i1(1.0) - 0.565159103992485).abs() < 1e-12);
        assert!((bessel_i1(2.0) - 1.590636854637329).abs() < 1e-12);
    }
    // #[ignore]
    #[test]
    fn test_bessel_i1_edge_case() {
        // Test for edge case (only works with u128)
        assert!((bessel_i1(10.0) - 2670.988303701255).abs() < 1e-6);
    }
    #[ignore]
    #[test]
    fn test_phi() {
        let catalyst = CatalystData::new(GeometryType::Slab(0.5), 0.5, Some(0.5), 1.0);
        let k = 4.0; 
        let deff_mix_a = 5.0e-4; 
        let phi = thiele_modulus(&catalyst, k, deff_mix_a); 
        // dbg!(phi);
        // assert_eq!(1, 0)
    }
    // #[ignore]
    #[test]
    fn test_phi_g() {
        let catalyst = CatalystData::new(GeometryType::Sphere(0.744/2.0 * 1.0e-2), 0.4886, Some(1.3), 600.0);
        let reaction = ReactionType::TypeI;
        let reaction_data = ReactionData::new(reaction, 4.5747e-5, 2.6700, 3.2769, &[8.5326, 8.5326, 0.0, 0.0], &[0.3172e-4, 0.2203e-4, 0.2009e-4, 0.3684e-4]).unwrap();
        // Note: Ca_eq should be smaller than Ca0
        let phi_g = generalized_thiele_modulus(&reaction_data, &catalyst);
        let eta = eff_factor(phi_g, &catalyst.geometry);
        assert!((eta - 0.9626).abs() < 1e-3)
    }
}