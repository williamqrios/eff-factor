use std::fmt::Display; 
use ndarray::prelude::*;

#[derive(Debug)]
pub enum LibError {
    NotSquare,
    NotInvertible,
    UnknownReaction,
    UnknownGeometry,
    InconsistentComponents,
    NullConcentration,
    InvalidInput, 
    InvalidDimensions,
    InvalidConcentration,
}

#[derive(Debug)]
pub enum ReactionType {
    TypeI, // A + B <-> C + D
    TypeII, // 2A <-> C + D
    TypeIII, // A + B <-> 2 C 
    TypeIV, // A <-> C + D 
    TypeV, // A + B <-> C 
    TypeVI, // A <-> C
    TypeVII // 2 A + B <-> C + D 
}

impl ReactionType {
    pub fn from_str(name: &str) -> Result<Self, LibError> {
        match name {
            "TypeI" => Ok(ReactionType::TypeI),
            "TypeII" => Ok(ReactionType::TypeII),
            "TypeIII" => Ok(ReactionType::TypeIII),
            "TypeIV" => Ok(ReactionType::TypeIV),
            "TypeV" => Ok(ReactionType::TypeV),
            "TypeVI" => Ok(ReactionType::TypeVI),
            "TypeVII" => Ok(ReactionType::TypeVII),
            _ => Err(LibError::UnknownReaction),
        }
    }
    pub fn num_components(&self) -> usize {
        use ReactionType::*;
        match self {
            TypeI => 4,
            TypeII => 3,
            TypeIII => 3,
            TypeIV => 3,
            TypeV => 3,
            TypeVI => 2,
            TypeVII => 4,
        }
    }
    pub fn stoichiometry(&self) -> Array<f64, Ix1> {
        use ReactionType::*;
        match self {
            TypeI => array![-1.0, -1.0, 1.0, 1.0],
            TypeII => array![-2.0, 1.0, 1.0],
            TypeIII => array![-1.0, -1.0, 2.0],
            TypeIV => array![-1.0, 1.0, 1.0],
            TypeV => array![-1.0, -1.0, 1.0],
            TypeVI => array![-1.0, 1.0],
            TypeVII => array![-2.0, -1.0, 1.0, 1.0],
        }
    }
    /// Diffusion flux ratios (J_j / J_i) for each reaction type. 
    pub fn flux_ratios(&self, x: &Array<f64, Ix1>) -> Array<f64, Ix2> {
        let n = x.len(); 
        let mut flux = Array::<f64, Ix2>::zeros((n, n));
        let nu = self.stoichiometry();
        if nu.sum() == 0.0 {
            // If sum(stoichiometric coefficients) = 0, J_j / J_i = nu_j / nu_i. 
            for i in 0..n {
                for j in 0..n 
                {
                    flux[[i, j]] = if nu[i] != 0.0 && nu[j] != 0.0 
                        { nu[j] / nu[i] } else { 0.0 };
                }
            }
        } else {
            // Otherwise, J_j / J_i = (N_j - x_j * N_total) / (N_i - x_i * N_total)
            let nt = nu.sum(); 
            for i in 0..n {
                for j in 0..n {
                        flux[[i, j]] = if nu[i] != 0.0 && nu[j] != 0.0 
                            { 
                                (nu[j] - x[j] * nt)  / (nu[i] - x[i] * nt) 
                            } else { 0.0 };
                }
            }
        }
        flux
    }
}


#[derive(Debug)]
/// Inner value -> characteristic length
/// The characteristic length of a slab is the semi-thickness 
/// The characteristic length of a sphere or cylinder is the radius 
/// The characteristic length of any other geometry type is the Volume divided by Surface area. 
pub enum GeometryType {
    Slab(f64), 
    Sphere(f64),
    Cylinder(f64), 
    Other(f64)
}
impl GeometryType {
    pub fn from_str(name: &str, cl: f64) -> Self {
        match name {
            "Slab" => GeometryType::Slab(cl),
            "Sphere" => GeometryType::Sphere(cl),
            "Cylinder" => GeometryType::Cylinder(cl),
            _ => GeometryType::Other(cl)
        }
    }
    pub fn get_value(&self) -> f64 {
        match self {
            GeometryType::Slab(cl) => *cl,
            GeometryType::Sphere(cl) => *cl,
            GeometryType::Cylinder(cl) => *cl,
            GeometryType::Other(cl) => *cl
        }
    }
}

#[derive(Debug)]
/// Aggregates catalyst information  
/// Epsilon: porosity (between 0 and 1) 
/// Tau: tortuosity (between 0 and 1), 
/// Rho_p: catalyst density. Units do not matter as long as they are consistent with the rest of the data. 
pub struct CatalystData {
    pub geometry: GeometryType,
    pub epsilon: f64, // Porosity 
    pub tau: f64, // Tortuosity 
    pub rho_p: f64, // Density 
}


impl CatalystData {
    pub fn new(geometry: GeometryType, epsilon: f64, tau: Option<f64>, rho_p: f64) -> Self {
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


#[derive(Debug)]
/// Aggregates reaction information (kinetics, equilibrium, concentrations, ...)
pub struct ReactionData {
    pub reaction_type: ReactionType, 
    pub k: f64, // Reaction rate constant
    pub Kc: f64, // Equilibrium constant
    pub Ca_eq: f64, // Equilibrium concentration for limiting reactant 
    pub conc: Vec<f64>,
    pub deff: Vec<f64>, // Effective diffusivities in multicomponent mixture
}

impl ReactionData {
    pub fn new(reaction_type: ReactionType, k: f64, Kc: f64, Ca_eq: f64, conc: &[f64], deff: &[f64]) -> Result<Self, LibError> {
        // Consistency checks
        // Length of vector vs. the type of reaction chosen 
        if conc.len() != reaction_type.num_components() || deff.len() != reaction_type.num_components() {
            return Err(LibError::InconsistentComponents);
        }
        // Negative values 
        if conc.iter().any(|&c| c < 0.0) || deff.iter().any(|&d| d <= 0.0) || Kc < 0.0 || k < 0.0 || Ca_eq < 0.0 {
            return Err(LibError::InvalidInput);
        }
        // Null values for concentration of the limiting reactant ("A")
        if conc[0] == 0.0 || Ca_eq == 0.0 {
            return Err(LibError::NullConcentration);
        }
        // Equilibrium concentration is larger than the current concentration for the limiting reactant 
        if conc[0] < Ca_eq {
            return Err(LibError::InvalidConcentration);
        }

        Ok(Self { 
            reaction_type, 
            k, 
            Kc, 
            Ca_eq, 
            conc: conc.to_vec(),
            deff: deff.to_vec(),
         })

    }
}

pub struct WilkeChangData {
    pub vbp: Vec<f64>, // Molar volume at normal boiling point (cm^3/mol) for pure components 
    pub mw: Vec<f64>, // Molecular weight (g/mol)
    pub phi: Vec<f64>, // Association factors
    pub visc: Vec<f64>, // Viscosities (cP)
}

impl WilkeChangData {
    pub fn new(vbp: &[f64], mw: &[f64], phi: &[f64], visc: &[f64], num_components: Option<usize>) -> Result<Self, LibError> {
        // Num components is added to function parameters for consistency checks. 
        if num_components.is_some() {
            if [vbp, mw, phi, visc].iter().any(|x| x.len() != num_components.unwrap()) {
                return Err(LibError::InconsistentComponents);
            }
        }
        if [vbp, mw, phi, visc].iter().any(|x| x.iter().any(|&c| c <= 0.0)) {
            return Err(LibError::InvalidInput);
        }
        Ok(Self { vbp: vbp.to_vec(), mw: mw.to_vec(), phi: phi.to_vec(), visc: visc.to_vec() })
    }
}


pub struct UnifacParams {
    pub R: Array<f64, Ix1>,
    pub Q: Array<f64, Ix1>,
    pub a: Array<f64, Ix2>,
    pub nu: Array<f64, Ix2>, // Occurrence of each group per molecule 
    pub temp: f64 // System temperature 
}

impl UnifacParams {
    pub fn new(R: Array<f64, Ix1>, Q: Array<f64, Ix1>, a: Array<f64, Ix2>, nu: Array<f64, Ix2>, temp: f64) -> Result<Self, LibError> {
        // Consistency checks: 
        // R must have length equal to Q (number of groups)
        // a must be a square matrix with dimensions num_groups x num_groups 
        // nu must be a matrix with dimensions num_components x num_groups
        let ng = R.len(); 
        if ng != Q.len() || a.shape()[0] != a.shape()[1] || a.shape()[0] != ng || nu.shape()[1] != ng {
            return Err(LibError::InvalidDimensions);
        }

        Ok(UnifacParams { R, Q, a, nu, temp })
    }
    pub fn get_r(&self) -> Array<f64, Ix1> {
        self.R.dot(&self.nu.t())
    }
    pub fn get_q(&self) -> Array<f64, Ix1> {
        self.Q.dot(&self.nu.t())
    }
}

#[cfg(test)]
mod tests {
    use super::*; 
    #[test]
    fn test_j_ratios_rx_vii() {
        let reaction = ReactionType::TypeVII; 
        let x = array![0.6698, 0.3302, 0.0, 0.0];
        let j_ratios = reaction.flux_ratios(&x); 
        let ans = array![
            [1.000, 0.5036, -0.7518, -0.7518],
            [1.9858, 1.000, -1.4929, -1.4929],
            [-1.3302, -0.6698, 1.000, 1.000],
            [-1.3302, -0.6698, 1.000, 1.000]
        ];
        let diff = (j_ratios - ans).abs().sum(); 
        assert!(diff < 1e-3, "diff = {diff}");

    }
}