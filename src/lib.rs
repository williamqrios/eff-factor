#![allow(non_snake_case)]
#![allow(unused)]
mod eta; mod diffusivities; mod effective_diffusivity; mod utils; mod models; mod unifac;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use ndarray::prelude::*;
use crate::models::*;
use crate::eta::{generalized_thiele_modulus, eff_factor};
use crate::diffusivities::wc_diffusivity;
use crate::effective_diffusivity::{ideal_diffusivity, new_diffusivity};

#[pyfunction]
/// Computes the effectiveness factor, user provides required information as well as the pre-computed effective (mixture) diffusivities. 
fn ef_precomputed_diffusivities(reaction_type: &str, catalyst_geometry: &str, catalyst_characteristic_length: f64, catalyst_porosity: f64, catalyst_density: f64, catalyst_tortuosity: Option<f64>, kinetic_constant: f64, equilibrium_constant: f64, equilibrium_Ca: f64, concentrations: Vec<f64>, diffusivities: Vec<f64>) -> PyResult<f64> {
    let reaction = ReactionType::from_str(reaction_type).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction type."))?;
    let geometry = GeometryType::from_str(catalyst_geometry, catalyst_characteristic_length);
    let catalyst = CatalystData::new(geometry, catalyst_porosity, catalyst_tortuosity, catalyst_density);
    let reaction_data = ReactionData::new(
        reaction,
        kinetic_constant,
        equilibrium_constant,
        equilibrium_Ca,
        &concentrations,
        &diffusivities
    ).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction data."))?;
    let generalized_phi = generalized_thiele_modulus(&reaction_data, &catalyst);
    let eta = eff_factor(generalized_phi, &catalyst.geometry);
    Ok(eta)
}


#[pyfunction]
/// Computes the effectiveness factor, user provides required information, but not diffusivities, which are estimated with the ideal model (assumes ideal mixture). Diffusivities in the limit of infinite dilution are computed with the Wilke-Chang equation, so the user must provide Wilke-Chang inputs for each of the components (molecular weight, association factor, molar volume at normal boiling point in cm3/mol, viscosities in cP), and the system temperature. 
fn ef_diffusivities_ideal(reaction_type: &str, catalyst_geometry: &str, catalyst_characteristic_length: f64, catalyst_porosity: f64, catalyst_density: f64, catalyst_tortuosity: Option<f64>, kinetic_constant: f64, equilibrium_constant: f64, equilibrium_Ca: f64, concentrations: Vec<f64>, molecular_weights: Vec<f64>, association_factors: Vec<f64>, molar_volumes_bp: Vec<f64>, viscosities: Vec<f64>, temperature: f64) -> PyResult<f64> {
    let reaction = ReactionType::from_str(reaction_type).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction type."))?;
    let geometry = GeometryType::from_str(catalyst_geometry, catalyst_characteristic_length);
    let catalyst = CatalystData::new(geometry, catalyst_porosity, catalyst_tortuosity, catalyst_density);
    let wc = WilkeChangData::new(
        &molar_volumes_bp,
        &molecular_weights,
        &association_factors,
        &viscosities,
        Some(reaction.num_components())
    ).map_err(|_| PyErr::new::<PyValueError, _>("Invalid Wilke-Chang data."))?;
    let d_infty = wc_diffusivity(&wc, temperature);
    let ct = concentrations.iter().sum::<f64>(); // total concentration
    let x = Array::from_vec(concentrations.clone()).mapv(|c| c / ct ); // mole fractions
    let diffusivities = ideal_diffusivity(&d_infty, &x, &reaction).to_vec();
    let reaction_data = ReactionData::new(
        reaction,
        kinetic_constant,
        equilibrium_constant,
        equilibrium_Ca,
        &concentrations,
        &diffusivities
    ).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction data."))?;
    let generalized_phi = generalized_thiele_modulus(&reaction_data, &catalyst);
    let eta = eff_factor(generalized_phi, &catalyst.geometry);
    Ok(eta)
}


#[pyfunction]
/// Computes the effectiveness factor, user provides required information, but not diffusivities, which are estimated with the new model using UNIFAC to compute activity coefficients. Requires the same inputs as the previous functions, as well as UNIFAC-specific parameters (R, Q, a, nu).
fn ef_diffusivities_new(reaction_type: &str, catalyst_geometry: &str, catalyst_characteristic_length: f64, catalyst_porosity: f64, catalyst_density: f64, catalyst_tortuosity: Option<f64>, kinetic_constant: f64, equilibrium_constant: f64, equilibrium_Ca: f64, concentrations: Vec<f64>, molecular_weights: Vec<f64>, association_factors: Vec<f64>, molar_volumes_bp: Vec<f64>, viscosities: Vec<f64>, temperature: f64, unifac_R: Vec<f64>, unifac_Q: Vec<f64>, unifac_a: Vec<Vec<f64>>, unifac_nu: Vec<Vec<f64>>) -> PyResult<f64> {
    let reaction = ReactionType::from_str(reaction_type).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction type."))?;
    let geometry = GeometryType::from_str(catalyst_geometry, catalyst_characteristic_length);
    let catalyst = CatalystData::new(geometry, catalyst_porosity, catalyst_tortuosity, catalyst_density);
    let wc = WilkeChangData::new(
        &molar_volumes_bp,
        &molecular_weights,
        &association_factors,
        &viscosities,
        Some(reaction.num_components())
    ).map_err(|_| PyErr::new::<PyValueError, _>("Invalid Wilke-Chang data."))?;
    let d_infty = wc_diffusivity(&wc, temperature);
    let ct = concentrations.iter().sum::<f64>(); // total concentration
    let x = Array::from_vec(concentrations.clone()).mapv(|c| c / ct ); // mole fractions
    let R = Array::from_vec(unifac_R);
    let Q = Array::from_vec(unifac_Q);
    let a = Array::from_shape_vec((unifac_a.len(), unifac_a[0].len()), unifac_a.concat()).map_err(|_| PyErr::new::<PyValueError, _>("Invalid UNIFAC \"a\" matrix."))?;
    let nu = Array::from_shape_vec((unifac_nu.len(), unifac_nu[0].len()), unifac_nu.concat()).map_err(|_| PyErr::new::<PyValueError, _>("Invalid UNIFAC groups per molecule matrix."))?;
    let unifac_params = UnifacParams::new(R, Q, a, nu, temperature).map_err(|_| PyErr::new::<PyValueError, _>("Invalid UNIFAC parameters."))?;
    let diffusivities = new_diffusivity(&d_infty, &x, &reaction, &unifac_params).map_err(|_| PyErr::new::<PyValueError, _>("Error computing diffusivities."))?.to_vec();
    let reaction_data = ReactionData::new(
        reaction,
        kinetic_constant,
        equilibrium_constant,
        equilibrium_Ca,
        &concentrations,
        &diffusivities
    ).map_err(|_| PyErr::new::<PyValueError, _>("Invalid reaction data."))?;
    let generalized_phi = generalized_thiele_modulus(&reaction_data, &catalyst);
    let eta = eff_factor(generalized_phi, &catalyst.geometry);
    Ok(eta)
}


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn effectiveness_factor(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(ef_precomputed_diffusivities, m)?)?;
    m.add_function(wrap_pyfunction!(ef_diffusivities_ideal, m)?)?;
    m.add_function(wrap_pyfunction!(ef_diffusivities_new, m)?)?;
    Ok(())
}

