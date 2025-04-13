use ndarray::prelude::*; 
use crate::diffusivities::{wc_diffusivity, ms_diffusivity};
use crate::utils::{round, wc_example, unifac_example, invert};  
use crate::models::{ReactionType, GeometryType, CatalystData, ReactionData, LibError, UnifacParams};
use crate::eta::{generalized_thiele_modulus, eff_factor};
use crate::unifac::unifac_gamma;

/// Effective diffusivity with ideal model  
pub fn ideal_diffusivity(d_infty: &Array<f64, Ix2>, x: &Array<f64, Ix1>, reaction_type: &ReactionType) -> Array<f64, Ix1> {
    let nu = reaction_type.stoichiometry(); 
    let n = d_infty.shape()[0];
    let d_ms = ms_diffusivity(d_infty, x);
    let mut d_ideal = Array::<f64, Ix1>::zeros(n);
    for i in 0..n {
        for j in 0..n {
            if i != j {
                d_ideal[i] += x[j]/d_ms[[i, j]] * (1.0 - nu[j] * x[i] / (nu[i] * x[j]));
            }
        }
    }
    d_ideal.mapv(|x| 1.0 / x)
}

/// Computes the "B" matrix for the new model.
pub fn b_matrix(d_ms: &Array<f64, Ix2>, x: &Array<f64, Ix1>) -> Array<f64, Ix2> {
    let n = d_ms.shape()[0];
    let d_ms_inv = d_ms.mapv(|d| if d != 0.0 { 1.0 / d } else { 0.0 });
    // Taking care of the diagonal part 
    let sum_bii = x * &d_ms_inv;
    let sum_bii = sum_bii.reversed_axes();
    let sum_bii = sum_bii.sum_axis(Axis(0));
    let mut bii = Array::<f64, Ix2>::zeros((n - 1, n - 1));
    for i in 0..(n-1) {
        bii[[i, i]] = x[i] * d_ms_inv[[i, n - 1]] + sum_bii[i];
    }
    // Off-diagonal elements 
    let mut bij = Array::<f64, Ix2>::zeros((n - 1, n - 1));
    for i in 0..(n-1) {
        for j in 0..(n-1) {
            if i != j {
                bij[[i, j]] = d_ms_inv[[i, j]] - d_ms_inv[[i, n - 1]];
            }
        }
    }
    let x_slice = x.slice(s![..n-1]);
    let b = x_slice.broadcast((n-1, n-1)).unwrap().to_owned().reversed_axes() * bij * (-1.0) + bii; 
    b
}

/// Computes Gamma matrix which depends on activity/fugacity coefficients. Paramater "func" takes in a function that computes the activity/fugacity coefficients for a given composition and system information.
pub fn gamma_matrix(x: &Array<f64, Ix1>, params: &UnifacParams) -> Array<f64, Ix2> {
    let n = x.len();
    let h = 1.0e-8; // step size
    let mut gamma = Array::<f64, Ix2>::zeros((n - 1, n - 1));
    for j in 0..(n-1) {
        let mut x_forward = x.clone(); 
        let mut x_backward = x.clone(); 
        x_forward[j] += h;
        x_backward[j] -= h;
        x_forward[n - 1] -= h; 
        x_backward[n - 1] += h;
        let gamma_forward: Array<f64, Ix1> = unifac_gamma(&x_forward, &params);
        let gamma_backward: Array<f64, Ix1> = unifac_gamma(&x_backward, &params);
        let dln_gamma: Array<f64, Ix1> = gamma_forward.mapv(|x| x.ln()) - gamma_backward.mapv(|x| x.ln());
        
        gamma.column_mut(j).assign(
            &(&x.slice(s![..(n-1)]) * &dln_gamma.slice(s![..(n-1)]) / (2.0 * h))
        );
    }

    gamma + Array::<f64, Ix2>::eye(n - 1)
}


/// Effective diffusivities computed with the new model
pub fn new_diffusivity(d_infty: &Array<f64, Ix2>, x: &Array<f64, Ix1>, reaction_type: &ReactionType, params: &UnifacParams) -> Result<Array<f64, Ix1>, LibError> {
    let j_ratio = reaction_type.flux_ratios(x); 
    let d_ms = ms_diffusivity(d_infty, x);
    let n = d_infty.shape()[0];
    let gamma = gamma_matrix(x, params);
    let b: Array2<f64> = b_matrix(&d_ms, x); 
    let d = invert(b)?.dot(&gamma);
    let j_slice = j_ratio.slice(s![..(n-1), ..(n-1)]);
    let deff_i_reciprocal = (invert(d)? * j_slice).sum_axis(Axis(1));
    let deff_i = deff_i_reciprocal.mapv(|x| 1.0 / x);
    let j_slice_last = j_ratio.slice(s![n-1, ..(n-1)]); 
    let deff_n_reciprocal= (deff_i_reciprocal * j_slice_last * (-1.0)).sum();
    let deff = ndarray::concatenate(Axis(0), &[deff_i.view(), array![1.0/deff_n_reciprocal].view()]).unwrap();
    Ok(deff) 
}


#[cfg(test)]
mod tests {
    use super::*; 
    #[test]
    fn test_j_ratio_type_i() 
    {
        let x = array![0.5, 0.5, 0.0, 0.0]; 
        let reaction = ReactionType::TypeI; 
        let flux = reaction.flux_ratios(&x);
        assert_eq!(flux, array![[1.0, 1.0, -1.0, -1.0], 
            [1.0, 1.0, -1.0, -1.0],
            [-1.0, -1.0, 1.0, 1.0], 
            [-1.0, -1.0, 1.0, 1.0]
            ]);
    }
    #[test]
    fn test_d_ideal() {
        let x = array![0.1898, 0.1898, 0.3102, 0.3102]; 
        let d_infty = wc_example(); 
        let reaction = ReactionType::TypeI;
        let d_ideal = ideal_diffusivity(&d_infty, &x, &reaction);
        let ans = array![0.4596, 0.5465, 0.4826, 0.5172] * 1.0e-4;
        let diff = (d_ideal - ans).abs().sum(); 
        assert!(diff < 1.0e-6);  
    }
    #[test]
    pub fn test_b() {
        let x = array![0.189824390411377, 0.189824390411377, 0.310175609588623, 0.310175609588623];
        let d_infty = wc_example(); 
        let d_ms = ms_diffusivity(&d_infty, &x);
        let b = 1.0e-4 * b_matrix(&d_ms, &x);
        let ans: Array<f64, Ix2> = array![
        [2.299986432352584,  -0.152574879930062,  -0.028521028281142],
        [-0.216033522566799,   2.021727176574570,  -0.023982473091819],
        [-0.258093173298355,  -0.146984836586150,   1.666873808759636],
        ];
        let diff = (b - ans).abs().sum(); 
        assert!(diff < 1.0e-5, "diff = {}", diff);  
    }
    #[test]
    pub fn test_new_model() {
        let x = array![0.189824390411377, 0.189824390411377, 0.310175609588623, 0.310175609588623];
        let d_infty = wc_example(); 
        let reaction = ReactionType::TypeI;
        let d_new = new_diffusivity(&d_infty, &x, &reaction, &unifac_example());
        
        assert!(d_new.is_ok());
        
        let ans = array![
        0.516428975536377,
        0.494102357100941,
        0.390200463315405,
        0.715583181721867
        ] * 1.0e-4; 
        
        let d_new = d_new.unwrap();
        let diff = (&d_new - ans).abs().sum(); 
        assert!(diff < 1.0e-6);
        let catalyst = CatalystData::new(GeometryType::Sphere(0.744/2.0 * 1.0e-1), 0.4886, Some(1.3), 0.600);
        let reaction = ReactionType::TypeI;
        let reaction_data = ReactionData::new(reaction, 4.5747e-5/60.0*1.0e6, 2.6700, 3.2450e-3,  &[3.2482e-3,    3.2482e-3,    5.2844e-3,    5.2844e-3], &d_new.to_vec()).unwrap();
        // Note: Ca_eq should be smaller than Ca0
        let phi_g = generalized_thiele_modulus(&reaction_data, &catalyst);
        let eta = eff_factor(phi_g, &catalyst.geometry);
        assert!(eta > 0.0, "eta = {}", eta);

    }
}