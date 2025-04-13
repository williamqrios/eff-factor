use ndarray::prelude::*; 
use crate::utils::{round, wc_example}; 
use crate::models::WilkeChangData;

/// Infinite-dilution diffusion coefficient for all pairs of molecules. 
pub fn wc_diffusivity(wc: &WilkeChangData, temp: f64) -> Array<f64, Ix2> 
{
    let length = wc.phi.len();
    let mut d = Array::<f64, _>::zeros((length, length).f());
    for i in 0..length {
        for j in 0..length {
            if i != j {
                d[[i, j]] = 7.4e-8 * (wc.phi[j] * wc.mw[j]).powf(0.5) * temp / (wc.visc[j] * wc.vbp[i].powf(0.6)); 
            }
        }
    }
    d
}


/// Computes Maxwell-Stefan diffusion coefficients matrix based on the input infinite dilution matrix and the vector of compositions. 
pub fn ms_diffusivity(d_infty: &Array<f64, Ix2>, x: &Array<f64, Ix1>) -> Array<f64, Ix2> {
    let n = d_infty.shape()[0]; 
    let mut d_ms = Array::<f64, Ix2>::zeros((n, n).f());
    for i in 0..n {
        for j in 0..n {
            if i != j {
                d_ms[[i, j]] = d_infty[[i, j]].powf(0.5 * (1.0 + x[j] - x[i])) * d_infty[[j, i]].powf(0.5 * (1.0 + x[i] - x[j])); 
            }
        }
    }
    d_ms
}



#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_wilke_chang_diffusivity() {
        let d = wc_example(); 
        dbg!(&d); 
        let d_rounded = d.mapv(|x| round(x * 1.0e3, 2)); 
        let ans = array![
            [0.0, 0.0397, 0.0789, 0.0395],
            [0.0299, 0.0, 0.0828, 0.0414],
            [0.0214, 0.0298, 0.0, 0.0296],
            [0.0605, 0.0843, 0.1676, 0.0]
        ];
        let diff = (d_rounded - ans).abs().sum();
        assert!(diff < 1.0e-3, "diff = {diff}");
    }

    #[test]
    fn test_maxwell_stefan() {
        let x = array![0.1898, 0.1898, 0.3102, 0.3102]; 
        let d_infty = wc_example(); 
        let d_ms = ms_diffusivity(&d_infty, &x).mapv(|x| round(x * 1.0e4, 3));
        let ans = array![
        [0.,  0.3443,    0.4442,    0.4760],
        [0.3443,         0.,    0.5283,    0.5661],
        [0.4442,    0.5283,         0.,    0.7047],
        [0.4760,    0.5661,    0.7047,         0.]];
        let diff = (d_ms - ans).abs().sum();
        assert!(diff < 1.0e-10);
    }
}