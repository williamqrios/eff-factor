use ndarray::prelude::*; 
use crate::diffusivities::wc_diffusivity;
use crate::models::{LibError, WilkeChangData, UnifacParams}; 

/// Round to n decimal places  
pub fn round(number: f64, n: i32) -> f64 {
    if number == 0.0 { return number }
    let exponent = f64::log10(f64::abs(number)).floor();
    // Transform number to value between 1 and 10. 
    let whole = number * (10_f64).powi(- exponent as i32); 
    // Apply rounding to desired digits 
    let multiplier = (10.0_f64).powi(n);
    let rounded = (whole * multiplier).round() / multiplier;
    // Recover the exponent part 
    rounded * (10_f64).powi(exponent as i32)
}


pub fn wc_example() -> Array<f64, Ix2> {
    let temp = 78.0 + 273.15; 
        let mass = [60.0520, 46.0684, 88.1051, 18.0153];
        let vbp =  [66.0, 60.85, 106.3, 18.8];
        let c1 = [-9.03, 7.875, 14.354, -52.843];
        let c2 = [1212.3, 781.98, -154.6, 3703.6];
        let c3 = [-0.322, -3.0418, -3.7887, 5.866]; 
        let c4 = [0., 0., 0., -5.879e-29]; 
        let c5 = [0., 0., 0., 10.]; 
        let phi = [1.0, 1.5, 1., 2.6];
        let mut mu = [0.0; 4]; 
        for i in 0..4 {
            mu[i] = f64::exp(c1[i] + c2[i]/temp + c3[i]*f64::ln(temp) + c4[i] * temp.powf(c5[i])) * 1.0e3; 
        }
        let wc = WilkeChangData::new(&vbp, &mass, &phi, &mu, None).unwrap();
        let d = wc_diffusivity(&wc, temp);
        d
}

pub fn unifac_example() -> UnifacParams {
    let R = array![0.9011, 0.6744, 1., 0.92, 1.9031, 1.3013]; 
    let Q = array![0.848, 0.540, 1.2, 1.4, 1.728, 1.224]; 
    let a = array![
        [0.0, 0., 986.5, 1318., 232.1, 663.5],
        [0., 0., 986.5, 1318., 232.1, 663.5], 
        [156.4, 156.4, 0., 353.5, 101.1, 199.],
        [300., 300., -229.1, 0., 72.87, -14.09],
        [114.8, 114.8, 245.4, 200.8, 0., 660.2],
        [315.3, 315.3, -151., -66.17, -256.3, 0.]
    ]; 
    let nu = array![[1., 0., 0., 0., 0., 1.], [1., 1., 1., 0., 0., 0.], [1., 1., 0., 0., 1., 0.], [0., 0., 0., 1., 0., 0.]];
    let temp = 78.0 + 273.15; 
    let params = UnifacParams::new(R, Q, a, nu, temp); 
    params.unwrap()
}

/// General function for inverting matrices.
pub fn invert(matrix: Array2<f64>) -> Result<Array2<f64>, LibError> {
    let n = matrix.shape()[0];
    if n != matrix.shape()[1] {
        return Err(LibError::NotSquare); // Not a square matrix
    }
    
    let result = match n {
        1 => matrix.mapv(|x| 1.0 / x).into_shape_clone((1, 1)).ok(),
        2 => invert_2x2(matrix),
        3 => invert_3x3(matrix),
        _ => None, // Only 1x1, 2x2 and 3x3 matrices are supported
    };

    match result {
        Some(inverse) => Ok(inverse),
        None => Err(LibError::NotInvertible), // Matrix is not invertible (null determinant)
    }
}


fn invert_3x3(matrix: Array2<f64>) -> Option<Array2<f64>> {
    // Compute the determinant
    let det = matrix[[0, 0]] * (matrix[[1, 1]] * matrix[[2, 2]] - matrix[[1, 2]] * matrix[[2, 1]])
            - matrix[[0, 1]] * (matrix[[1, 0]] * matrix[[2, 2]] - matrix[[1, 2]] * matrix[[2, 0]])
            + matrix[[0, 2]] * (matrix[[1, 0]] * matrix[[2, 1]] - matrix[[1, 1]] * matrix[[2, 0]]);
    if det < 1.0e-23 {
        return None; 
    }
    // Compute the cofactor matrix
    let cofactor = Array2::from_shape_vec((3, 3), vec![
        matrix[[1, 1]] * matrix[[2, 2]] - matrix[[1, 2]] * matrix[[2, 1]],
        -(matrix[[1, 0]] * matrix[[2, 2]] - matrix[[1, 2]] * matrix[[2, 0]]),
        matrix[[1, 0]] * matrix[[2, 1]] - matrix[[1, 1]] * matrix[[2, 0]],
        -(matrix[[0, 1]] * matrix[[2, 2]] - matrix[[0, 2]] * matrix[[2, 1]]),
        matrix[[0, 0]] * matrix[[2, 2]] - matrix[[0, 2]] * matrix[[2, 0]],
        -(matrix[[0, 0]] * matrix[[2, 1]] - matrix[[0, 1]] * matrix[[2, 0]]),
        matrix[[0, 1]] * matrix[[1, 2]] - matrix[[0, 2]] * matrix[[1, 1]],
        -(matrix[[0, 0]] * matrix[[1, 2]] - matrix[[0, 2]] * matrix[[1, 0]]),
        matrix[[0, 0]] * matrix[[1, 1]] - matrix[[0, 1]] * matrix[[1, 0]],
    ]).unwrap();

    // Transpose the cofactor matrix to get the adjugate matrix
    let adjugate = cofactor.t();

    // Divide each element of the adjugate matrix by the determinant
    let inverse = adjugate.mapv(|x| x / det);

    Some(inverse)
}

fn invert_2x2(matrix: Array2<f64>) -> Option<Array2<f64>> {
    // Compute the determinant
    let det = matrix[[0, 0]] * matrix[[1, 1]] - matrix[[0, 1]] * matrix[[1, 0]];
    if det.abs() < 1.0e-23 {
        return None;
    }

    // Compute the inverse matrix
    let inverse = Array2::from_shape_vec((2, 2), vec![
        matrix[[1, 1]] / det,
        -matrix[[0, 1]] / det,
        -matrix[[1, 0]] / det,
        matrix[[0, 0]] / det,
    ]).unwrap();

    Some(inverse)
}

/// I_0 function
pub fn bessel_i0(z: f64) -> f64 {
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

pub fn bessel_i1(z: f64) -> f64 {
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