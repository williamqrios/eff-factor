/// Thermodynamic models to compute fugacity/activity coefficients. 
use ndarray::prelude::*; 
use crate::utils::unifac_example;
use crate::models::UnifacParams;


pub fn unifac_gamma(x: &Array<f64, Ix1>, params: &UnifacParams) -> Array<f64, Ix1> {
    let x = x.mapv(|val| if val <= 0.0 { 1.0e-20 } else { val }); 
    let r = params.get_r();
    let q = params.get_q();

    let num_groups = params.R.len(); 
    let num_components = x.len();

    // Molecular volume fraction
    let Phi = (&r * &x)/r.dot(&x);
    // Molecular surface fraction
    let theta = (&q * &x)/q.dot(&x);

    // Combinatorial term of the activity coefficient 
    let ln_gamma_comb = &Phi.mapv(f64::ln) - &x.mapv(f64::ln)
    + (1.0 - &Phi / &x)
    - 5.0 * &q * (&Phi.mapv(f64::ln) - &theta.mapv(f64::ln) + (1.0 - &Phi / &theta));

    // Residual term computations below 

    // Group mole fractions
    let X = x.dot(&params.nu)/(x.dot(&params.nu)).sum();
    
    // Surface area fraction 
    let Theta = &X * &params.Q/(&X * &params.Q).sum(); 
    
    // Molecular energy 
    let Psi = params.a.mapv(|val| (-val/params.temp).exp()); 

    // Creating an internal utility function here, which will get reused later on.
    fn get_sum_term(ng: usize, Th: &Array<f64, Ix1>, Ps: &Array<f64, Ix2>) -> Array<f64, Ix1> {
        let mut sum_term = Array::<f64, Ix1>::zeros(ng);
        for i in 0..ng {
            for j in 0..ng {
                sum_term[i] += Th[j] * Ps[[i,j]] / (&Th.dot(&Ps.column(j)));
            }
        }
        sum_term
    }

    let st = get_sum_term(num_groups, &Theta, &Psi);

    let ln_Gamma = &params.Q * ( 1.0 - &Theta.dot(&Psi).mapv(|val| val.ln()) - st );

    let mut ln_Gamma_pure = Array::<f64, Ix2>::zeros((num_components, num_groups));

    for j in 0..num_components {
        let logical_nu = params.nu.row(j).mapv(|val| if val > 0.0 { 1.0 } else { 0.0 });
        // Copmponent made up of only one group
        if logical_nu.sum() == 1.0 {
            ln_Gamma_pure.row_mut(j).assign(
                &Array::<f64, Ix1>::zeros(num_groups)
            );
            continue;
        }
        let new_X = &params.nu.row(j)/params.nu.row(j).sum();
        let new_Theta = &new_X * &params.Q/(&new_X * &params.Q).sum();
        let new_st = get_sum_term(num_groups, &new_Theta, &Psi);
        let new_Q = &params.Q * logical_nu;
        ln_Gamma_pure.row_mut(j).assign(
            &(new_Q * (1.0 - new_Theta.dot(&Psi).mapv(|val| val.ln()) - new_st))
        );
    }

    let ln_gamma_res = (&params.nu * (ln_Gamma - ln_Gamma_pure)).sum_axis(Axis(1));

    let gamma = (ln_gamma_comb + ln_gamma_res).mapv(|val| val.exp());

    gamma
}


#[cfg(test)]
mod tests {
    use std::array;

    use super::*; 
    #[test]
    fn test_activity_coefficients() {
        let x = array![0.189824390411377, 0.189824390411377, 0.310175609588623, 0.310175609588623];
        let params = unifac_example();
        let gamma = unifac_gamma(&x, &params);
        let ans = array![
            0.836300435649362,
            1.169665823772968,
            1.677618672740010,
            2.169669575188481
        ];
        let diff = (gamma - ans).abs().sum();
        assert!(diff < 1.0e-6, "diff = {diff}");
    }
}
