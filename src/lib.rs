mod eta; 
use pyo3::prelude::*;

 

/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
#[pymodule]
fn string_sum(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sum_as_string_positive_numbers() {
        let result = sum_as_string(2, 3).unwrap();
        assert_eq!(result, "5");
    }

    #[test]
    fn test_sum_as_string_zero() {
        let result = sum_as_string(0, 0).unwrap();
        assert_eq!(result, "0");
    }

    #[test]
    fn test_sum_as_string_large_numbers() {
        let result = sum_as_string(1_000_000, 2_000_000).unwrap();
        assert_eq!(result, "3000000");
    }

    #[test]
    fn test_sum_as_string_mixed_numbers() {
        let result = sum_as_string(0, 42).unwrap();
        assert_eq!(result, "42");
    }
}
