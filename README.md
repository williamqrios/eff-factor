
Running Python from Rust:
1. Using conda or pip, set up an environment in the project's root directory and install `maturin`
```
conda create --prefix ./envs maturin 
```
2. Activate the environment and build the package with `maturin develop`
```
conda activate ./envs 
maturin develop 
```
3. Use `testing.py` to test out functionality from Python, by first importing the package
```
import eff_factor
```
4. To test some functionality from within Rust (on Linux)
```
LD_LIBRARY_PATH=./envs/lib cargo test --no-default-features 
```