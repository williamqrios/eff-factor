## Effectiveness factor calculator
### Project setup
**Requirements**
- Rust development environment (`rustc`, `cargo`, `rustup`). See [official instructions](https://rust-lang.org/tools/install/)
- Python 3.8+ 
- Any Python package/environment manager (examples use [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions))

**Compilation, testing, and usage**
1. Using conda or pip, create an environment in the project's root directory and install `maturin`
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
import effectiveness_factor as ef
```
4. To test some functionality from within Rust (on Linux)
```
LD_LIBRARY_PATH=./envs/lib cargo test --no-default-features 
```
5. To test some functionality and display the output even if the test passes 
```
LD_LIBRARY_PATH=./envs/lib cargo test --no-default-features -- --nocapture
```

**Using the calculator from a local web server**
1. Install flask 
```
(envs) $ conda install flask 
```  
2. From the root directory start the server (`--debug` enables hot reloading)
```
(envs) $ flask --app ./webdev/app.py run --debug 
```
3. Navigate to the specified address (port `5000` by default)