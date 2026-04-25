# LaWannTeX

LaWannTeX is a Julia project that reads Wannier90 Hamiltonian files and renders the extracted hopping parameters as TeX tables and equations.

The Korean version of README is available at [README(ko).md](README(ko).md).

## Usage and Required Files

The entry point is `main.jl` in the repository root, and the implementation lives under `src/`.

Run LaWannTeX with a TOML input file:

```bash
julia main.jl --input examples/graphene/input.toml
```

The input TOML must contain a top-level `[files]` section. The usual Wannier90 inputs are:

- `wannier90.win`: used to read the lattice and projection labels.
- `wannier90_hr.dat`: used to read the tight-binding hopping parameters.

Instead of reading structure and orbital metadata from a `.win` file, you can write them directly in the same TOML file with a top-level `[structure]` section and one or more `[[orbitals]]` entries.

## Input TOML Format

- `[files]` requires `hr`, the Wannier90 Hamiltonian file to read, and `out`, the TeX output path.
- `[render]` usually contains `atol` and `include_second_quant`. If you want the largest hoppings in sorted order without selecting specific orbitals, you can also set `max_rblocks`.
- If `[selection]` is omitted, LaWannTeX uses the legacy sorted rendering mode by default.
- In `[selection]`, set `mode = "orbital"` to keep only hopping terms between selected Wannier indices up to the `max_nn` shell.
- `max_nn` is computed as a global non-zero distance class over the selected basis, similar to `neighbors(n)` in [Kwant](https://kwant2.uber.space). `max_nn = 0` keeps only raw home-cell zero-distance onsite terms.
- Structure and orbital metadata can come from `wannier90.win` or from manually written `[structure]` and `[[orbitals]]` sections.
- Shorthand orbital labels are supported in both `[[orbitals]].orbitals` and `wannier90.win` projections. Currently, `s`, `p`, `d`, and `f` expand to Wannier90-order real orbital labels.

## Examples

Example inputs are available under the `examples/` directory.

### 1. Using `wannier90.win`

```toml
[files]
hr = "wannier90_hr.dat"
win = "wannier90.win"
out = "reports/hamiltonian.tex"

[render]
atol = 1.0e-6
max_rblocks = 12
include_second_quant = true

[selection]
mode = "legacy"
```

### 2. Writing Structure and Orbitals Manually

```toml
[structure]
lattice = [
  [2.46, 0.0, 0.0],
  [1.23, 2.130422493, 0.0],
  [0.0, 0.0, 20.0],
]

[[orbitals]]
name = "C"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["pz"]
spins = ["none"]

[[orbitals]]
name = "C"
position_frac = [0.3333333333, 0.3333333333, 0.0]
orbitals = ["pz"]
spins = ["none"]

[files]
hr = "wannier90_hr.dat"
out = "reports/hamiltonian.tex"

[render]
atol = 1.0e-6
max_rblocks = 12
include_second_quant = true

[selection]
mode = "legacy"
```

### 3. Using Orbital Selection Mode

```toml
[files]
hr = "wannier90_hr.dat"
out = "reports/orbital_model.tex"

[render]
atol = 1.0e-6
include_second_quant = true

[selection]
mode = "orbital"
wannier_indices = [1, 2]
max_nn = 1

[structure]
lattice = [
  [2.46, 0.0, 0.0],
  [1.23, 2.130422493, 0.0],
  [0.0, 0.0, 20.0],
]

[[orbitals]]
name = "C"
position_frac = [0.0, 0.0, 0.0]
orbitals = ["pz"]
spins = ["none"]

[[orbitals]]
name = "C"
position_frac = [0.3333333333, 0.3333333333, 0.0]
orbitals = ["pz"]
spins = ["none"]
```
