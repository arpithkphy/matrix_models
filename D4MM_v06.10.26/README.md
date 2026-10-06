# D4MM — bosonic D=4 matrix model by Langevin dynamics

A small C++ program that generates configurations of the **bosonic D=4 matrix quantum mechanics** at finite temperature and measures
three observables on every saved configuration:

| observable | definition | meaning |
|---|---|---|
| `P` (Polyakov loop) | `\|tr D\| / N`, `D = diag(exp(i θ_j))` | order parameter of the confinement/deconfinement transition |
| `R2` (extent) | `λ/(N L) Σ_{t,i} tr (X^i_t)²` | size of the matrix configuration, summed over the 3 matrices |
| `E` (energy) | `−3λ/(4 N L) Σ_{t,i,j} tr [X^i_t, X^j_t]²` | energy per N² (virial estimator, exact at λ = 1) |

Only **configuration generation and measurement** are included; there is no analysis code. The output is a plain time series per run
(see "Output"), to be analysed with your own tools.

## The model

Three Hermitian traceless N×N matrices `X^i_t` (i = 1,2,3) on a periodic time lattice of L sites with spacing `a = β/L`, `β = 1/T`,
and the holonomy `D = diag(e^{iθ_1},…,e^{iθ_N})` of the gauge field on the last link (static gauge, `Σθ_j = 0`):

```
S = (N/λ) Σ_t a tr[ (1/(2a²)) (U_t X^i_{t+1} U_t† − X^i_t)²  −  (1/4) Σ_{i,j} [X^i_t, X^j_t]² ]  −  Σ_{j≠k} ln |sin((θ_j − θ_k)/2)|
```

`U_t = 1` except `U_{L−1} = D`. The last term is the Faddeev–Popov (Vandermonde) term of the angles. Units: λ = 1, lengths in `λ^{-1/3}`.
Sampling is by the Langevin equation (Euler–Maruyama step `ε`):

```
X     ← X     − ε ∂S/∂X     + √ε η ,   ⟨η_ij η_kl⟩ = 2 δ_il δ_jk     (then the trace is removed)
θ     ← θ     − ε ∂S/∂θ     + √ε ξ ,   ⟨ξ ξ⟩ = 2                      (then the angles are shifted to sum to zero)
```

All matrices and angles stay real/Hermitian: the action is real, so no complexification is needed.

## Requirements

A C++17 compiler (`g++` or `clang++`) and Python 3 (standard library only, for `run.py`). Tested with Apple clang on macOS.

## Quick start

```bash
python3 run.py list                                  # the parameter sets
python3 run.py run seeds --tier smoke                # a few seconds: checks that everything works
python3 run.py run scan_L20 --tier standard --jobs 4 # 17 temperatures, N=4, L=20  (about 8 minutes on 4 cores)
python3 run.py status scan_L20 --tier standard       # which runs are finished
```

`run.py` copies `src/` into a folder per run, sets N, L and T (they are compile-time constants), compiles, and runs. Results go to
`runs/<set>/<tier>/<name>/`. A finished run is skipped when you call the command again (`--force` repeats it, `--dry-run` only shows
what would run). Your own sets: add an entry to `PARAM_SETS` at the top of `run.py`.

| tier | thermalisation steps | production steps | measure every | N=4, L=20 takes |
|---|---|---|---|---|
| `smoke` | 2 000 | 20 000 | 20 | about 5 s |
| `standard` | 50 000 | 1 600 000 | 100 | about 100 s |
| `long` | 100 000 | 10 000 000 | 100 | about 11 min |

### One run by hand

```bash
# 1. edit the constants N (NCOLOR), L, TEMP in src/model.h
# 2. write params.txt (see below) in an empty working folder
g++ -O2 -std=c++17 src/*.cpp -o d4mm
./d4mm            # reads ./params.txt, writes ./meas/*.txt
```

`params.txt` holds whitespace-separated numbers: `step  thermalization_steps  production_steps  measure_every  [adaptive_factor]  [seed]`.
The last two are optional (defaults 0 = adaptive step off, seed 1). An example is in this folder.

## Output (per run, folder `meas/`)

One line per measured configuration; the first column is the Langevin time since the start of production.

| file | columns |
|---|---|
| `polyakov.txt` | `time  P` |
| `extent.txt` | `time  R2  R2_matrix1  R2_matrix2  R2_matrix3` (the second column is the sum of the last three) |
| `energy.txt` | `time  E` |

`run.log` holds the set-up and a summary of the adaptive step; `run.json` all settings and the run status; `src/` the exact code compiled.

### Using the data

* The first part of every run still remembers the starting configuration: discard it (the `smoke` runs are too short to be physical).
* Successive lines are correlated. Near the transition (T ≈ 1.1) the autocorrelation time is of the order of 100–250 measurements at
  `standard` settings, so a run of 16 000 measurements holds only a few dozen independent samples. Use blocking/jackknife errors.
* The step `ε` is a systematic error (first order in `ε`); the default is 2·10⁻⁴, at most 2 % of the lattice spacing.

## Adaptive step

Before each update the largest force `Kmax` (over all matrices and angles) is compared with 3 × its running average. If it is larger,
that single step is shortened (`ε_n = ε · 3 · average / Kmax`, noise included). This protects against rare huge forces, for example when
two holonomy angles come very close (the Faddeev–Popov force is ∝ 1/distance). It does not change normal steps. Set the fifth number
of `params.txt` (or `ADAPT_FACTOR` in `run.py`) to 0 to switch it off. `run.log` reports how many steps were shortened (typically
0.01–0.1 %, rarely down to ε/1000).

## Check your build

`standard` runs, N=4, L=20, λ=1 (mean ± error from 1.6 M steps, 90 % of each run used):

| T | P | R2 | E |
|---|---|---|---|
| 0.2 | 0.218(5) | 1.204(4) | 1.260(6) |
| 1.0 | 0.485(27) | 1.495(26) | 1.638(36) |
| 2.5 | 0.895(13) | 2.62(7) | 3.59(9) |

A fit of `P(T) = A arctan(B (T − T_c)) + D` to the 17 temperatures of `scan_L20` gives `T_c = 1.078(24)`.

## Source files

```
src/model.h       constants (N, L, T, λ), field variables, matrix class, all declarations
src/matrix.cpp    N×N complex matrix and its operations
src/random.cpp    uniform and Gaussian random numbers, starting matrices, Langevin noise matrix
src/forces.cpp    ∂S/∂X (matrices) and ∂S/∂θ (angles)
src/evolve.cpp    one Langevin update, adaptive step
src/measure.cpp   the observables P, R2, E
src/params.cpp    reads params.txt
src/main.cpp      starting configuration, thermalisation, production loop
run.py            parameter sets and the parallel run script
params.txt        example parameter file
```

## Limits

N, L, T and λ are compile-time constants; the gauge group is SU(N) in the static diagonal gauge; the random numbers come from the C
library `rand()` (the stream differs between platforms); the energy is the virial estimator; the model is bosonic only.

## Created by

Arpith Kumar, Postdoctoral Researcher, Bielefeld University (10.2026)
