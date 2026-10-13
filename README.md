# Exact Analysis of List-Based Caches with Delayed Hits

MATLAB code accompanying the paper **"Exact Analysis of List-Based Caches with Delayed Hits"** by Roberto Sala, Andrea Marin, Danilo Ardagna, and Giuliano Casale (SIGMETRICS 2027).

The code models and optimizes a **multi-list cache backed by a fetching system** of queueing stations. The model captures *delayed hits* (requests that arrive for an item already being fetched) and provides:

- an **exact** CTMC solution for small systems,
- an **exact product-form** stationary distribution and normalization constants,
- a fast **fixed-point iteration (FPI)** approximation that scales to large catalogs,
- **optimization** of the cache access probabilities (and routing probabilities) for hit rate, latency, or server utilization,
- a **discrete-event simulator**, including trace-driven simulation on real traces.

## Model at a glance

| Symbol (in code) | Meaning |
|---|---|
| `n` | Number of items in the catalog |
| `m = [m_1, ..., m_h]` | Sizes of the `h` cache lists |
| `replacement_policy` | `'fifo'` or `'rr'` (random replacement) |
| `S` | Number of fetching stations |
| `F` | Number of phases per station (`F(s)` for station `s`) |
| `station_type` | Per-station discipline: `"INF"` (infinite server) or `"PS"` (processor sharing) |
| `p_access` | Item popularity (Zipf with exponent `zipf_rate`) |
| `access_probabilities` (`n x h`) | Probability that a hit/insertion at list `j` is accepted for item `i` (the decision variables `a_{i,j}`) |
| `access_fetching` (`n x 1`) | Probability that a miss triggers a fetch |
| `alpha`, `T` | Phase-type fetch-time parameters (initial vector and sub-generator) for each station and item |
| `R` (`(S+1) x (S+1) x n`) | Routing matrix among fetching stations for each item (row/column `S+1` is the exit) |
| `gamma`, `eta` | Constants of the product-form solution |

## Repository contents

The three scripts in this repo are entry points. They rely on helper functions that are called from the scripts (see [Helper functions](#helper-functions)) and are expected to be on the MATLAB path (each script runs `addpath(genpath(pwd))`).

### `test_fetching_mixed.m` — validation of the model
Validates the analytical results against the exact CTMC and the simulator on a small system (`n = 5`, `m = [1,2]`, FIFO, one `INF` and one `PS` station, random phase-type parameters from `generateData`). It:

1. Builds the CTMC state space and generator `Q`, and solves it exactly (`ctmc_solve`).
2. Compares the exact solution with the **product-form** distribution.
3. Checks the normalization constant computed through the closed form (Sec. 4.3) and the two recursive schemes (Sec. 4.2, `recursive_mixed` / `recursive_bis_mixed`).
4. Checks miss, delayed-hit and hit ratios against Eqs. (17)–(20), and the identity in Eq. (23) for `PS` stations.
5. Runs the simulator (`1e6` events) and the **FPI** approximation.
6. Computes the exact latency (Eqs. (8)–(10)) and compares it with the FPI approximation.
7. Contains a final, independent section that runs the **real-trace simulation** (`n = 1000`, 3 stations, random replacement). This section starts with `clear all` and needs the trace data under `simulator/reduced_chunks_1000_items/`.

Errors are printed to the console (`Computed error on ...`, `Max computed error on ...`).

### `minimization_problem.m` — hit rate vs. latency trade-off (Section 6.1 example)
Small example (`n = 5`, `m = [3]`, random replacement, one `INF` + one `PS` station, Zipf exponent 1.4) showing that **maximizing hit rate and minimizing latency are different objectives**. It:

1. Optimizes `access_probabilities` with `fmincon` for (a) maximum hit rate (`compute_hit_fast`) and (b) minimum latency (`compute_latency`).
2. **Scenarios 1 & 2:** samples `1e4` random access-probability vectors and plots latency vs. hit rate, highlighting the max-hit-rate, min-latency, and initial points. Output: `hit_latency_plot_fast_less_popular.png`.
3. **Scenario 3:** sweeps the service rate of item 1 over `10^-2 ... 10^2` and uses `MultiStart` to track the optimal `a_{1,1}` for each objective. Output: `optimal_access_probabilities_item1.png`.

### `load_balancing.m` — large-scale optimization with the FPI approximation
Larger system (`n = 100`, `m = [10,10]`, FIFO, two `PS` stations, Zipf exponent 1.0) where popular items (`i <= sum(m)`) and unpopular ones have different fetch rates and routing. Uses the FPI to evaluate latency and optimizes in three stages with `MultiStart` + `fmincon`:

| Stage | Decision variables | Objective |
|---|---|---|
| Initial | none (`a = 1`) | baseline latency and occupancies |
| Stage 1 | `access_probabilities` | minimize maximum station utilization (`compute_max_utilization`) |
| Stage 2 | `access_probabilities` and routing `R(1,1,i)` | minimize maximum utilization (`compute_max_utilization_with_R`) |
| Stage 3 | `access_probabilities` and routing `R(1,1,i)` | minimize latency (`compute_latency_FPI`) |

For each stage it reports the latency (`latency`, `latency_util`, `latency_util_r`, `latency_best_latency`) and cache occupancies, and saves the plots:

- `optimal_access_probabilities_stage{1,2,3}.png`
- `optimal_R_stage{1,2,3}.png`

## Helper functions

The scripts call the following functions, which must be on the MATLAB path:

| Area | Functions |
|---|---|
| Data generation | `generateData`, `compute_alpha_T_simulator` |
| CTMC (exact) | `generate_states`, `generate_QMatrix_mixed`, `ctmc_solve`, `split_probability_distribution` |
| Product form | `compute_constants`, `compute_gamma`, `compute_eta`, `compute_probability_distribution_mixed`, `compute_norm_const_mixed`, `recursive_mixed`, `recursive_bis_mixed`, `oner` |
| FPI approximation | `fetching_iteration_mixed`, `estimate_mean_number_awaiting_item` |
| Metrics / objectives | `compute_di`, `compute_hit_fast`, `compute_latency`, `compute_latency_FPI`, `compute_max_utilization`, `compute_max_utilization_with_R` |
| Simulation | `simulate_trace_mixed_queue_fetching_system`, `simulate_trace_mixed_real_trace` |

## Requirements

- MATLAB (R2021a or newer; the scripts use the `parpool(IdleTimeout=...)` name-value syntax)
- Parallel Computing Toolbox (`parpool`, `parfor`, `UseParallel`)
- Optimization Toolbox (`fmincon`)
- Global Optimization Toolbox (`MultiStart`, `RandomStartPointSet`)

## Getting started

```matlab
% from the repository root
test_fetching_mixed      % validate the model on a small system
minimization_problem     % hit rate vs. latency example
load_balancing           % 3-stage optimization on n = 100 items
```

Each script starts a parallel pool if none is running. Figures are saved as PNGs in the working directory.

Notes:

- `load_balancing.m` runs `MultiStart` with `fmincon` over `n*h (+ n)` variables and can take a long time; reduce `nStartPoints` or `n` for a quick run.
- `test_fetching_mixed.m` uses random parameters and access probabilities, so small numerical differences between runs are expected.
- The real-trace section needs the trace chunks and the `popularities_*.mat` / `sizes_*.mat` files under `simulator/reduced_chunks_1000_items/`.

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{sala2027exact,
  title     = {Exact Analysis of List-Based Caches with Delayed Hits},
  author    = {Sala, Roberto and Marin, Andrea and Ardagna, Danilo and Casale, Giuliano},
  booktitle = {SIGMETRICS},
  year      = {2027}
}
```
