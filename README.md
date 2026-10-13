# Exact Analysis of List-Based Caches with Delayed Hits

MATLAB code accompanying the paper **"Exact Analysis of List-Based Caches with Delayed Hits"** by Roberto Sala, Andrea Marin, Danilo Ardagna, and Giuliano Casale (SIGMETRICS 2027).

## About the paper

*Delayed hits* arise in latency-sensitive caching systems (CDNs, edge computing) when a request finds an item that is not yet in the cache but is already being fetched: the request is withheld until the fetch completes and is then served as a delayed hit, rather than triggering a new miss.

The paper derives the first **exact steady-state solution** of a delayed hit system under the independent reference model (IRM):

- The cache is a **multi-list cache** with `h` lists of sizes `m_1, ..., m_h`, using **RR** (random replacement) or **FIFO** replacement, optionally controlled by per-item **access probabilities** `a_{i,j}` (Markovian policies).
- The retrieval of a missed item is modeled as a **visit of the fetch job to a queueing network** of infinite-server (IS) and processor-sharing (PS) stations with phase-type sojourn times, which captures retrieval workflows, transfer-time distributions, and bandwidth sharing.
- The joint CTMC of the cache and the retrieval system admits a **product-form** stationary distribution. The proof uses Whittle's partial balance method. The result is insensitive to the phase-type distributions beyond their means.
- A **recurrence relation** for the normalizing constant, and an exact reformulation in terms of normalizing constants of caches *without* delayed hits, allow exact metrics for small systems.
- To cope with state-space explosion, a **fixed point iteration (FPI)** approximates hit, miss, and delayed-hit ratios and the mean latency for large systems (thousands of items).
- The FPI is used to **optimize** Markovian replacement policies for hit rate, latency, and station load balancing.

## Model at a glance

| Paper | Code | Meaning |
|---|---|---|
| `n` | `n` | Number of items |
| `m = (m_1, ..., m_h)` | `m` | Sizes of the `h` cache lists (`h = numel(m)`) |
| RR / FIFO | `replacement_policy` | `'rr'` or `'fifo'` |
| `r + 1` stations | `S` | Number of retrieval stations |
| IS / PS | `station_type` | `"INF"` (IS) or `"PS"` per station |
| phases `f_{i,s}` | `F` | Number of PH phases per station |
| `lambda_i` | `p_access` | Item popularity |
| `a_{i,j}` | `access_probabilities` (`n x h`) | Access probabilities (decision variables in the optimization examples) |
| | `access_fetching` (`n x 1`) | Probability that a miss starts a fetch (all ones in the experiments) |
| PH parameters | `alpha`, `T` | PH initial vector and sub-generator for each station and item |
| `R_i` | `R` (`(S+1) x (S+1) x n`) | Routing matrix of the fetch job of each item (first row: entry from the miss state, last column: exit) |
| `gamma_{i,j}` | `gamma` | Cache access factors |
| `eta_{s,i}` | `eta` | Fetching demands |
| `phi_{s,i}`, `xi_j`, `pi_{i,j}` | `phi`, `xi`, `pi_vals` | Outputs of the FPI (`fetching_iteration_mixed`) |

## Repository contents

The three scripts below are the entry points. Each one starts a parallel pool if none is running and calls `addpath(genpath(pwd))`.

### `test_fetching_mixed.m` — validation of the analytical results
Small system (`n = 5`, `m = [1,2]`, FIFO, one IS and one PS station, random PH parameters from `generateData`). It checks, against the exact CTMC solution (`ctmc_solve`) and the simulator, that:

1. the **product-form** stationary distribution (Theorem 4.1) matches the CTMC solution;
2. the normalizing constant computed through the exact recurrence (Sec. 4.2, `recursive_mixed`, `recursive_bis_mixed`) and through the reformulation with caches without delayed hits (Sec. 4.3, `compute_norm_const_mixed`) matches the CTMC;
3. hit, miss, and delayed-hit ratios match Eqs. (17)–(20);
4. the identity in Eq. (23) holds for PS stations;
5. the exact latency (Eqs. (8)–(10)) can be compared with the FPI approximation (Sec. 4.5) and with the simulator (`simulate_trace_mixed_queue_fetching_system`).

A final, self-contained section (it starts with `clear all`) runs the **trace-driven simulation** on the Wikipedia 2019 trace (Sec. 5.3). It requires the trace chunks and the `popularities_*.mat` / `sizes_*.mat` files under `simulator/*`.

### `minimization_problem.m` — maximizing hit rate is not minimizing latency (Sec. 6.1)
`n = 5`, `m = [3]` (single list), RR, one IS and one PS station, Zipf exponent 1.4. Uses `fmincon` (and `MultiStart` for Scenario 3) to optimize `a_{i,j}`:

- **Scenario 1 (identical fetch times):** both objectives share the same optimum. Samples `1e4` random access-probability vectors and plots latency against hit rate. Output: `hit_latency_plot_fast_less_popular.png`.
- **Scenario 2 (slow fetches for less popular items):** the two optimal points lie at opposite extremes of the (hit rate, latency) space. To reproduce, set the fetch rates of the least popular items to `0.01` in `T`.
- **Scenario 3 (sensitivity to the fetch time of item 1):** sweeps the fetch rate of item 1 over `10^-2 ... 10^2` and tracks the optimal `a_{1,1}` for each objective. Output: `optimal_access_probabilities_item1.png`.

> The default `T` in the script (fetch rates `1.0` for items 1–3 and `0.01` for items 4–5) corresponds to a configuration with slower fetches for the less popular items. Adjust `T` to reproduce Scenario 1.

### `load_balancing.m` — load balancing of retrieval stations with the FPI (Sec. 6.2)
`n = 100`, `m = [10, 10]`, FIFO, two PS stations, Zipf exponent 1.0. The first 20 items are fetched much faster at station 1 than at station 2 and are initially routed 75–25% to the two stations, the remaining items 25–75%. The goal is to minimize the maximum station occupancy, in three stages using `MultiStart` + `fmincon`:

| Stage | Decision variables | Objective |
|---|---|---|
| Initial | none (`a_{i,j} = 1`) | baseline occupancies and latency |
| Stage 1 | access probabilities `a_{i,j}` | minimize max occupancy (`compute_max_utilization`) |
| Stage 2 | `a_{i,j}` and routing split `psi` | minimize max occupancy (`compute_max_utilization_with_R`) |
| Stage 3 | `a_{i,j}` and routing split `psi` | minimize latency (`compute_latency_FPI`) |

It prints the latency and station occupancies of each stage (compare with Table 3 of the paper) and saves `optimal_access_probabilities_stage{1,2,3}.png` and `optimal_R_stage{1,2,3}.png` (Figures 6 and 7).

## Helper functions

The scripts rely on the following functions, which must be on the MATLAB path:

| Area | Functions |
|---|---|
| Data generation | `generateData`, `compute_alpha_T_simulator` |
| CTMC (exact) | `generate_states`, `generate_QMatrix_mixed`, `ctmc_solve`, `split_probability_distribution` |
| Product form | `compute_constants`, `compute_gamma`, `compute_eta`, `compute_probability_distribution_mixed`, `compute_norm_const_mixed`, `recursive_mixed`, `recursive_bis_mixed`, `oner` |
| FPI approximation | `fetching_iteration_mixed`, `estimate_mean_number_awaiting_item` |
| Metrics and objectives | `compute_di`, `compute_hit_fast`, `compute_latency`, `compute_latency_FPI`, `compute_max_utilization`, `compute_max_utilization_with_R` |
| Simulation | `simulate_trace_mixed_queue_fetching_system`, `simulate_trace_mixed_real_trace` |

## Requirements

- MATLAB R2021a or newer (the scripts use the `parpool(IdleTimeout=...)` syntax)
- Parallel Computing Toolbox
- Optimization Toolbox (`fmincon`)
- Global Optimization Toolbox (`MultiStart`, `RandomStartPointSet`)

## Getting started

```matlab
% from the repository root
test_fetching_mixed      % validate the product-form results and the FPI on a small system
minimization_problem     % Section 6.1: hit rate vs. latency
load_balancing           % Section 6.2: load balancing with the FPI
```

Notes:

- The exact CTMC solution is only feasible for very small systems (about 10 items): the number of states grows exponentially with the number of items, stations, and cache size.
- `load_balancing.m` runs `MultiStart` over `n*h (+ n)` variables and can be slow. Reduce `nStartPoints` for a quick run.
- `test_fetching_mixed.m` uses random parameters, so results vary slightly between runs.
- Figures are saved as PNG files in the working directory.

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