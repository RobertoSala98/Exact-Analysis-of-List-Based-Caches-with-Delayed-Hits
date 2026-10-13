%% Setup the environment
clear all
clc

if isempty(gcp('nocreate'))
    pool = parpool(IdleTimeout=1440);
end

addpath(genpath(pwd));

%% System Parameters
n = 5;
m = [1,2]; h = numel(m);
replacement_policy = 'fifo';
S = 2;
F = [3,3]; F = F(1:S);
station_type = ["INF","PS"]; station_type = station_type(1:S);
sym_computation = false;

zipf_rate = 1.0;
ranks = 1:n;
H = sum(1 ./ (ranks.^zipf_rate));
p_access = (1 ./ (ranks.^zipf_rate)) / H;

[alpha, T, R] = generateData(n, S, F, 16, 0.65);
for ii = 1:n
    for ss = 1:S
        alpha{ss}(:,:,ii) = alpha{ss}(:,:,ii) / sum(alpha{ss}(:,:,ii));
    end
    R(1,1:S,ii) = R(1,1:S,ii) / sum(R(1,1:S,ii));
    R(1,S+1,ii) = 0;
end

access_probabilities = rand(n, h);
access_fetching = ones(n, 1);

%% CTMC states and infinitesimal generator matrix
States = generate_states(n, m, replacement_policy, S, F);
fprintf("Your system has %d states.\n", size(States, 1))
Q = generate_QMatrix_mixed(n, m, replacement_policy, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, States, sym_computation);

%% Exact stationary distribution solving the CTMC
pi = ctmc_solve(Q);
[miss_ratios, delayed_hits_ratios, hit_ratios] = split_probability_distribution(n, m, S, F, States, sym_computation, pi);

norm_const = 1;
idx_cache = 1;
for cc = 1:sum(m)
    while cc > sum(m(1:idx_cache))
        idx_cache = idx_cache + 1;
    end
    if isequal(replacement_policy, 'fifo')
        norm_const = norm_const * p_access(sum(m) - cc + 1)^idx_cache * access_fetching(sum(m) - cc + 1);

        for g = 1:idx_cache
            norm_const = norm_const * access_probabilities(sum(m) - cc + 1, g);
        end
    else
        norm_const = norm_const * p_access(cc)^idx_cache * access_fetching(cc);

        for g = 1:idx_cache
            norm_const = norm_const * access_probabilities(cc, g);
        end
    end
end
norm_const = norm_const / pi(1);

%% Exact product-form solution
constants = compute_constants(n, S, F, p_access, alpha, T, R, sym_computation, access_fetching);
[pi_est, norm_const_est] = compute_probability_distribution_mixed(n, m, S, p_access, station_type, access_probabilities, access_fetching, States, sym_computation, constants);
fprintf("Computed error on stationary probability distribution: %d.\n", max(abs(pi' - pi_est)))

%% Normalization constant computed as in 4.3
gamma = compute_gamma(n, h, p_access, access_probabilities, access_fetching);
eta = compute_eta(n, S, constants);
exact_norm_const = compute_norm_const_mixed(n, m, S, gamma, eta, station_type, replacement_policy);
fprintf("Computed error on exact normalization constant: %d\n", abs(norm_const - exact_norm_const))

%% Normalization constant computed as in 4.2
norm_const_comp = recursive_mixed(n, m, S, F, p_access, constants, [1:n], station_type, access_probabilities, access_fetching, sym_computation, zeros(1,S), replacement_policy);
fprintf("Computed error on normalization constant: %d\n", abs(norm_const - norm_const_comp))

norm_const_comp_bis = recursive_bis_mixed(n, m, S, F, p_access, constants, [1:n], station_type, access_probabilities, access_fetching, sym_computation, zeros(1,S), zeros(1,n), replacement_policy);
fprintf("Computed error on normalization constant bis: %d\n", abs(norm_const - norm_const_comp_bis))

%% Probabilities computed using (17)--(20), Exact model
max_err_miss = 0;
max_err_delayed = 0;
max_err_hit = 0;

for item_idx = 1:n

    p_access_red = [];
    constants_red = [];
    idxs_red = [];

    for kk = 1:n
        if kk ~= item_idx
            p_access_red = [p_access_red p_access(kk)];
            idxs_red = [idxs_red kk];
        end
    end

    for ss = 1:S
        constants_red{ss} = constants{ss};
        constants_red{ss}(:,:,item_idx) = [];
    end

    access_probabilities_red = access_probabilities;
    access_probabilities_red(item_idx, :) = [];
    access_fetching_red = access_fetching;
    access_fetching_red(item_idx) = [];

    PS_values = zeros(1,S);

    miss = miss_ratios(item_idx);
    computed_miss = recursive_mixed(n-1, m, S, F, p_access_red, constants_red, idxs_red, station_type, access_probabilities_red, access_fetching_red, sym_computation, PS_values, replacement_policy) / norm_const_comp;
    max_err_miss = max(max_err_miss, abs(double(miss - computed_miss)));

    computed_delayed = [];
    for ss = 1:S
        computed_delayed = [computed_delayed, {}];

        PS_values_red = PS_values;
        if isequal(station_type(ss), "PS")
            PS_values_red(ss) = PS_values_red(ss) + 1;
        end

        for ff = 1:F(ss)
            delayed = delayed_hits_ratios{ss}(1,ff,item_idx);
            computed_delayed = constants{ss}(1,ff,item_idx) * recursive_mixed(n-1, m, S, F, p_access_red, constants_red, idxs_red, station_type, access_probabilities_red, access_fetching_red, sym_computation, PS_values_red, replacement_policy) / norm_const_comp;
            max_err_delayed = max(max_err_delayed, abs(double(delayed - computed_delayed)));
        end

    end

    for list_num = 1:length(m)
        hit = hit_ratios(list_num, item_idx);

        computed_hit = 1;
        for g = 1:list_num
            computed_hit = computed_hit * access_probabilities(item_idx, g);
        end
        computed_hit = computed_hit * p_access(item_idx)^list_num * access_fetching(item_idx) * recursive_mixed(n-1, oner(m,list_num), S, F, p_access_red, constants_red, idxs_red, station_type, access_probabilities_red, access_fetching_red, sym_computation, PS_values, replacement_policy) / norm_const_comp;

        if isequal(replacement_policy, 'fifo')
            computed_hit = computed_hit * m(list_num);
        end
        max_err_hit = max(max_err_hit, abs(hit - computed_hit));
    end
end

fprintf("Max computed error on hit ratios: %d\n", double(max_err_hit))
fprintf("Max computed error on delayed hits ratios: %d\n", double(max_err_delayed))
fprintf("Max computed error on miss ratios: %d\n", double(max_err_miss))

%% Simulator
num_events = 1e6;
[alpha_sim, T_sim] = compute_alpha_T_simulator(n, S, alpha, T);
[miss_ratios, delayed_hits_ratios, hit_ratios, mean_req_fetching, mean_fetch_time, mean_squared_fetch_time, mean_number_awaiting_completition, mean_arrived_jobs_during_fetch, mean_latency, utilization] = simulate_trace_mixed_queue_fetching_system(n, m, replacement_policy, S, F, station_type, access_probabilities, access_fetching, p_access, alpha_sim, T_sim, R, num_events);

%% FPI method
[phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta', m, station_type);

%% Numerical verification of (23)
err = 0;

for ii = 1:n
    for ss = 1:S
        if isequal(station_type(ss), "PS")

            p_access_red = [];
            constants_red = [];
            idxs_red = [];

            for kk = 1:n
                if kk ~= ii
                    p_access_red = [p_access_red p_access(kk)];
                    idxs_red = [idxs_red kk];
                end
            end

            for tt = 1:S
                constants_red{tt} = constants{tt};
                constants_red{tt}(:,:,ii) = [];
            end

            access_probabilities_red = access_probabilities;
            access_probabilities_red(ii, :) = [];
            access_fetching_red = access_fetching;
            access_fetching_red(ii) = [];

            PS_values = zeros(1,S);

            del_hit = sum(delayed_hits_ratios{ss}(:,:,ii));

            States_red = generate_states(n-1, m, replacement_policy, S, F);
            [pi_red, norm_const_red] = compute_probability_distribution_mixed(n-1, m, S, p_access_red, station_type, access_probabilities_red, access_fetching_red, States_red, sym_computation, constants_red);
            [miss_ratios_red, delayed_hits_ratios_red, hit_ratios_red] = split_probability_distribution(n-1, m, S, F, States_red, sym_computation, pi_red);

            comp = sum(constants{ss}(:,:,ii)) * (1 + sum(sum(delayed_hits_ratios_red{ss}))) * miss_ratios(ii);

            err = max(err, abs(del_hit - comp));

        end
    end
end

fprintf("Max computed error on (23): %d\n", err)

%% Exact latency
% (8)
di = compute_di(n, h, p_access, States, Q, pi, miss_ratios, hit_ratios);

% (9)
dh_ratios = zeros(1,n);
for ii = 1:n
    for ss = 1:S 
        dh_ratios(ii) = dh_ratios(ii) + sum(delayed_hits_ratios{ss}(:,:,ii));
    end
end
mean_req_retrieval_system = dh_ratios + di;

% (10)
latency = sum(mean_req_retrieval_system)/sum(p_access);

%% Approximated latency (with FPI)
parfor item = 1:n
    phi_item = sum(phi(item,:));
    mean_number_awaiting_approx(item) = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, item, phi_item);
end

latency_approx = sum(sum(phi',1) + mean_number_awaiting_approx) / (sum(p_access));

%% Real trace simulation
clear all
clc

n = 1000;

trace_folder = "simulator/reduced_chunks_1000_items/";
file_names = "chunks";
last_idx = 2655;

m = [n/20, n/20];
replacement_policy = "rr";
S = 3;
F = [1,1,1];
station_type = ["INF", "PS", "PS"];

h = numel(m);
F = F(1:S);

tot_time = 1814399;
load("simulator/reduced_chunks_1000_items/popularities_reduced_1000_items.mat"); lambdas = double(popularities) / tot_time; lambdas = lambdas(1:n);
load("simulator/reduced_chunks_1000_items/sizes_reduced_1000_items.mat"); times = double(sizes) / (10 * 1024 * 1024); times = times(1:n);

for ii = 1:n
    if mod(ii, 1e5) == 0
        fprintf("1, %d \n", ii)
    end
    for ss = 1:S
        alpha{ss}(:,:,ii) = 1;
    end
    T{1}(:,:,ii) = -1/(times(ii) * 10);
    T{2}(:,:,ii) = -1/(times(ii) * 1);
    T{3}(:,:,ii) = -1/(times(ii) * 1);
    R(:,:,ii) = [0.5, 0.25, 0.25 0; 0, 0, 0, 1; 0.1, 0, 0.3, 0.6; 0.1, 0.3, 0, 0.6];
end

factors = [10, 1, 1];

alpha_sim = cell(S, n);
T_sim = cell(S, n);
for ii = 1:n
    if mod(ii, 1e5) == 0
        fprintf("2, %d \n", ii)
    end
    for s = 1:S
        alpha_sim{s, ii} = alpha{s}(:,:,ii);
        T_sim{s, ii} = T{s}(:,:,ii);
    end
end

access_probabilities = ones(n, h);
access_fetching = ones(n, 1);

R_sim = [0.5, 0.25, 0.25 0; 0, 0, 0, 1; 0.1, 0, 0.3, 0.6; 0.1, 0.3, 0, 0.6];

fprintf("Events: %d \n", sum(popularities(1:n)));

recover = 0;
[miss_ratios, delayed_hits_ratios, hit_ratios] = simulate_trace_mixed_real_trace(trace_folder, file_names, last_idx, n, m, replacement_policy, S, F, station_type, access_probabilities, access_fetching, times, factors, R_sim, recover);

