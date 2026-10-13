%% Setup the environment
clear all
clc

if isempty(gcp('nocreate'))
    pool = parpool(IdleTimeout=1440);
end

addpath(genpath(pwd));

%% System Parameters
n = 100;
m = [10,10];
replacement_policy = 'fifo';
S = 2;
F = [1,1];
station_type = ["PS", "PS"];
sym_computation = false;

h = numel(m);
F = F(1:S);

zipf_rate = 1.0;
ranks = 1:n;
H = sum(1 ./ (ranks.^zipf_rate));
p_access = (1 ./ (ranks.^zipf_rate)) / H;

alpha = [{}, {}];
T = [{}, {}];
for ii = 1:n
    for ss = 1:S
        alpha{ss}(:,:,ii) = 1;
    end

    if ii <= sum(m)
        T{1}(:,:,ii) = -1.0;
        T{2}(:,:,ii) = -0.01;
        R(:,:,ii) = [0.75, 0.25, 0; 0, 0, 1; 0, 0, 1];
    else
        T{1}(:,:,ii) = -0.01;
        T{2}(:,:,ii) = -0.01;
        R(:,:,ii) = [0.25, 0.75, 0; 0, 0, 1; 0, 0, 1];
    end
end

access_probabilities = ones(n, h);
access_fetching = ones(n, 1);

constants = compute_constants(n, S, F, p_access, alpha, T, R, sym_computation, access_fetching);

%% Initial scenario

gamma = zeros(n,h);
for ii = 1:n
    for hh = 1:h
        gamma(ii,hh) = p_access(ii)^hh * access_fetching(ii);
        for g = 1:hh
            gamma(ii,hh) = gamma(ii,hh) * access_probabilities(ii,g);
        end
    end
end

eta_tmp =[];
for ss = 1:S
    for kk = 1:n
        eta_tmp(ss,kk) = sum(constants{ss}(:,:,kk));
    end
end
eta = eta_tmp';

[phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type);

parfor item = 1:n
    phi_item = sum(phi(item,:));
    mean_number_awaiting_approx(item) = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, item, phi_item);
end

latency = sum(sum(phi',1) + mean_number_awaiting_approx) / (sum(p_access))

occupancies = sum(phi)

%% Stage 1

x0 = [];
for hh = 1:h
    x0 = [x0; access_probabilities(:,hh)];
end

fun = @(x) compute_max_utilization(x, n, m, S, station_type, p_access, constants);

lb = zeros(n*h, 1);  % Column vector
ub = ones(n*h, 1);   % Column vector

x0 = reshape(access_probabilities, [], 1);  % n*h x 1

options = optimoptions('fmincon', ...
    'MaxIterations', 1e6, ...
    'OptimalityTolerance', 1e-6, ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 1e6, ...
    'Display','iter');

problem = createOptimProblem('fmincon', ...
    'objective', fun, ...
    'x0', x0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

ms = MultiStart('StartPointsToRun','all','Display','iter');

nStartPoints = 10;  % You can increase for better global search
startPoints = RandomStartPointSet('NumStartPoints', nStartPoints);

[x,fval] = run(ms, problem, startPoints);

access_probabilities_util = zeros(n, h);
for hh = 1:h
    access_probabilities_util(:,hh) = x(n*(hh-1)+1:n*hh);
end

gamma = zeros(n,h);
for ii = 1:n
    for hh = 1:h
        gamma(ii,hh) = p_access(ii)^hh * access_fetching(ii);
        for g = 1:hh
            gamma(ii,hh) = gamma(ii,hh) * access_probabilities_util(ii,g);
        end
    end
end

eta_tmp =[];
for ss = 1:S
    for kk = 1:n
        eta_tmp(ss,kk) = sum(constants{ss}(:,:,kk));
    end
end
eta = eta_tmp';
%
[phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type);

del_hit_ratios = sum(phi,2);

parfor item = 1:n
    phi_item = sum(phi(item,:));
    mean_number_awaiting_approx(item) = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R, station_type, access_probabilities_util, access_fetching, item, phi_item);
end

latency_util = sum(sum(phi',1) + mean_number_awaiting_approx) / (sum(p_access))

occupancies_util = sum(phi)

%% Stage 2

lb = zeros(n*h + n, 1);
ub = ones(n*h + n, 1);

x0 = reshape(access_probabilities, [], 1);

idx = numel(x0);
for ii = 1:n
    idx = idx + 1;
    x0(idx,1) = R(1,1,ii);  % Append R values
end

fun = @(x) compute_max_utilization_with_R(x, n, m, S, F, station_type, p_access, alpha, T);

options = optimoptions('fmincon', ...
    'MaxIterations', 1e6, ...
    'OptimalityTolerance', 1e-6, ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 1e6, ...
    'Display','iter');  % Optional, can override display

problem = createOptimProblem('fmincon', ...
    'objective', fun, ...
    'x0', x0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

ms = MultiStart('StartPointsToRun','all','Display','iter');

nStartPoints = 10;
startPoints = RandomStartPointSet('NumStartPoints', nStartPoints);

[x,fval] = run(ms, problem, startPoints);

access_probabilities_util_r = zeros(n, h);
R_util_r = zeros(S+1, S+1, n);
R_util_r(2:S+1, S+1, :) = 1;

for hh = 1:h
    access_probabilities_util_r(:,hh) = x(n*(hh-1)+1:n*hh);
end

idx = prod(size(access_probabilities_util_r));

for ii = 1:n
    idx = idx + 1;
    R_util_r(1, 1, ii) = x(idx);
    R_util_r(1, 2, ii) = 1 - x(idx);
end

%
gamma = zeros(n,h);
for ii = 1:n
    for hh = 1:h
        gamma(ii,hh) = p_access(ii)^hh * access_fetching(ii);
        for g = 1:hh
            gamma(ii,hh) = gamma(ii,hh) * access_probabilities_util_r(ii,g);
        end
    end
end

constants = compute_constants(n, S, F, p_access, alpha, T, R_util_r, false, ones(1,n));

eta_tmp =[];
for ss = 1:S
    for kk = 1:n
        eta_tmp(ss,kk) = sum(constants{ss}(:,:,kk));
    end
end
eta = eta_tmp';
%
[phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type);

parfor item = 1:n
    phi_item = sum(phi(item,:));
    mean_number_awaiting_approx(item) = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R_util_r, station_type, access_probabilities_util_r, access_fetching, item, phi_item);
end

latency_util_r = sum(sum(phi',1) + mean_number_awaiting_approx) / (sum(p_access))

occupancies_util_r = sum(phi)

%% Stage 3

lb = zeros(n*h + n, 1);
ub = ones(n*h + n, 1);

x0 = reshape(access_probabilities, [], 1);  % Stack columns into a column vector

idx = numel(x0);
for ii = 1:n
    idx = idx + 1;
    x0(idx,1) = R(1,1,ii);  % Append R values
end

fun = @(x) compute_latency_FPI(x, n, m, S, F, station_type, p_access, alpha, T);

options = optimoptions('fmincon', ...
    'MaxIterations', 1e5, ...
    'OptimalityTolerance', 1e-4, ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 1e5, ...
    'Display','iter');

problem = createOptimProblem('fmincon', ...
    'objective', fun, ...
    'x0', x0, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', options);

ms = MultiStart('StartPointsToRun','all','Display','iter');

nStartPoints = 10;
startPoints = RandomStartPointSet('NumStartPoints', nStartPoints);

[x, fval] = run(ms, problem, startPoints);
%
access_probabilities_best_latency = zeros(n, h);
R_best_latency = zeros(S+1, S+1, n);
R_best_latency(2:S+1, S+1, :) = 1;

for hh = 1:h
    access_probabilities_best_latency(:,hh) = x(n*(hh-1)+1:n*hh);
end

idx = prod(size(access_probabilities_best_latency));
for ii = 1:n
    idx = idx + 1;
    R_best_latency(1, 1, ii) = x(idx);
    R_best_latency(1, 2, ii) = 1 - x(idx);
end

gamma = zeros(n,h);
for ii = 1:n
    for hh = 1:h
        gamma(ii,hh) = p_access(ii)^hh * access_fetching(ii);
        for g = 1:hh
            gamma(ii,hh) = gamma(ii,hh) * access_probabilities_best_latency(ii,g);
        end
    end
end

constants = compute_constants(n, S, F, p_access, alpha, T, R_best_latency, false, ones(1,n));

eta_tmp =[];
for ss = 1:S
    for kk = 1:n
        eta_tmp(ss,kk) = sum(constants{ss}(:,:,kk));
    end
end
eta = eta_tmp';
%
[phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type);

parfor item = 1:n
    phi_item = sum(phi(item,:));
    mean_number_awaiting_approx(item) = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R_best_latency, station_type, access_probabilities_best_latency, access_fetching, item, phi_item);
end

latency_best_latency = sum(sum(phi',1) + mean_number_awaiting_approx) / (sum(p_access))

occupancies_best_latency = sum(phi)

x_util_r = reshape(access_probabilities_util_r, [], 1);

idx = numel(x_util_r);
for ii = 1:n
    idx = idx + 1;
    x_util_r(idx,1) = R_util_r(1,1,ii);  % Append R values
end

%% Plots
RR_init = [];
RR_util_r = [];
RR_best_latency = [];

for ii = 1:n
    RR_init = [RR_init, R(1,1,ii)];
    RR_util_r = [RR_util_r, R_util_r(1,1,ii)];
    RR_best_latency = [RR_best_latency, R_best_latency(1,1,ii)];
end



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, access_probabilities_util(:,1), markerSize, 'x', 'LineWidth', 1.5); 
hold on; 
scatter(1:n, access_probabilities_util(:,2), markerSize, 'x', 'LineWidth', 1.5);

xlabel('item index', 'FontSize', 24)
ylabel('optimal $\mathrm{a}_{i,j}$ in stage 1', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

legend({'$a_{\cdot,1}$', '$a_{\cdot,2}$'}, 'FontSize', 24, 'Location', 'northeast', 'Interpreter', 'latex');

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_access_probabilities_stage1.png', '-dpng', '-r300');



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, access_probabilities_util_r(:,1), markerSize, 'x', 'LineWidth', 1.5); 
hold on; 
scatter(1:n, access_probabilities_util_r(:,2), markerSize, 'x', 'LineWidth', 1.5);

xlabel('item index', 'FontSize', 24)
ylabel('optimal $\mathrm{a}_{i,j}$ in stage 2', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

legend({'$a_{\cdot,1}$', '$a_{\cdot,2}$'}, 'FontSize', 24, 'Location', 'northeast', 'Interpreter', 'latex');

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_access_probabilities_stage2.png', '-dpng', '-r300');



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, access_probabilities_best_latency(:,1), markerSize, 'x', 'LineWidth', 1.5); 
hold on; 
scatter(1:n, access_probabilities_best_latency(:,2), markerSize, 'x', 'LineWidth', 1.5);

xlabel('item index', 'FontSize', 24)
ylabel('optimal $\mathrm{a}_{i,j}$ in stage 3', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

legend({'$a_{\cdot,1}$', '$a_{\cdot,2}$'}, 'FontSize', 24, 'Location', 'northeast', 'Interpreter', 'latex');

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_access_probabilities_stage3.png', '-dpng', '-r300');



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, RR_init, markerSize, 'x', 'LineWidth', 1.5); 

xlabel('item index', 'FontSize', 24)
ylabel('initial $\psi_{1,\cdot}$', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_R_stage1.png', '-dpng', '-r300');



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, RR_util_r, markerSize, 'x', 'LineWidth', 1.5); 

xlabel('item index', 'FontSize', 24)
ylabel('optimal $\psi_{1,\cdot}$ in stage 2', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_R_stage2.png', '-dpng', '-r300');



figure

markerSize = 100; % Adjust this number to make markers bigger

scatter(1:n, RR_best_latency, markerSize, 'x', 'LineWidth', 1.5); 

xlabel('item index', 'FontSize', 24)
ylabel('optimal $\psi_{1,\cdot}$ in stage 3', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

grid on
set(gca, 'FontSize', 24);

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_R_stage3.png', '-dpng', '-r300');