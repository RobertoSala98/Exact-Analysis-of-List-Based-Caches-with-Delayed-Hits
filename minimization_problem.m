% Example in section 6.1

%% Setup the environment
clear all
clc

if isempty(gcp('nocreate'))
    pool = parpool(IdleTimeout=1440);
end

addpath(genpath(pwd));

%% Parameters setting
n = 5;
m = [3];
replacement_policy = 'rr';
S = 2;
F = [1,1];
station_type = ["INF", "PS"];
sym_computation = false;

h = numel(m);
F = F(1:S);
station_type = station_type(1:S);

replacement_policy = 'rr';

zipf_rate = 1.4;
ranks = 1:n;
H = sum(1 ./ (ranks.^zipf_rate));
p_access = (1 ./ (ranks.^zipf_rate)) / H;

alpha = [{}, {}];
T = [{}, {}];
for ii = 1:n
    for ss = 1:S
        alpha{ss}(:,:,ii) = 1;
    end
    R(:,:,ii) = 1/(1+S) * ones(S+1, S+1);
    R(1,:,ii) = [1/S*ones(1,S), 0];
end

T{1}(:,:,1) = -1.0;
T{1}(:,:,2) = -1.0;
T{1}(:,:,3) = -1.0;
T{1}(:,:,4) = -0.01;
T{1}(:,:,5) = -0.01;
T{2}(:,:,1) = -1.0;
T{2}(:,:,2) = -1.0;
T{2}(:,:,3) = -1.0;
T{2}(:,:,4) = -0.01;
T{2}(:,:,5) = -0.01;

access_probabilities_init = ones(n, h);
access_fetching = ones(n, 1);

constants_without_pfetch = compute_constants(n, S, F, p_access, alpha, T, R, false, ones(n,1));

%% Optimal hit ratio and latency
x0 = [];
for hh = 1:h
    x0 = [x0; access_probabilities_init(:,hh)];
end

fprintf("Starting non-Hit ratio with policy %s: %.6f \n", replacement_policy, compute_hit_fast(x0, n, m, replacement_policy, S, F, station_type, p_access, constants_without_pfetch))

fun = @(x) compute_hit_fast(x, n, m, replacement_policy, S, F, station_type, p_access, constants_without_pfetch);

lb = zeros(1, n*h);
ub = ones(1, n*h);

options = optimoptions('fmincon','MaxIterations',10000,'OptimalityTolerance',1e-10,'UseParallel',true,'MaxFunctionEvaluations',1e5);
[x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

access_probabilities_opt_hit = zeros(n, h);
for hh = 1:h
    access_probabilities_opt_hit(:,hh) = x(n*(hh-1)+1:n*hh);
end

access_probabilities_opt_hit

fprintf("Final latency %s: %.16f \n", replacement_policy, compute_latency(x, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch))
fprintf("Final non-Hit ratio: %.16f \n", compute_hit_fast(x, n, m, replacement_policy, S, F, station_type, p_access, constants_without_pfetch))

%
fprintf("\n")

x0 = [];
for hh = 1:h
    x0 = [x0; access_probabilities_init(:,hh)];
end

fprintf("Starting latency %s: %.6f \n", replacement_policy, compute_latency(x0, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch))

fun = @(x) compute_latency(x, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch);

lb = zeros(1, n*h);
ub = ones(1, n*h);

options = optimoptions('fmincon','MaxIterations',10000,'OptimalityTolerance',1e-10,'UseParallel',true,'MaxFunctionEvaluations',1e5);
[x,fval] = fmincon(fun,x0,[],[],[],[],lb,ub,[],options);

access_probabilities_opt_latency = zeros(n, h);
for hh = 1:h
    access_probabilities_opt_latency(:,hh) = x(n*(hh-1)+1:n*hh);
end

access_probabilities_opt_latency

fprintf("Final latency %s: %.16f \n", replacement_policy, compute_latency(x, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch))
fprintf("Final non-Hit ratio: %.16f \n", compute_hit_fast(x, n, m, replacement_policy, S, F, station_type, p_access, constants_without_pfetch))

%% Scenarios 1 and 2

numEv = 1e4;
hit_rates = zeros(1,numEv);
latencies = zeros(1,numEv);

parfor Ev = 1:numEv
    access_probabilities = rand(n,h);
    hit_rates(Ev) = 1 - compute_hit_fast(access_probabilities, n, m, ...
        replacement_policy, S, F, station_type, p_access, constants_without_pfetch);
    latencies(Ev) = compute_latency(access_probabilities, n, m, ...
        replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch);
end

figure;
scatter(hit_rates, latencies, 10, [0.6 0.6 0.6], 'x');  % grey, smaller size
xlabel('Hit rate', 'FontSize', 24); 
ylabel('Latency', 'FontSize', 24);
xlim([0.2, 1.0])
hold on

% Compute special points
x1 = 1 - compute_hit_fast(access_probabilities_opt_hit, n, m, ...
    replacement_policy, S, F, station_type, p_access, constants_without_pfetch);
y1 = compute_latency(access_probabilities_opt_hit, n, m, ...
    replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch);

x2 = 1 - compute_hit_fast(access_probabilities_opt_latency, n, m, ...
    replacement_policy, S, F, station_type, p_access, constants_without_pfetch);
y2 = compute_latency(access_probabilities_opt_latency, n, m, ...
    replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch);

x3 = 1 - compute_hit_fast(access_probabilities_init, n, m, ...
    replacement_policy, S, F, station_type, p_access, constants_without_pfetch);
y3 = compute_latency(access_probabilities_init, n, m, ...
    replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch);

% Plot the special points
scatter(x1, y1, 160, '+', 'MarkerEdgeColor','red', 'LineWidth',3.0)
scatter(x2, y2, 160, 'x', 'MarkerEdgeColor','blue', 'LineWidth',3.0)
scatter(x3, y3, 160, 'x', 'MarkerEdgeColor','magenta', 'LineWidth',3.0)

% Add labels with bigger fonts
text(x1-0.05, y1, 'Max Hit Rate  ', 'Color','red', 'FontSize',24, ...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
text(x2, y2, 'Min Latency ', 'Color','blue', 'FontSize',24, ...
    'VerticalAlignment','bottom','HorizontalAlignment','right');
text(x3-0.05, y3, 'Initial value ', 'Color','magenta', 'FontSize',24, ...
    'VerticalAlignment','top','HorizontalAlignment','right');

% Optional: set axis tick font size larger
set(gca, 'FontSize', 24);

hold off
grid on;

fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'hit_latency_plot_fast_less_popular.png', '-dpng', '-r300');

%% Scenario 3

item = 1;

numEv = 100;
T_ = 10.^linspace(-2,2,numEv);

best_access_hit = zeros(1,numEv);
best_access_lat = zeros(1,numEv);

hit_rates = zeros(1,numEv);
latencies = zeros(1,numEv);

time_h = [];
time_l = [];

lb = zeros(1, n*h);
ub = ones(1, n*h);

options = optimoptions('fmincon', ...
    'MaxIterations', 10000, ...
    'OptimalityTolerance', 1e-10, ...
    'UseParallel', true, ...
    'MaxFunctionEvaluations', 1e5, ...
    'Display','off'); % to suppress verbose output

% Multi-start configuration
ms = MultiStart('UseParallel', true, 'Display', 'off');
num_starts = 20; % number of random starting points (adjust as needed)

for Ev = 1:numEv
    % Update simulation constants for this iteration
    T_sim = T;
    T_sim{1}(:,:,item) = -T_(Ev);
    T_sim{2}(:,:,item) = -T_(Ev);

    constants = compute_constants(n, S, F, p_access, alpha, T_sim, R, false, ones(n,1));

    % Objective functions
    hit_fun = @(x) compute_hit_fast(x, n, m, replacement_policy, S, F, station_type, p_access, constants);
    lat_fun = @(x) compute_latency(x, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants);

    % Define optimization problems
    problem_hit = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', hit_fun, ...
        'lb', lb, 'ub', ub, 'options', options);

    problem_lat = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', lat_fun, ...
        'lb', lb, 'ub', ub, 'options', options);

    % Run MultiStart
    tic
    [x_hit, fval_hit] = run(ms, problem_hit, num_starts);
    el = toc;
    time_h = [time_h, el]

    tic
    [x_lat, fval_lat] = run(ms, problem_lat, num_starts);
    el = toc;
    time_l = [time_l, el]

    % Store results
    best_access_hit(Ev) = x_hit(item);
    best_access_lat(Ev) = x_lat(item);

    hit_rates(Ev) = 1 - fval_hit;
    latencies(Ev) = fval_lat;

    if Ev == 10
        mean(time_h)
        mean(time_l)
    end
end

% === Parameters ===
windowSize = 10; % Adjust for mean average smoothing (set to 1 to disable)

% Compute mean average (moving average)
best_access_hit_avg = movmean(best_access_hit, windowSize);
best_access_lat_avg = movmean(best_access_lat, windowSize);

% === Plotting ===
figure
yyaxis left
semilogx(T_, best_access_hit_avg, '-or', 'LineWidth', 1.5)
ylabel('$a_{1,1}$ for optimal hit rate', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

ax = gca; % current axes
ax.YColor = 'r'; % left y-axis color = red

yyaxis right
semilogx(T_, best_access_lat_avg, '-ob', 'LineWidth', 1.5)
ylabel('$a_{1,1}$ for optimal latency', 'FontSize', 24, 'Interpreter', 'latex')
ylim([0 1])

ax.YColor = 'b'; % right y-axis color = blue

xlabel('$\mu_1$', 'FontSize', 24, 'Interpreter', 'latex')
set(gca, 'FontSize', 24)

grid on

% === Save figure ===
fig = gcf;
set(fig, 'PaperPositionMode', 'auto');
print(fig, 'optimal_access_probabilities_item1.png', '-dpng', '-r300');
