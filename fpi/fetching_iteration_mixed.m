function [phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type, max_iter, tol, verbose)
    % FPI for delayed-hit caches
    % Inputs:
    %   gamma: n x h,  eta: n x S,  m: h x 1,  station_type: S x 1 string ("INF"/"PS")
    % Outputs:
    %   phi: n x S,  xi: 1 x h,  pi_vals: n x h,  pi_0: n x 1

    if nargin < 5 || isempty(max_iter)
        max_iter = 1e5; 
    end
    
    if nargin < 6 || isempty(tol)      
        tol = 1e-6;     
    end

    if nargin < 7                      
        verbose = false; 
    end

    [n, h] = size(gamma);
    S = size(eta, 2);
    m = m(:)';
    station_type = string(station_type(:)');
    isINF = station_type == "INF";
    isPS  = station_type == "PS";

    % ---- initialization ----
    pi_vals = ones(n, h) .* m/n;
    pi_0 = (1 - sum(pi_vals,2)) / (S+1);
    phi = repmat(pi_0, 1, S);
    xi = m ./ max(sum(gamma,1), eps);
    xi = max(xi, tol);

    for iter = 1:max_iter
        phi_prev = phi; 
        pi_prev = pi_vals; 
        pi0_prev = pi_0; 
        xi_prev = xi;

        % ---- beta_i = 1 + lam*eta_0 + sum_s lam*eta_s*(1 + sum_{k~=i} phi_{s,k}) ----
        G = 1 + sum(phi_prev(:,isPS), 1) - phi_prev(:,isPS);      % n x #PS
        beta = 1 + sum(eta(:,isINF), 2) + sum(eta(:,isPS) .* G, 2);

        % ---- theta ----
        theta = gamma ./ beta;                                    % n x h

        % ---- xi: solve m_j = sum_i pi_ij ----
        xi = solve_xi(theta, m, xi_prev);

        % ---- pi_ij ----
        txi     = theta .* xi;                                    % n x h
        pi_vals = txi ./ (1 + sum(txi, 2));

        % ---- pi_i0 ----
        pi_0 = (1 - sum(pi_vals, 2)) ./ beta;

        % ---- phi ----
        phi(:,isINF) = eta(:,isINF) .* pi_0;
        phi(:,isPS)  = eta(:,isPS) .* G .* pi_0;

        % ---- convergence (max relative change) ----
        rel = @(a,b) max(abs(a(:)-b(:)) ./ max(abs(b(:)), eps));
        max_diff = max([rel(phi,phi_prev), rel(pi_vals,pi_prev), rel(pi_0,pi0_prev)]);

        if verbose && (mod(iter,10)==0 || iter<=5)
            fprintf('Iter %d: max rel change = %.3e\n', iter, max_diff);
        end
        if max_diff < tol
            if verbose, fprintf('Converged after %d iterations\n', iter); end
            break;
        end
    end
    if verbose && iter == max_iter
        fprintf('Maximum iterations reached. Max change: %.2e\n', max_diff);
    end

    if verbose
        fprintf('Probability sums per item (should be 1):\n');
        disp(max(abs(sum(phi,2) + sum(pi_vals,2) + pi_0 - 1)));
    end
end

% =====================================================================
function xi = solve_xi(theta, m, xi0)
    % Damped Newton on u = log(xi) for F_j(xi) = xi_j*sum_i theta_ij/D_i - m_j,
    % with D_i = 1 + sum_l theta_il*xi_l.
    m = m(:);  u = log(xi0(:));
    resid = @(u) local_F(theta, m, exp(u));
    [F, ~] = resid(u);
    for it = 1:200
        if norm(F, inf) < 1e-12*max(1, max(m)), break; end
        xi = exp(u);
        D  = 1 + theta*xi;                                        % n x 1
        A  = theta' * (theta ./ D.^2);                            % h x h
        Jxi = diag(sum(theta ./ D, 1)) - xi .* A;                 % dF_j/dxi_l
        Ju  = Jxi .* xi';                                         % dF_j/du_l
        du  = -(Ju \ F);
        t = 1;                                                    % backtracking line search
        while true
            un = u + t*du;
            Fn = resid(un);
            if norm(Fn) < norm(F) || t < 1e-10, break; end
            t = t/2;
        end
        u = un;  F = Fn;
    end
    xi = exp(u)';                                                 % 1 x h
end

function [F, D] = local_F(theta, m, xi)
    D = 1 + theta*xi;
    F = xi .* sum(theta ./ D, 1)' - m;
end