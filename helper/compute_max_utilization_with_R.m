function max_utilization = compute_max_utilization_with_R(x, n, m, S, F, station_type, p_access, alpha, T)

    h = numel(m);
    access_probabilities = zeros(n, h);
    R = zeros(S+1, S+1, n);

    for hh = 1:h
        access_probabilities(:,hh) = x(n*(hh-1)+1:n*hh);
    end
    
    idx = prod(size(access_probabilities));
    R = zeros(S+1, S+1, n);
    R(2:S+1, S+1, :) = 1;

    for ii = 1:n
        idx = idx + 1;
        R(1, 1, ii) = x(idx);
        R(1, 2, ii) = 1 - x(idx);
    end

    gamma = zeros(n,h);
    for ii = 1:n
        for hh = 1:h
            gamma(ii,hh) = p_access(ii)^hh;
            for g = 1:hh
                gamma(ii,hh) = gamma(ii,hh) * access_probabilities(ii,g);
            end
        end
    end

    constants = compute_constants(n, S, F, p_access, alpha, T, R, false, ones(1,n));
    
    eta_tmp =[];
    for ss = 1:S
        for kk = 1:n
            eta_tmp(ss,kk) = sum(constants{ss}(:,:,kk));
        end
    end
    eta = eta_tmp';

    [phi, xi, pi_vals, pi_0] = fetching_iteration_mixed(gamma, eta, m, station_type);

    max_utilization = max(sum(phi));

end


