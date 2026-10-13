function max_utilization = compute_max_utilization(x, n, m, S, station_type, p_access, constants)

    h = numel(m);
    access_probabilities = zeros(n, h);

    for hh = 1:h
        access_probabilities(:,hh) = x(n*(hh-1)+1:n*hh);
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


