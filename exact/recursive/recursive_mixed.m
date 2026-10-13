function norm_const = recursive_mixed(n, m, S, F, p_access, eta, idxs, station_type, access_probabilities, access_fetching, sym_computation, PS_values, replacement_policy)

    h=numel(m);

    if sum(m)>n || min(m)<0
        norm_const=0;
        return
    end

    if n==0
        norm_const = 1;
        return
    end

    kk = randi(length(idxs));

    p_access_without_k = [];
    eta_without_k = [];

    access_probabilities_without_k = access_probabilities;
    access_probabilities_without_k(kk, :) = [];
    access_fetching_without_k = access_fetching;
    access_fetching_without_k(kk) = [];

    for ss = 1:S
        eta_without_k{ss} = eta{ss};
        eta_without_k{ss}(:,:,kk) = [];
    end

    idxs_without_k = [];
    for jj = 1:n
        if jj ~= kk
            p_access_without_k = [p_access_without_k p_access(jj)];
            idxs_without_k = [idxs_without_k idxs(jj)];
        end
    end
    
    norm_const = recursive_mixed(n-1, m, S, F, p_access_without_k, eta_without_k, idxs_without_k, station_type, access_probabilities_without_k, access_fetching_without_k, sym_computation, PS_values, replacement_policy);

    for ss = 1:S
        PS_values_updated = PS_values;
        if isequal(station_type(ss), "PS")
            PS_values_updated(ss) = PS_values_updated(ss) + 1;
        else
            PS_values_updated(ss) = 1;
        end

        in_fetching_recursive_term = recursive_mixed(n-1, m, S, F, p_access_without_k, eta_without_k, idxs_without_k, station_type, access_probabilities_without_k, access_fetching_without_k, sym_computation, PS_values_updated, replacement_policy);
        norm_const = norm_const + sum(eta{ss}(:,:,kk)) * max(1, PS_values_updated(ss)) * in_fetching_recursive_term;
    end

    for j=1:h
        if m(j)>0
            recursive_term = 1;
            for g = 1:j
                recursive_term = recursive_term * access_probabilities(kk, g);
            end
            recursive_term = recursive_term * p_access(kk)^j * access_fetching(kk) * recursive_mixed(n-1, oner(m,j), S, F, p_access_without_k, eta_without_k, idxs_without_k, station_type, access_probabilities_without_k, access_fetching_without_k, sym_computation, PS_values, replacement_policy);
            
            if isequal(replacement_policy, 'fifo')
                recursive_term = recursive_term * m(j);
            end

            norm_const = norm_const + recursive_term;
        end
    end
end