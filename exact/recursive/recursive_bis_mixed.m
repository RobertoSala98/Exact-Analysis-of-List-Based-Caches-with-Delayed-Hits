function norm_const = recursive_bis_mixed(n, m, S, F, p_access, eta, idxs, station_type, access_probabilities, access_fetching, sym_computation, PS_values, non_fetching_terms, replacement_policy)

    if or(n == sum(m), sum(non_fetching_terms)==n)
        norm_const = recursive_no_fetching(n, m, replacement_policy, p_access, access_probabilities, access_fetching);
        return
    end

    if sum(m)>n || min(m)<0
        norm_const=0;
        return
    end
    
    non_fetching_idxs = find(non_fetching_terms == 0);
    kk = non_fetching_idxs(randi(length(non_fetching_idxs)));

    access_probabilities_without_k = access_probabilities;
    access_probabilities_without_k(kk, :) = [];
    access_fetching_without_k = access_fetching;
    access_fetching_without_k(kk) = [];

    eta_without_k = [];
    for ss = 1:S
        eta_without_k{ss} = eta{ss};
        eta_without_k{ss}(:,:,kk) = [];
    end

    p_access_without_k = [];
    idxs_without_k = [];
    non_fetching_terms_without_k = [];
    for jj = 1:n
        if jj ~= kk
            p_access_without_k = [p_access_without_k p_access(jj)];
            idxs_without_k = [idxs_without_k idxs(jj)];
            non_fetching_terms_without_k = [non_fetching_terms_without_k non_fetching_terms(jj)];
        end
    end

    non_fetching_terms_updated = non_fetching_terms;
    non_fetching_terms_updated(kk) = 1;

    norm_const = recursive_bis_mixed(n, m, S, F, p_access, eta, idxs, station_type, access_probabilities, access_fetching, sym_computation, PS_values, non_fetching_terms_updated, replacement_policy);

    for ss = 1:S
        PS_values_updated = PS_values;
        if isequal(station_type(ss), 'PS')
            PS_values_updated(ss) = PS_values_updated(ss) + 1;
        else
            PS_values_updated(ss) = 1;
        end

        in_fetching_recursive_term = recursive_bis_mixed(n-1, m, S, F, p_access_without_k, eta_without_k, idxs_without_k, station_type, access_probabilities_without_k, access_fetching_without_k, sym_computation, PS_values_updated, non_fetching_terms_without_k, replacement_policy);
        norm_const = norm_const + sum(eta{ss}(:,:,kk)) * max(1, PS_values_updated(ss)) * in_fetching_recursive_term;
    end
    
end