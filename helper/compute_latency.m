function latency = compute_latency(x, n, m, replacement_policy, S, F, station_type, p_access, alpha, T, R, constants_without_pfetch)
    
    h = numel(m);
    access_probabilities = zeros(n, h);
    access_fetching = ones(n,1);

    for hh = 1:h
        access_probabilities(:,hh) = x(n*(hh-1)+1:n*hh);
    end

    constants = cell(1,S);
    for ss = 1:S
        for ii = 1:n
            for ff = 1:F(ss)
                constants{ss}(1,ff,ii) = constants_without_pfetch{ss}(1,ff,ii);
            end
        end
    end
   
    %{
    norm_const = recursive_bis_mixed(n, m, S, F, p_access, constants, [1:n], station_type, access_probabilities, access_fetching, false, zeros(1,S), zeros(1,n), replacement_policy);
    
    [miss_ratios, delayed_hits_ratios, hit_ratios] = compute_probability_distribution_fast_mixed(n, m, S, F, p_access, constants, station_type, access_probabilities, access_fetching, replacement_policy);

    del_hit_ratios = zeros(1,n);
    for ii = 1:n
        for ss = 1:S
            del_hit_ratios(ii) = del_hit_ratios(ii) + sum(delayed_hits_ratios{ss}(:,:,ii));
        end
    end

    latency = sum(del_hit_ratios ./ miss_ratios .* (del_hit_ratios + miss_ratios));
    %}

    States = generate_states(n, m, replacement_policy, S, F);
    Q = generate_QMatrix_mixed(n, m, replacement_policy, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, States, false);
    
    [pi, norm_const] = compute_probability_distribution_mixed(n, m, S, p_access, station_type, access_probabilities, access_fetching, States, false, constants);
    [miss_ratios, delayed_hits_ratios, hit_ratios] = split_probability_distribution(n, m, S, F, States, false, pi);
    pi = pi';

    for item = 1:n
        A1 = zeros(size(Q,1));
        
        nS = size(States,1);
        
        in_cache = false(nS,1);
        in_fetch = false(nS,1);
        
        for s = 1:nS
            state_s = States(s,:);
        
            for mm = 1:h
                if any(state_s{mm} == item)
                    in_cache(s) = true;
                    break
                end
            end
        
            if ~in_cache(s) && state_s{h+1}(item) > 0
                in_fetch(s) = true;
            end
        end
        
        cache_idx = find(in_cache);
        fetch_idx = find(in_fetch);
        used_idx  = union(cache_idx, fetch_idx);
        miss_idx  = setdiff(1:nS, used_idx);
        
        [~, loc_fetch] = ismember(fetch_idx, used_idx);
        A1 = sparse(Q(miss_idx, used_idx));
        
        pi_e = (pi(miss_idx) / (p_access(item) * miss_ratios(item)));
        pi_e = pi_e * A1;
        
        D0 = Q(fetch_idx, fetch_idx);
        
        first_moment = sum(pi_e(loc_fetch) * (-D0)^(-1));
        second_moment = sum(2 * pi_e(loc_fetch)  * (-D0)^(-2));
        
        phi_item = 1 - miss_ratios(item) - sum(hit_ratios(:,item));
        
        d_item(item) = phi_item * p_access(item) * second_moment / (2 * first_moment);
    end

    phi_ratios = 1 - (miss_ratios + sum(hit_ratios,1));
    latency = sum(phi_ratios + d_item) / sum(p_access);

end


