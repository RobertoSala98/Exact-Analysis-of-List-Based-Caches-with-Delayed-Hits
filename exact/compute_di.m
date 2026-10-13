function mean_awaiting_completition = compute_di(n, h, p_access, States, Q, pi, miss_ratios, hit_ratios)

    parfor item = 1:n
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
    
        mean_awaiting_completition(item) = phi_item * p_access(item) * second_moment / (2 * first_moment);
    end

end