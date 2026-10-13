function norm_const = compute_norm_const_mixed(n, m, S, gamma, eta, station_type, replacement_policy)

    H = compute_cache_erec(n, m, S, gamma, eta, m, [], station_type);
    norm_const = sum(H);

    if isequal(replacement_policy, 'rr')
        norm_const = norm_const / prod(factorial(m));
    end

end

function H = compute_cache_erec(n, m, S, gamma, eta, list_idx, H, station_type)
    if numel(list_idx) == numel(m) + S
        val = cache_erec([gamma, eta'], list_idx);
        for s = 1:S
            if isequal(station_type(s), "INF")
                val = val / factorial(list_idx(numel(m)+s));
            end
        end
        H(end+1) = val;
    else
        rem = n - sum(list_idx);
        for jj = 0:rem
            H = compute_cache_erec(n, m, S, gamma, eta, [list_idx, jj], H, station_type);
        end
    end
    return
end