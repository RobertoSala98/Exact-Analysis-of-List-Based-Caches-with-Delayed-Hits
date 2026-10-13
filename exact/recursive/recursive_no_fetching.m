function norm_const = recursive_no_fetching(n, m, replacement_policy, p_access, access_probabilities, access_fetching)

    h=length(m);

    if sum(m)==0
        norm_const = 1.0;
        return
    end

    if sum(m)>n || min(m)<0
        norm_const = 0.0;
        return
    end

    if n==1 && sum(m)==1
        norm_const = p_access(1)^find(m) * access_fetching(1);
        for g = 1:find(m)
            norm_const = norm_const * access_probabilities(1, g);
        end
        return
    end

    kk = randi(n);

    p_access_without_k = p_access;
    p_access_without_k(kk) = [];

    access_probabilities_without_k = access_probabilities;
    access_probabilities_without_k(kk, :) = [];
    access_fetching_without_k = access_fetching;
    access_fetching_without_k(kk) = [];
    
    norm_const = recursive_no_fetching(n-1, m, replacement_policy, p_access_without_k, access_probabilities_without_k, access_fetching_without_k);
    
    for j=1:h
        if m(j)>0
            recursive_term = 1;
            for g = 1:j
                recursive_term = recursive_term * access_probabilities(kk, g);
            end
            recursive_term = recursive_term * p_access(kk)^j * access_fetching(kk) * recursive_no_fetching(n-1, oner(m,j), replacement_policy, p_access_without_k, access_probabilities_without_k, access_fetching_without_k);

            if isequal(replacement_policy, 'fifo')
                recursive_term = recursive_term * m(j);
            end
            norm_const = norm_const + recursive_term;
        end
    end

end