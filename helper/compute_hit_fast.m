function hit_ratio = compute_hit_fast(x, n, m, replacement_policy, S, F, station_type, p_access, constants_without_pfetch)

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

    norm_const = recursive_bis_mixed(n, m, S, F, p_access, constants, [1:n], station_type, access_probabilities, access_fetching, false, zeros(1,S), zeros(1,n), replacement_policy);
    
    hit_ratios = zeros(n,h);
    for item_idx = 1:n

        p_access_red = [];
        constants_red = [];
        idxs_red = [];
        
        for kk = 1:n
            if kk ~= item_idx
                p_access_red = [p_access_red p_access(kk)];
                idxs_red = [idxs_red kk];
            end
        end
    
        for ss = 1:S
            constants_red{ss} = constants{ss};
            constants_red{ss}(:,:,item_idx) = [];
        end
    
        access_probabilities_red = access_probabilities;
        access_probabilities_red(item_idx, :) = [];
        access_fetching_red = access_fetching;
        access_fetching_red(item_idx) = [];
        
        PS_values = zeros(1,S);

        for list_num = 1:length(m)
            computed_hit = 1;
            for g = 1:list_num
                computed_hit = computed_hit * access_probabilities(item_idx, g);
            end
            computed_hit = computed_hit * p_access(item_idx)^list_num * access_fetching(item_idx) * recursive_mixed(n-1, oner(m,list_num), S, F, p_access_red, constants_red, idxs_red, station_type, access_probabilities_red, access_fetching_red, false, PS_values, replacement_policy) / norm_const;
            
            if isequal(replacement_policy, 'fifo')
                computed_hit = computed_hit * m(list_num);
            end
            
            hit_ratios(item_idx, list_num) = computed_hit;
        end

    end

    hit_ratio = 1 - sum(hit_ratios' * p_access');

end


