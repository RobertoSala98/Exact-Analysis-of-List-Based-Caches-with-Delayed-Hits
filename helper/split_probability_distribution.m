function [miss_ratios, delayed_hits_ratios, hit_ratios] = split_probability_distribution(n, m, S, F, States, sym_computation, pi)

    if sym_computation
        hit_ratios = sym(zeros(length(m),n));
        delayed_hits_ratios = [];
        for ss = 1:S
            delayed_hits_ratios = [delayed_hits_ratios, {sym(zeros(1,F(ss),n))}];
        end
        miss_ratios = sym(zeros(1,n));
    else
        hit_ratios = zeros(length(m),n);
        delayed_hits_ratios = [];
        for ss = 1:S
            delayed_hits_ratios = [delayed_hits_ratios, {zeros(1,F(ss),n)}];
        end
        miss_ratios = zeros(1,n);
    end
    
    s = size(States,1);
    h = numel(m);

    for i=1:s
        state = States(i,:);

        elem_in_cache = [state{1:h}];
        fetch_index = [find(state{h+1} > 0)];
        all_elems = [fetch_index elem_in_cache];
        out_in_cache = setdiff(1:n, all_elems);

        for l=1:h
            for k=1:m(l)
                hit_ratios(l,state{l}(k)) = hit_ratios(l,state{l}(k)) + pi(i);
            end
        end
    
        for j=1:length(fetch_index)
            if state{h+1}(fetch_index(j)) > 0
                delayed_hits_ratios{state{h+1}(fetch_index(j))}(1, state{h+2}(fetch_index(j)), fetch_index(j)) = delayed_hits_ratios{state{h+1}(fetch_index(j))}(1, state{h+2}(fetch_index(j)), fetch_index(j)) + pi(i);
            end
        end
    
        for j=1:length(out_in_cache)
            miss_ratios(out_in_cache(j)) = miss_ratios(out_in_cache(j)) + pi(i);
        end
    end