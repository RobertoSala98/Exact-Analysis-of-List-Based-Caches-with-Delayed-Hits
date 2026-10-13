function sequences = generate_combinations(n, m)

    total_len = sum(m);
    elements = 1:n;
    
    if total_len > n
        error('Not enough elements in the set to partition without repetition.');
    end
    
    base_combos = nchoosek(elements, total_len);
    sequences = [];
    for i = 1:size(base_combos, 1)
        if mod(i, 10) == 0
            fprintf('');
        end
        combo = base_combos(i, :);
        all_perms = perms(combo);
        for j = 1:size(all_perms, 1)
            p = all_perms(j, :);
            split_seq = mat2cell(p, 1, m);
            sequences = [sequences; split_seq];
        end
    end

end