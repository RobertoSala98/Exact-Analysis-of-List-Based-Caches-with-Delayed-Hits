function sequences = generate_disjoint_combinations(n, m)
    total_len = sum(m);
    if total_len > n
        error('Not enough elements to create disjoint groups.');
    end

    elements = 1:n;
    sequences = generate_recursive(elements, m);
end

function result = generate_recursive(pool, m)
    if isempty(m)
        result = {{}};
        return;
    end

    k = m(1);
    rest_m = m(2:end);

    this_combos = nchoosek(pool, k);
    result = {};

    for i = 1:size(this_combos, 1)
        group = this_combos(i, :);
        remaining = setdiff(pool, group);

        sub_results = generate_recursive(remaining, rest_m);

        for j = 1:length(sub_results)
            result{end+1, 1} = [{group}, sub_results{j}];
        end
    end
end
