function assignments = assign_values_to_vector(n, complement, k, h)
    if k > h
        error('Starting value k must be less than or equal to maximum value h.');
    end

    num_vars = length(complement);

    if num_vars == 0
        assignments = zeros(1, n);
        return;
    end

    % Generate all combinations from k to h for each variable
    ranges = repmat({k:h}, 1, num_vars);
    combos = cell(1, num_vars);
    [combos{:}] = ndgrid(ranges{:});
    val_matrix = zeros(numel(combos{1}), num_vars);

    for i = 1:num_vars
        val_matrix(:, i) = combos{i}(:);
    end

    num_combos = size(val_matrix, 1);
    assignments = zeros(num_combos, n);

    for i = 1:num_combos
        vec = zeros(1, n);
        vec(complement) = val_matrix(i, :);
        assignments(i, :) = vec;
    end
end
