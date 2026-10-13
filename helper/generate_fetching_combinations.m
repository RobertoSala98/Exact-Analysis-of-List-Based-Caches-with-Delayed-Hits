function combos = generate_fetching_combinations(n, fetching_queue, F)
    % Check input dimensions
    if length(fetching_queue) ~= n
        error('fetching_queue must be of length n');
    end

    % Create cell array of ranges for each variable
    ranges = cell(1, n);
    for i = 1:n
        idx = fetching_queue(i);
        if idx == 0
            ranges{i} = 0;
        else
            ranges{i} = 1:F(idx);
        end
    end

    % Generate ndgrid and store into cell array 'grid'
    grid = cell(1, n);  % <- Fix: predefine the cell array
    [grid{:}] = ndgrid(ranges{:});

    % Convert grid to matrix of combinations
    combos = zeros(numel(grid{1}), n);
    for i = 1:n
        combos(:, i) = grid{i}(:);
    end
end
