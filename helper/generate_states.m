function States = generate_states(n, m, replacement_policy, S, F)

    States = [];

    % Generate cache state sequences
    if isequal(replacement_policy, 'rr')
        sequences = generate_disjoint_combinations(n, m);
    else
        sequences = generate_combinations(n, m);
    end

    num_seqs = size(sequences, 1);
    tempStates = cell(num_seqs, 1);  % Cell to hold results per iteration

    parfor i = 1:num_seqs
        localStates = [];

        if isequal(replacement_policy, 'rr')
            seq = sequences(i);
            seq = seq{1};
        else
            seq = sequences(i,:);
        end

        used = [seq{:}];
        complement = setdiff(1:n, used);

        assignments = assign_values_to_vector(n, complement, 0, S);
        num_assignments = size(assignments, 1);

        for aa = 1:num_assignments
            comb = generate_fetching_combinations(n, assignments(aa,:), F);
            num_comb = size(comb, 1);

            for cc = 1:num_comb
                localStates = [localStates; [seq assignments(aa,:) comb(cc,:)]];
            end
        end

        tempStates{i} = localStates;
    end

    % Combine all partial results
    States = vertcat(tempStates{:});
end
