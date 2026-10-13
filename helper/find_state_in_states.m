function idx = find_state_in_states(new_state, States)
    nRows = size(States, 1);
    idx = [];  % default if no match is found
    
    for i = 1:nRows
        match = true;

        % Compare each cell's contents
        for j = 1:numel(new_state)
            % Check both size and contents
            if ~isequal(new_state{j}, States{i, j})
                match = false;
                break;
            end
        end

        if match
            idx = i;  % Found matching row
            return;   % Exit early on first match
        end
    end
end
