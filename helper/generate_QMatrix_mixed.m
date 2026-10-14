function Q = generate_QMatrix_mixed(n, m, replacement_policy, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, States, sym_computation)
    s = size(States, 1);
    h = numel(m);

    if sym_computation
        Q = sym(zeros(s, s));
        for i = 1:s
            Q = compute_row(Q, i, States, h, n, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, m, replacement_policy);
        end
    else
        Q = zeros(s, s);
        Q_rows = cell(s, 1);

        %wb = waitbar(0, 'Building Q matrix...');
        %dq = parallel.pool.DataQueue;
        %afterEach(dq, @(~) update_wb(wb, s));  % closure captures wb and s by value

        parfor i = 1:s
            Q_rows{i} = compute_row_i(i, States, h, n, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, m, replacement_policy, s);
            %send(dq, i);
        end
        %close(wb);

        for i = 1:s
            Q(i, :) = Q_rows{i};
        end
    end
end

%% Local function — lives outside the main function's end
function update_wb(wb, s)
    persistent p;
    if isempty(p), p = 0; end
    p = p + 1;
    waitbar(p / s, wb, sprintf('Row %d / %d', p, s));
end

%% -----------------------------------------------------------------------
function row = compute_row_i(i, States, h, n, S, F, p_access, alpha, T, R, ...
    station_type, access_probabilities, access_fetching, m, replacement_policy, s)
% Computes the full i-th row of Q and returns it as a 1-x-s vector.

    row = zeros(1, s);
    act_state = States(i, :);

    if h > 0
        last_list = act_state{h};
    else
        last_list = [];
    end
    act_PS_queues    = act_state{h+1};
    act_fetching_state = act_state{h+2};

    % Count PS items per station
    PS_values = ones(1, S);
    for ss = 1:S
        if isequal(station_type(ss), "PS")
            PS_values(ss) = length(find(act_PS_queues == ss));
        end
    end

    % ----------------------------------------------------------------
    if isequal(replacement_policy, 'const')

        % Arrivals
        for j = 1:n
            if act_fetching_state(j) == 0
                for rr = 1:S
                    for jj_ = 1:F(rr)
                        new_state = act_state;
                        new_state{h+1}(j) = rr;
                        new_state{h+2}(j) = jj_;
                        idx_state = find_state_in_states(new_state, States);
                        row(idx_state) = row(idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1, jj_, j);
                    end
                end
            end
        end

        % Fetching state update
        for j = 1:n
            if act_PS_queues(j) > 0
                for jj_ = 1:F(act_PS_queues(j))
                    if ~isequal(jj_, act_fetching_state(j))
                        new_state = act_state;
                        new_state{h+2}(j) = jj_;
                        idx_state = find_state_in_states(new_state, States);
                        row(idx_state) = row(idx_state) + T{act_PS_queues(j)}(act_fetching_state(j), jj_, j) / PS_values(act_PS_queues(j));
                    end
                end

                for rr = 1:S
                    for jj_ = 1:F(rr)
                        if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * alpha{rr}(1, jj_, j);
                        end
                    end

                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;
                    idx_state = find_state_in_states(new_state, States);
                    row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * (1 - sum(alpha{rr}(1, :, j)));
                end

                new_state = act_state;
                new_state{h+1}(j) = 0;
                new_state{h+2}(j) = 0;
                idx_state = find_state_in_states(new_state, States);
                row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, S+1, j);
            end
        end

    % ----------------------------------------------------------------
    elseif isequal(replacement_policy, 'rr')

        % Arrivals
        for j = 1:n
            if ~any(last_list == j) && act_fetching_state(j) == 0

                act_list_idx = 0;
                for k = 1:h-1
                    if any(act_state{k} == j)
                        act_list_idx = k;
                        break;
                    end
                end

                if act_list_idx > 0     % hit in cache j = 1, ..., h-1
                    next_list_idx = act_list_idx + 1;
                    next_list = act_state{next_list_idx};
                    for jj_ = 1:length(next_list)
                        new_state = act_state;
                        jj = next_list(jj_);
                        new_state{act_list_idx}(new_state{act_list_idx} == j) = jj;
                        new_state{next_list_idx}(new_state{next_list_idx} == jj) = j;
                        new_state{act_list_idx} = sort(new_state{act_list_idx});
                        new_state{next_list_idx} = sort(new_state{next_list_idx});
                        idx_state = find_state_in_states(new_state, States);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, next_list_idx) * p_access(j) / length(next_list);
                        end
                    end
                else    % miss
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);

                            row(idx_state) = row(idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1, jj_, j);
                        end

                        next_list = act_state{1};
                        for jj_ = 1:length(next_list)
                            new_state = act_state;
                            jj = next_list(jj_);
                            new_state{1}(new_state{1} == jj) = j;
                            new_state{1} = sort(new_state{1});
                            idx_state = find_state_in_states(new_state, States);
                            if ~isempty(idx_state)
                                row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1, :, j))) / length(next_list);
                            end
                        end
                    end

                    next_list = act_state{1};
                    for jj_ = 1:length(next_list)
                        new_state = act_state;
                        jj = next_list(jj_);
                        new_state{1}(new_state{1} == jj) = j;
                        new_state{1} = sort(new_state{1});
                        idx_state = find_state_in_states(new_state, States);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j) / length(next_list);
                        end
                    end
                end
            end
        end

        % Fetching state update
        for j = 1:n
            if act_PS_queues(j) > 0
                for jj_ = 1:F(act_PS_queues(j))
                    if ~isequal(jj_, act_fetching_state(j))
                        new_state = act_state;
                        new_state{h+2}(j) = jj_;
                        idx_state = find_state_in_states(new_state, States);
                        row(idx_state) = row(idx_state) + T{act_PS_queues(j)}(act_fetching_state(j), jj_, j) / PS_values(act_PS_queues(j));
                    end
                end

                for rr = 1:S
                    for jj_ = 1:F(rr)
                        if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * alpha{rr}(1, jj_, j);
                        end
                    end

                    next_list = act_state{1};
                    for jj_ = 1:length(next_list)
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;
                        miss_state = new_state;
                        jj = next_list(jj_);
                        new_state{1}(new_state{1} == jj) = j;
                        new_state{1} = sort(new_state{1});
                        idx_state = find_state_in_states(new_state, States);
                        rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * (1 - sum(alpha{rr}(1, :, j))) / length(next_list);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                            idx_miss = find_state_in_states(miss_state, States);
                            row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                        else
                            idx_miss = find_state_in_states(miss_state, States);
                            row(idx_miss) = row(idx_miss) + rate;
                        end
                    end
                end

                next_list = act_state{1};
                if length(next_list) > 0
                    for jj_ = 1:length(next_list)
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;
                        miss_state = new_state;
                        jj = next_list(jj_);
                        new_state{1}(new_state{1} == jj) = j;
                        new_state{1} = sort(new_state{1});
                        idx_state = find_state_in_states(new_state, States);
                        rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, S+1, j) / length(next_list);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                            idx_miss = find_state_in_states(miss_state, States);
                            row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                        else
                            idx_miss = find_state_in_states(miss_state, States);
                            row(idx_miss) = row(idx_miss) + rate;
                        end
                    end
                else
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;
                    idx_miss = find_state_in_states(new_state, States);
                    row(idx_miss) = row(idx_miss) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, S+1, j);
                end
            end
        end

    % ----------------------------------------------------------------
    elseif isequal(replacement_policy, 'fifo')

        % Arrivals
        for j = 1:n
            if ~any(last_list == j) && act_fetching_state(j) == 0

                act_list_idx = 0;
                for k = 1:h-1
                    if any(act_state{k} == j)
                        act_list_idx = k;
                        break;
                    end
                end

                if act_list_idx > 0
                    next_list_idx = act_list_idx + 1;
                    next_list = act_state{next_list_idx};
                    new_state = act_state;
                    jj = next_list(m(next_list_idx));
                    new_state{next_list_idx}(2:end) = new_state{next_list_idx}(1:end-1);
                    new_state{next_list_idx}(1) = j;
                    new_state{act_list_idx}(new_state{act_list_idx} == j) = jj;
                    idx_state = find_state_in_states(new_state, States);
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, next_list_idx) * p_access(j);
                    end
                else
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1, jj_, j);
                        end

                        next_list = act_state{1};
                        new_state = act_state;
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;
                        idx_state = find_state_in_states(new_state, States);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1, :, j)));
                        end
                    end

                    next_list = act_state{1};
                    new_state = act_state;
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;
                    idx_state = find_state_in_states(new_state, States);
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j);
                    end
                end
            end
        end

        % Fetching state update
        for j = 1:n
            if act_PS_queues(j) > 0
                for jj_ = 1:F(act_PS_queues(j))
                    if ~isequal(jj_, act_fetching_state(j))
                        new_state = act_state;
                        new_state{h+2}(j) = jj_;
                        idx_state = find_state_in_states(new_state, States);
                        row(idx_state) = row(idx_state) + T{act_PS_queues(j)}(act_fetching_state(j), jj_, j) / PS_values(act_PS_queues(j));
                    end
                end

                for rr = 1:S
                    for jj_ = 1:F(rr)
                        if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * alpha{rr}(1, jj_, j);
                        end
                    end

                    next_list = act_state{1};
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;
                    miss_state = new_state;
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;
                    idx_state = find_state_in_states(new_state, States);
                    rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * (1 - sum(alpha{rr}(1, :, j)));
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                        idx_miss = find_state_in_states(miss_state, States);
                        row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                    else
                        idx_miss = find_state_in_states(miss_state, States);
                        row(idx_miss) = row(idx_miss) + rate;
                    end
                end

                next_list = act_state{1};
                new_state = act_state;
                new_state{h+1}(j) = 0;
                new_state{h+2}(j) = 0;
                miss_state = new_state;
                new_state{1}(2:end) = new_state{1}(1:end-1);
                new_state{1}(1) = j;
                idx_state = find_state_in_states(new_state, States);
                rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, S+1, j);
                if ~isempty(idx_state)
                    row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                    idx_miss = find_state_in_states(miss_state, States);
                    row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                else
                    idx_miss = find_state_in_states(miss_state, States);
                    row(idx_miss) = row(idx_miss) + rate;
                end
            end
        end

    % ----------------------------------------------------------------
    elseif isequal(replacement_policy, 'lru')

        % Arrivals
        for j = 1:n
            if ~any(last_list == j) && act_fetching_state(j) == 0

                act_list_idx = 0;
                for k = 1:h-1
                    if any(act_state{k} == j)
                        act_list_idx = k;
                        break;
                    end
                end

                if act_list_idx > 0
                    next_list_idx = act_list_idx + 1;
                    next_list = act_state{next_list_idx};
                    new_state = act_state;
                    jj = next_list(m(next_list_idx));
                    j_idx = find(new_state{act_list_idx} == j);
                    new_state{next_list_idx}(2:end) = new_state{next_list_idx}(1:end-1);
                    new_state{next_list_idx}(1) = j;
                    new_state{act_list_idx} = [jj, new_state{act_list_idx}(1:j_idx-1), new_state{act_list_idx}(j_idx+1:end)];
                    idx_state = find_state_in_states(new_state, States);
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, next_list_idx) * p_access(j);
                    end

                    act_list = act_state{act_list_idx};
                    new_state = act_state;
                    act_list(act_list == j) = [];
                    new_state{act_list_idx}(1) = j;
                    new_state{act_list_idx}(2:end) = act_list;
                    
                    idx_state = find_state_in_states(new_state, States);
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + (1 - access_probabilities(j, next_list_idx)) * p_access(j);
                    end
                else
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1, jj_, j);
                        end

                        next_list = act_state{1};
                        new_state = act_state;
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;
                        idx_state = find_state_in_states(new_state, States);
                        if ~isempty(idx_state)
                            row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1, :, j)));
                        end
                    end

                    next_list = act_state{1};
                    new_state = act_state;
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;
                    idx_state = find_state_in_states(new_state, States);
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j);
                    end
                end
            end

            if any(last_list == j)
                new_state = act_state;
                old_list = act_state{h};
                item_idx = find(old_list == j);
                new_state{h} = [j, old_list(1:item_idx-1), old_list(item_idx+1:end)];
                idx_state = find_state_in_states(new_state, States);
                row(idx_state) = row(idx_state) + p_access(j);
            end
        end

        % Fetching state update
        for j = 1:n
            if act_PS_queues(j) > 0
                for jj_ = 1:F(act_PS_queues(j))
                    if ~isequal(jj_, act_fetching_state(j))
                        new_state = act_state;
                        new_state{h+2}(j) = jj_;
                        idx_state = find_state_in_states(new_state, States);
                        row(idx_state) = row(idx_state) + T{act_PS_queues(j)}(act_fetching_state(j), jj_, j) / PS_values(act_PS_queues(j));
                    end
                end

                for rr = 1:S
                    for jj_ = 1:F(rr)
                        if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            idx_state = find_state_in_states(new_state, States);
                            row(idx_state) = row(idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * alpha{rr}(1, jj_, j);
                        end
                    end

                    next_list = act_state{1};
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;
                    miss_state = new_state;
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;
                    idx_state = find_state_in_states(new_state, States);
                    rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, rr, j) * (1 - sum(alpha{rr}(1, :, j)));
                    if ~isempty(idx_state)
                        row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                        idx_miss = find_state_in_states(miss_state, States);
                        row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                    else
                        idx_miss = find_state_in_states(miss_state, States);
                        row(idx_miss) = row(idx_miss) + rate;
                    end
                end

                next_list = act_state{1};
                new_state = act_state;
                new_state{h+1}(j) = 0;
                new_state{h+2}(j) = 0;
                miss_state = new_state;
                new_state{1}(2:end) = new_state{1}(1:end-1);
                new_state{1}(1) = j;
                idx_state = find_state_in_states(new_state, States);
                rate = (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1, S+1, j);
                if ~isempty(idx_state)
                    row(idx_state) = row(idx_state) + access_probabilities(j, 1) * rate;
                    idx_miss = find_state_in_states(miss_state, States);
                    row(idx_miss) = row(idx_miss) + (1 - access_probabilities(j, 1)) * rate;
                else
                    idx_miss = find_state_in_states(miss_state, States);
                    row(idx_miss) = row(idx_miss) + rate;
                end
            end
        end
    end

    % Diagonal element
    row(i) = -sum(row) + row(i);
end
