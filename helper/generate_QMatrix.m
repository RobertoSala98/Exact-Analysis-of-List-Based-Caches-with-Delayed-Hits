function Q = generate_QMatrix(n, m, replacement_policy, S, F, p_access, alpha, T, R, ProcessSharing, access_probabilities, access_fetching, States, sym_computation)

    s = size(States,1);
    h = numel(m);

    if sym_computation
        Q = sym(zeros(s,s));
    else
        Q = zeros(s,s);
    end


    if isequal(replacement_policy, 'const')
        
        for i = 1:s
            act_state = States(i,:);
            
            act_PS_queues = act_state{h+1};
            act_fetching_state = act_state{h+2};

            if ProcessSharing
                PS_values = zeros(1, S);
                for ss = 1:S
                    PS_values(ss) = length(find(act_PS_queues == ss));    
                end
            else
                PS_values = ones(1, S);
            end

            % Arrivals
            for j=1:n
                
                if act_fetching_state(j) == 0    
                    for rr = 1:S
                        % Enter in a fetching state for the specific PS queue, based on the routing probability
                        for jj_ = 1:F(rr)
                            new_state = act_state;
                            new_state{h+1}(j) = rr;
                            new_state{h+2}(j) = jj_;
                            
                            idx_state = find_state_in_states(new_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1,jj_,j);
                        end
                    end
                end
            end

            % Fetching state update
            for j=1:n
                if act_PS_queues(j) > 0

                    % Could change state in the same PS queue
                    for jj_ = 1:F(act_PS_queues(j))
                        if ~isequal(jj_, act_fetching_state(j)) 
                            new_state = act_state;
                            new_state{h+2}(j) = jj_;
                            
                            idx_state = find_state_in_states(new_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + T{act_PS_queues(j)}(act_fetching_state(j),jj_,j) / PS_values(act_PS_queues(j));
                        end
                    end

                    % Could exit the actual PS queue and enter another
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * alpha{rr}(1, jj_, j);
                            end
                        end

                        % Or miss
                        next_list = act_state{1};
                        
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;

                        idx_state = find_state_in_states(new_state, States);
                        Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));
                    end

                    % Could exit the actual PS queue and directly miss
                    next_list = act_state{1};
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;
                    
                    idx_state = find_state_in_states(new_state, States);
                    Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);
                end
            end
            
            % Update the diagonal element
            Q(i,i) = -sum(Q(i,:)) + Q(i,i);

        end
    end


    if isequal(replacement_policy, 'rr')
        
        for i = 1:s
            act_state = States(i,:);
            if h > 0
                last_list = act_state{h};
            else
                last_list = [];
            end
            act_PS_queues = act_state{h+1};
            act_fetching_state = act_state{h+2};
            
            % Count PS items
            if ProcessSharing
                PS_values = zeros(1, S);
                for ss = 1:S
                    PS_values(ss) = length(find(act_PS_queues == ss));    
                end
            else
                PS_values = ones(1, S);
            end

            % Arrivals
            for j=1:n
                
                % If the item is not in the last list of the cache, nor in fetching
                if ~any(last_list == j) && act_fetching_state(j) == 0
                    
                    % Find where the item is
                    act_list_idx = 0;
                    for k=1:h-1
                        if any(act_state{k} == j)
                            act_list_idx = k;
                            break;
                        end
                    end
                    
                    % If the item is in any list of the cache, but not the last one
                    if act_list_idx > 0
                        next_list_idx = act_list_idx + 1;
                        next_list = act_state{next_list_idx};
                        
                        for jj_ = 1:length(next_list)
                            new_state = act_state;
                            jj = next_list(jj_);     % Item to be evicted from the next list

                            new_state{act_list_idx}(new_state{act_list_idx} == j) = jj;
                            new_state{next_list_idx}(new_state{next_list_idx} == jj) = j;

                            new_state{act_list_idx} = sort(new_state{act_list_idx});
                            new_state{next_list_idx} = sort(new_state{next_list_idx});

                            idx_state = find_state_in_states(new_state, States);
                            if ~isempty(idx_state)
                                Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, next_list_idx) * p_access(j) / length(next_list);
                            end
                        end
                    
                    % Otherwise enter in fetching state 
                    else

                        % Select a PS queue from the routing probabilities
                        for rr = 1:S
                            % Enter in a fetching state for the specific PS queue, based on the routing probability
                            for jj_ = 1:F(rr)
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                    
                                Q(i,idx_state) = Q(i,idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1,jj_,j);
                            end

                            % Or directly in cache
                            next_list = act_state{1};
                            
                            for jj_ = 1:length(next_list)
                                new_state = act_state;
                                jj = next_list(jj_);     % Item to be evicted from the next list
    
                                new_state{1}(new_state{1} == jj) = j;
                                new_state{1} = sort(new_state{1});
    
                                idx_state = find_state_in_states(new_state, States);
                                if ~isempty(idx_state)
                                    Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1,:,j))) / length(next_list);
                                end
                            end
                        end
                        
                        % Or directly enter in cache
                        next_list = act_state{1};

                        for jj_ = 1:length(next_list)
                            new_state = act_state;
                            jj = next_list(jj_);     % Item to be evicted from the next list

                            new_state{1}(new_state{1} == jj) = j;
                            new_state{1} = sort(new_state{1});

                            idx_state = find_state_in_states(new_state, States);

                            if ~isempty(idx_state)
                                Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j) / length(next_list);
                            end
                        end
                    end
                end
            end

            % Fetching state update
            for j=1:n
                if act_PS_queues(j) > 0

                    % Could change state in the same PS queue
                    for jj_ = 1:F(act_PS_queues(j))
                        if ~isequal(jj_, act_fetching_state(j)) 
                            new_state = act_state;
                            new_state{h+2}(j) = jj_;
                            
                            idx_state = find_state_in_states(new_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + T{act_PS_queues(j)}(act_fetching_state(j),jj_,j) / PS_values(act_PS_queues(j));
                        end
                    end

                    % Could exit the actual PS queue and enter another
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)))
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * alpha{rr}(1, jj_, j);
                            end
                        end

                        % Or enter in cache
                        next_list = act_state{1};
                        for jj_ = 1:length(next_list)
                            new_state = act_state;
                            new_state{h+1}(j) = 0;
                            new_state{h+2}(j) = 0;

                            miss_state = new_state;

                            jj = next_list(jj_);     % Item to be evicted from the next list

                            new_state{1}(new_state{1} == jj) = j;
                            new_state{1} = sort(new_state{1});

                            idx_state = find_state_in_states(new_state, States);
                            if ~isempty(idx_state)
                                Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j))) / length(next_list);

                                idx_state = find_state_in_states(miss_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j))) / length(next_list);
                            else
                                idx_state = find_state_in_states(miss_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j))) / length(next_list);
                            end
                        end
                    end

                    % Could exit the actual PS queue and enter directly in cache
                    next_list = act_state{1};
                    for jj_ = 1:length(next_list)
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;

                        miss_state = new_state;

                        jj = next_list(jj_);     % Item to be evicted from the next list

                        new_state{1}(new_state{1} == jj) = j;
                        new_state{1} = sort(new_state{1});

                        idx_state = find_state_in_states(new_state, States);
                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j) / length(next_list);
        
                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j) / length(next_list);
                        else
                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j) / length(next_list);
                        end
                    end
                end
            end
            
            % Update the diagonal element
            Q(i,i) = -sum(Q(i,:)) + Q(i,i);

        end
    end


    if isequal(replacement_policy, 'fifo')
        
        % Count PS items
        for i = 1:s
            act_state = States(i,:);
            if h > 0
                last_list = act_state{h};
            else
                last_list = [];
            end
            act_PS_queues = act_state{h+1};
            act_fetching_state = act_state{h+2};

            if ProcessSharing
                PS_values = zeros(1, S);
                for ss = 1:S
                    PS_values(ss) = length(find(act_PS_queues == ss));    
                end
            else
                PS_values = ones(1, S);
            end

            % Arrivals
            for j=1:n
                
                % If the item is not in the last list of the cache, nor in fetching
                if ~any(last_list == j) && act_fetching_state(j) == 0
                    
                    % Find where the item is
                    act_list_idx = 0;
                    for k=1:h-1
                        if any(act_state{k} == j)
                            act_list_idx = k;
                            break;
                        end
                    end
                    
                    % If the item is in any list of the cache, but not the last one
                    if act_list_idx > 0
                        next_list_idx = act_list_idx + 1;
                        next_list = act_state{next_list_idx};
                        
                        new_state = act_state;
                        jj = next_list(m(next_list_idx));     % Item to be evicted from the next list
                        
                        new_state{next_list_idx}(2:end) = new_state{next_list_idx}(1:end-1);
                        new_state{next_list_idx}(1) = j;
                        new_state{act_list_idx}(new_state{act_list_idx} == j) = jj;

                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, next_list_idx) * p_access(j);
                        end
                    
                    % Otherwise enter in fetching state 
                    else

                        % Select a PS queue from the routing probabilities
                        for rr = 1:S
                            % Enter in a fetching state for the specific PS queue, based on the routing probability
                            for jj_ = 1:F(rr)
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1,jj_,j);
                            end

                            % Or directly in cache
                            next_list = act_state{1};
                            
                            new_state = act_state;
                            jj = next_list(m(1));     % Item to be evicted from the next list
                            
                            new_state{1}(2:end) = new_state{1}(1:end-1);
                            new_state{1}(1) = j;
    
                            idx_state = find_state_in_states(new_state, States);

                            if ~isempty(idx_state)
                                Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1,:,j)));
                            end
                        end
                        
                        % Or directly enter in cache
                        next_list = act_state{1};
                        
                        new_state = act_state;
                        jj = next_list(m(1));     % Item to be evicted from the next list
                        
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;

                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j);
                        end
                    end
                end
            end

            % Fetching state update
            for j=1:n
                if act_PS_queues(j) > 0

                    % Could change state in the same PS queue
                    for jj_ = 1:F(act_PS_queues(j))
                        if ~isequal(jj_, act_fetching_state(j)) 
                            new_state = act_state;
                            new_state{h+2}(j) = jj_;
                            
                            idx_state = find_state_in_states(new_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + T{act_PS_queues(j)}(act_fetching_state(j),jj_,j) / PS_values(act_PS_queues(j));
                        end
                    end

                    % Could exit the actual PS queue and enter another
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)) )
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * alpha{rr}(1, jj_, j);
                            end
                        end

                        % Or enter in cache
                        next_list = act_state{1};
                        
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;

                        miss_state = new_state;

                        jj = next_list(m(1));     % Item to be evicted from the next list
                        
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;

                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));

                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));
                        else
                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));
                        end
                    end

                    % Could exit the actual PS queue and enter directly in cache
                    next_list = act_state{1};
                        
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;

                    miss_state = new_state;

                    jj = next_list(m(1));     % Item to be evicted from the next list
                    
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;

                    idx_state = find_state_in_states(new_state, States);

                    if ~isempty(idx_state)
                        Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);

                        idx_state = find_state_in_states(miss_state, States);
                        Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);
                    else
                        idx_state = find_state_in_states(miss_state, States);
                        Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);
                    end
                end
            end
            
            % Update the diagonal element
            Q(i,i) = -sum(Q(i,:)) + Q(i,i);

        end
    end


    if isequal(replacement_policy, 'lru')
        
        for i = 1:s

            if mod(i, 10) == 0
                fprintf('');
            end

            act_state = States(i,:);
            if h > 0
                last_list = act_state{h};
            else
                last_list = [];
            end
            act_PS_queues = act_state{h+1};
            act_fetching_state = act_state{h+2};

            if ProcessSharing
                PS_values = zeros(1, S);
                for ss = 1:S
                    PS_values(ss) = length(find(act_PS_queues == ss));    
                end
            else
                PS_values = ones(1, S);
            end

            % Arrivals
            for j=1:n
                
                % If the item is not in the last list of the cache, nor in fetching
                if ~any(last_list == j) && act_fetching_state(j) == 0
                    
                    % Find where the item is
                    act_list_idx = 0;
                    for k=1:h-1
                        if any(act_state{k} == j)
                            act_list_idx = k;
                            break;
                        end
                    end
                    
                    % If the item is in any list of the cache, but not the last one
                    if act_list_idx > 0
                        next_list_idx = act_list_idx + 1;
                        next_list = act_state{next_list_idx};
                        
                        new_state = act_state;
                        jj = next_list(m(next_list_idx));     % Item to be evicted from the next list
                        j_idx = find(new_state{act_list_idx} == j);
                        
                        new_state{next_list_idx}(2:end) = new_state{next_list_idx}(1:end-1);
                        new_state{next_list_idx}(1) = j;

                        new_state{act_list_idx} = [jj, new_state{act_list_idx}(1:j_idx-1), new_state{act_list_idx}(j_idx+1:end)];
                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, next_list_idx) * p_access(j);
                        end
                    
                    % Otherwise enter in fetching state 
                    else

                        % Select a PS queue from the routing probabilities
                        for rr = 1:S
                            % Enter in a fetching state for the specific PS queue, based on the routing probability
                            for jj_ = 1:F(rr)
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + access_fetching(j,1) * p_access(j) * R(1, rr, j) * alpha{rr}(1,jj_,j);
                            end

                            % Or directly in cache
                            next_list = act_state{1};
                            
                            new_state = act_state;
                            jj = next_list(m(1));     % Item to be evicted from the next list
                            
                            new_state{1}(2:end) = new_state{1}(1:end-1);
                            new_state{1}(1) = j;
    
                            idx_state = find_state_in_states(new_state, States);

                            if ~isempty(idx_state)
                                Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, rr, j) * (1 - sum(alpha{rr}(1,:,j)));
                            end
                        end
                        
                        % Or directly enter in cache
                        next_list = act_state{1};
                        
                        new_state = act_state;
                        jj = next_list(m(1));     % Item to be evicted from the next list
                        
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;

                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * access_fetching(j,1) * p_access(j) * R(1, S+1, j);
                        end
                    end
                end

                if any(last_list == j)
                    new_state = act_state;
                    
                    old_list = act_state{h};
                    item_idx = find(old_list == j);
                    new_state{h} = [j, old_list(1:item_idx-1), old_list(item_idx+1:end)];

                    idx_state = find_state_in_states(new_state, States);

                    Q(i,idx_state) = Q(i,idx_state) + p_access(j);
                end
            end

            % Fetching state update
            for j=1:n
                if act_PS_queues(j) > 0

                    % Could change state in the same PS queue
                    for jj_ = 1:F(act_PS_queues(j))
                        if ~isequal(jj_, act_fetching_state(j)) 
                            new_state = act_state;
                            new_state{h+2}(j) = jj_;
                            
                            idx_state = find_state_in_states(new_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + T{act_PS_queues(j)}(act_fetching_state(j),jj_,j) / PS_values(act_PS_queues(j));
                        end
                    end

                    % Could exit the actual PS queue and enter another
                    for rr = 1:S
                        for jj_ = 1:F(rr)
                            if ~(isequal(act_PS_queues(j), rr) && isequal(jj_, act_fetching_state(j)) )
                                new_state = act_state;
                                new_state{h+1}(j) = rr;
                                new_state{h+2}(j) = jj_;
                                
                                idx_state = find_state_in_states(new_state, States);
                                Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * alpha{rr}(1, jj_, j);
                            end
                        end

                        % Or enter in cache
                        next_list = act_state{1};
                        
                        new_state = act_state;
                        new_state{h+1}(j) = 0;
                        new_state{h+2}(j) = 0;

                        miss_state = new_state;

                        jj = next_list(m(1));     % Item to be evicted from the next list
                        
                        new_state{1}(2:end) = new_state{1}(1:end-1);
                        new_state{1}(1) = j;

                        idx_state = find_state_in_states(new_state, States);

                        if ~isempty(idx_state)
                            Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));

                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));
                        else
                            idx_state = find_state_in_states(miss_state, States);
                            Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,rr,j) * (1 - sum(alpha{rr}(1,:,j)));
                        end
                    end

                    % Could exit the actual PS queue and enter directly in cache
                    next_list = act_state{1};
                        
                    new_state = act_state;
                    new_state{h+1}(j) = 0;
                    new_state{h+2}(j) = 0;

                    miss_state = new_state;

                    jj = next_list(m(1));     % Item to be evicted from the next list
                    
                    new_state{1}(2:end) = new_state{1}(1:end-1);
                    new_state{1}(1) = j;

                    idx_state = find_state_in_states(new_state, States);

                    if ~isempty(idx_state)
                        Q(i,idx_state) = Q(i,idx_state) + access_probabilities(j, 1) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);

                        idx_state = find_state_in_states(miss_state, States);
                        Q(i,idx_state) = Q(i,idx_state) + (1 - access_probabilities(j, 1)) * (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);
                    else
                        idx_state = find_state_in_states(miss_state, States);
                        Q(i,idx_state) = Q(i,idx_state) + (-sum(T{act_PS_queues(j)}(act_fetching_state(j), :, j))) / PS_values(act_PS_queues(j)) * R(act_PS_queues(j)+1,S+1,j);
                    end
                end
            end
            
            % Update the diagonal element
            Q(i,i) = -sum(Q(i,:)) + Q(i,i);

        end
    end

end