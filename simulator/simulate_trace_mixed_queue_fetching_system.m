function [miss_ratios, delayed_hits_ratios, hit_ratios, mean_req_fetching, mean_fetch_time, mean_squared_fetch_time, mean_number_awaiting_completition, mean_arrived_jobs_during_fetch, mean_latency, utilization] = simulate_trace_mixed_queue_fetching_system(n, m, replacement_policy, S, F, station_type, access_probabilities, access_fetching, p_access, alpha, T, R, num_events, csv_name, seed, verbose)
    
    if nargin < 14
        store_csv = false;
    else
        store_csv = true;
    end
    
    if nargin < 15
        seed = "default";
    end
    rng(seed)

    if nargin < 16
        verbose = false;
    end
    
    h = numel(m);
    mt = sum(m);

    % Job arrivals
    init_num_arrival = mt^2;

    lambda_arrival = sum(p_access);
    inter_arrival_times = exprnd(1/lambda_arrival, num_events + init_num_arrival, 1);
    arrival_times = cumsum(inter_arrival_times);
    item_types = randsample(n, num_events + init_num_arrival, true, p_access/sum(p_access));
    
    % Initial condition
    act_time = 0.0;
    cache = cell(1, h);
    cache_state = zeros(n, h); 
    
    item_num_init = n; 
    for kk = 1:h
        cache{kk} = zeros(1, m(kk)); 
        for jj = 1:m(kk)
            if item_num_init > 0
                cache{kk}(jj) = item_num_init;
                cache_state(item_num_init, kk) = jj;
                item_num_init = item_num_init - 1;
            end
        end
    end
    
    fetching_times_val = inf(n, 1);
    fetching_state = zeros(n, 4);       % [entering_station, entering_phase, next_station, next_phase]

    PS_values = ones(1, S);
    for ss = 1:S
        if isequal(station_type(ss), "PS")
            PS_values(ss) = 0;
        end
    end
    
    miss_time = zeros(1, n);
    delayed_hits_time = cell(1, S); 
    for s_idx = 1:S
        delayed_hits_time{s_idx} = zeros(1, F(s_idx), n);
    end
    hit_time = zeros(h, n);
    %dht = zeros(n,S);
    tot_time = 0;
    
    act_arrivals = 0;
    actual_queueing_items = zeros(1, n);
    queueing_times = zeros(1, n);
    mean_fetch_time = zeros(1,n);
    mean_squared_fetch_time = zeros(1,n);
    arrived_jobs_during_fetch = zeros(1,n);
    starting_squared_fetch_time = zeros(1,n);
    fetching_events = zeros(1,n);
    total_requests = zeros(1,n);
    mean_number_awaiting_completition = zeros(1,n);

    mean_latency = 0.0;

    raw_waiting_times = cell(n,1);

    arrivals_at_cache = zeros(1,n);
    started_analysing = false;

    utilization = zeros(1,S);

    if verbose
        %f_waitbar = waitbar(0, 'Trace simulation...');
        fid = fopen('log.csv', 'w');
        c = onCleanup(@() fclose(fid));

        fprintf(fid,'t');
        fprintf(fid,',m1');
        fprintf(fid,',m2');
        fprintf(fid,',d1');
        fprintf(fid,',d2');
        fprintf(fid,',d3');
        fprintf(fid,',d4');
        fprintf(fid,',Event');
        fprintf(fid,'\n');
        
        % build format once: 1 timestamp + (M+N) values, then newline
        fmt = ['%.6f', repmat(',%.9g',1,mt+n), ',%s', '\n'];
    end

    %{
    miss_time_prev = miss_time;
    hit_time_prev = hit_time;
    dht_prev = dht;
    %}
    
    % Iterations
    while act_arrivals < num_events || ~started_analysing
        [next_fetching_time, fetching_event_item_idx] = min(fetching_times_val);

        if act_arrivals + 1 > length(arrival_times) && ~started_analysing
            inter_arrival_times(act_arrivals + 1) = exprnd(1/lambda_arrival, 1, 1);
            arrival_times(act_arrivals + 1) = arrival_times(act_arrivals) + inter_arrival_times(act_arrivals + 1);
            item_types(act_arrivals + 1) = randsample(n, 1, true, p_access/sum(p_access));
        end

        next_arrival_time = arrival_times(act_arrivals + 1);
    
        % Stats update
        if sum(arrivals_at_cache) >= init_num_arrival
        
            if ~started_analysing
                new_arrivals = num_events - (size(inter_arrival_times,1) - act_arrivals);
                inter_arrival_times = [inter_arrival_times; exprnd(1/lambda_arrival, new_arrivals, 1)];
                arrival_times = cumsum(inter_arrival_times);
                item_types = [item_types; randsample(n, new_arrivals, true, p_access/sum(p_access))];
                num_events = size(inter_arrival_times,1);
                started_analysing = true;
            end
        
            update_time = min(next_arrival_time, next_fetching_time) - act_time;
            tot_time = tot_time + update_time;
        
            % ==========================================================
            % Vectorized hit / delayed-hit / miss accounting
            % ==========================================================
        
            cache_mask = cache_state > 0;
            has_cache  = any(cache_mask, 2);
            [~, cache_level] = max(cache_mask, [], 2);
        
            % ---- CACHE HITS ------------------------------------------
            hit_idx = find(has_cache);
            if ~isempty(hit_idx)
                lin_idx = sub2ind(size(hit_time), cache_level(hit_idx), hit_idx);
                hit_time(lin_idx) = hit_time(lin_idx) + update_time;
            end
        
            % ---- NO CACHE -------------------------------------------
            no_cache_idx = find(~has_cache);
        
            if ~isempty(no_cache_idx)
                station = fetching_state(no_cache_idx,1);
                phase   = fetching_state(no_cache_idx,2);
        
                % ---- DELAYED HITS ------------------------------------
                delayed_mask = station > 0;
                delayed_idx  = no_cache_idx(delayed_mask);
        
                for k = 1:numel(delayed_idx)
                    ii = delayed_idx(k);
                    st = fetching_state(ii,1);
                    ph = fetching_state(ii,2);
                    delayed_hits_time{st}(1,ph,ii) = delayed_hits_time{st}(1,ph,ii) + update_time;
                end
        
                % ---- MISSES ------------------------------------------
                miss_idx = no_cache_idx(~delayed_mask);
                if ~isempty(miss_idx)
                    miss_time(miss_idx) = miss_time(miss_idx) + update_time;
                end
            end

            queueing_times = queueing_times + actual_queueing_items * update_time;
            mean_fetch_time = mean_fetch_time + (actual_queueing_items > 0) * update_time;
            mean_number_awaiting_completition = mean_number_awaiting_completition + (actual_queueing_items - (actual_queueing_items > 0)) * update_time;
        
            % ==========================================================
            % Vectorized total-time-per-item computation
            % ==========================================================
        
            % hit + miss contribution
            tot_time_ii = sum(hit_time, 1).' + miss_time;
        
            % delayed hits contribution
            for ss = 1:S
                tot_time_ii = tot_time_ii + squeeze(sum(sum(delayed_hits_time{ss},1),2));
                if isequal(station_type(ss), "PS")
                    if PS_values(ss) > 0
                        utilization(ss) = utilization(ss) + update_time;
                    end
                end
            end
        end

        act_time = min(next_arrival_time, next_fetching_time);

        %if mod(act_arrivals, 1e5) == 0
        %    fprintf("%d \n", act_arrivals)
        %end

        %{
        if mod(act_arrivals, 1e4) == 0
            for s = 1:S
                dht(:,s) = sum(delayed_hits_time{s}(:,:,:));
            end
    
            miss_diff = abs(miss_time ./ act_time * p_access' - miss_time_prev ./ act_time * p_access') / (miss_time / act_time * p_access');
            hit_diff = abs(sum(hit_time,1) ./ act_time * p_access' - sum(hit_time_prev,1) ./ act_time * p_access') / (sum(hit_time,1) ./ act_time * p_access');
            dh_diff = abs(sum(dht,2)' ./ act_time *p_access' - sum(dht_prev,2)' ./ act_time *p_access') / (sum(dht,2)' ./ act_time *p_access');
        
            fprintf("%d, %.8f -> %.4f %.4f %.4f\n", act_arrivals, max(max(miss_diff, hit_diff), dh_diff), miss_diff, dh_diff, hit_diff)

            miss_time_prev = miss_time;
            hit_time_prev = hit_time;
            dht_prev = dht;
        end
        %}
        
        %{
        if mod(act_arrivals, 1e3) == 0
            clc
            act_arrivals
            miss_time * p_access' / act_time
            sum(hit_time * p_access') / act_time
        end
        %}
    
        if next_arrival_time < next_fetching_time           % new arrival
            act_arrivals = act_arrivals + 1;
    
            %if verbose
            %    waitbar(act_arrivals/num_events, f_waitbar);
            %end
    
            item = item_types(act_arrivals);
            item_cache_level = find(cache_state(item,:) > 0, 1);
    
            if isempty(item_cache_level)        % Item not in cache
                if fetching_state(item, 1) == 0         % Not currently fetching

                    if rand(1) <= access_fetching(item)      % Can enter in fetching phase
    
                        entering_station = pickFromDiscreteDistribution(R(1,:,item));
                        entering_PH_phase = 0;
                        if entering_station <= S
                            entering_PH_phase = pickPHorCache(alpha{entering_station,item});
                        end
        
                        if isequal(entering_station, S+1) || isequal(entering_PH_phase, 0)    % Directly enters cache (level 1)

                            if actual_queueing_items(item) ~= 0
                                fprintf("ERROR HERE!")
                            end

                            total_requests(item) = total_requests(item) + 1;

                            if rand(1) <= access_probabilities(item, 1)

                                if isequal(replacement_policy, 'rr')
                                    exit_item_idx_in_cache = randi(m(1));
                                    item_that_exits = cache{1}(exit_item_idx_in_cache);
                                    cache{1}(exit_item_idx_in_cache) = item;
            
                                    cache_state(item_that_exits, 1) = 0;
                                    cache_state(item, 1) = exit_item_idx_in_cache;
                                    
                                    if ~started_analysing
                                        arrivals_at_cache(item) = arrivals_at_cache(item) + 1;
                                    end
                                end
                                if isequal(replacement_policy, 'fifo') || isequal(replacement_policy, 'lru')
                                    item_that_exits = cache{1}(m(1));
                                    cache{1}(2:m(1)) = cache{1}(1:m(1)-1);
                                    cache{1}(1) = item;
            
                                    cache_state(item_that_exits, 1) = 0;
                                    
                                    if m(1) > 1
                                        cache_state(cache{1}(2:m(1)), 1) = (2:m(1))';
                                    end
                                    cache_state(item, 1) = 1;
        
                                    if ~started_analysing
                                        arrivals_at_cache(item) = arrivals_at_cache(item) + 1;
                                    end
                                end
        
                                if verbose
                                    %fprintf("Time %.6f: Item %d: cache level 1, Item %d: outside the cache.\n", act_time, item, item_that_exits)
                                    event = sprintf("Item %d enters in cache", item);
                                    row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                                    fprintf(fid, fmt, row);
                                end

                            end
        
                        else        % Enters fetching phase
                            
                            if actual_queueing_items(item) == 0
                                if item == 1
                                    arr_cache = cell2mat(cache);
                                    arr_qn = actual_queueing_items(:)';
                                    queued_ts = [];
                                end

                                actual_queueing_items(item) = 1;
                                fetching_events(item) = fetching_events(item) + 1;
                                starting_squared_fetch_time(item) = act_time;

                                raw_waiting_times{item} = [act_time];
                            else
                                fprintf("ERROR HERE 2")
                            end

                            [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, T{entering_station, item});
                            fetching_state(item, :) = [entering_station, entering_PH_phase, entering_station, next_fetching_state];
        
                            if isequal(station_type(entering_station), "PS")
                                PS_values(entering_station) = PS_values(entering_station) + 1;
    
                                active_items_in_station_mask = (fetching_state(:, 1) == entering_station) & (fetching_times_val ~= inf) & (1:n)' ~= item;
                                if any(active_items_in_station_mask)
                                    items_to_update = find(active_items_in_station_mask);
                                    residual_time = fetching_times_val(items_to_update) - act_time;
                                    factor = PS_values(entering_station) / (PS_values(entering_station) - 1);
                                    residual_time_modified = residual_time * factor;
                                    fetching_times_val(items_to_update) = residual_time_modified + act_time;
                                end
                            end
                            fetching_times_val(item) = act_time + sojourn_time * PS_values(entering_station);
        
                            if verbose
                                %fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from outside.\n", act_time, item, entering_station, entering_PH_phase)
                                event = sprintf("Item %d enters fetching station %d", item, entering_station);
                                row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                                fprintf(fid, fmt, row);
                            end
                        end

                    end

                else        % Item is already fetching, do nothing.

                    if actual_queueing_items(item) == 0
                        fprintf("ERROR HERE 3")
                    end

                    actual_queueing_items(item) = actual_queueing_items(item) + 1;
                    if item == 1
                        queued_ts = [queued_ts, act_time];
                    end

                    raw_waiting_times{item}(end+1) = act_time;

                    if verbose
                        %fprintf("Time %.6f: Item %d: already in fetching, station %d, phase %d.\n", act_time, item, fetching_state(item, 1), fetching_state(item, 2))
                        event = sprintf("Item %d is queued in the fetching system", item);
                        row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                        fprintf(fid, fmt, row);
                    end
                end
            else        % Item is already in cache

                if actual_queueing_items(item) ~= 0
                    fprintf("ERROR HERE 4")
                end

                total_requests(item) = total_requests(item) + 1;

                if item_cache_level < h         % Item can jump to the next level
                    target_cache_level = item_cache_level + 1;
                    original_pos_in_current_cache = cache_state(item, item_cache_level);

                    if rand(1) <= access_probabilities(item, target_cache_level)
                    
                        if isequal(replacement_policy, 'rr')
                            exit_item_idx_in_target_cache = randi(m(target_cache_level));
                            item_that_exits = cache{target_cache_level}(exit_item_idx_in_target_cache);
                            
                            cache{target_cache_level}(exit_item_idx_in_target_cache) = item;
                            cache_state(item, target_cache_level) = exit_item_idx_in_target_cache;
                            cache_state(item, item_cache_level) = 0; % Item no longer in its original cache level
        
                            cache{item_cache_level}(original_pos_in_current_cache) = item_that_exits;
                            cache_state(item_that_exits, item_cache_level) = original_pos_in_current_cache;
                            cache_state(item_that_exits, target_cache_level) = 0;
                        end
                        if isequal(replacement_policy, 'fifo')
                            item_that_exits = cache{target_cache_level}(m(target_cache_level));
                            
                            cache{target_cache_level}(2:m(target_cache_level)) = cache{target_cache_level}(1:m(target_cache_level)-1);
                            cache{target_cache_level}(1) = item;
                            cache_state(item, target_cache_level) = 1;
                            cache_state(item, item_cache_level) = 0;
                            
                            if m(target_cache_level) > 1
                                cache_state(cache{target_cache_level}(2:m(target_cache_level)), target_cache_level) = (2:m(target_cache_level))';
                            end
        
                            cache{item_cache_level}(original_pos_in_current_cache) = item_that_exits;
                            cache_state(item_that_exits, item_cache_level) = original_pos_in_current_cache;
                            cache_state(item_that_exits, target_cache_level) = 0;
                        end
                        if isequal(replacement_policy, 'lru')
                            item_that_exits = cache{target_cache_level}(m(target_cache_level));
                            
                            cache{target_cache_level}(2:m(target_cache_level)) = cache{target_cache_level}(1:m(target_cache_level)-1);
                            cache{target_cache_level}(1) = item;
                            cache_state(item, target_cache_level) = 1;
                            cache_state(item, item_cache_level) = 0;
    
                            if m(target_cache_level) > 1
                                cache_state(cache{target_cache_level}(2:m(target_cache_level)), target_cache_level) = (2:m(target_cache_level))';
                            end
    
                            cache{item_cache_level} = [item_that_exits, cache{item_cache_level}(1:original_pos_in_current_cache-1), cache{item_cache_level}(original_pos_in_current_cache+1:end)];
                            cache_state(item_that_exits, item_cache_level) = 1;
                            cache_state(item_that_exits, target_cache_level) = 0;
    
                            cache_state(cache{item_cache_level}(2:m(item_cache_level)), item_cache_level) = 2:m(item_cache_level);
                        end
                        if verbose
                            %fprintf("Time %.6f: Item %d: cache level %d, Item %d: cache level %d.\n", act_time, item, target_cache_level, item_that_exits, target_cache_level-1)
                            event = sprintf("Arrival item %d", item);
                            row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                            fprintf(fid, fmt, row);
                        end

                    end

                else        % Item is already in the last cache level (h)

                    if isequal(replacement_policy, 'lru')
                        original_pos_in_current_cache = cache_state(item, item_cache_level);
                        cache{h} = [item, cache{h}(1:original_pos_in_current_cache-1), cache{h}(original_pos_in_current_cache+1:end)];
                        
                        cache_state(cache{h}(:), h) = 1:m(h);
                    end

                    if verbose
                        %fprintf("Time %.6f: Item %d: already at cache level %d.\n", act_time, item, h)
                        event = sprintf("Arrival item %d", item);
                        row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                        fprintf(fid, fmt, row);
                    end
                end
            end
        else        % Process Fetching Completion Event

            item_id_completed = fetching_event_item_idx;
            fetching_times_val(item_id_completed) = inf;

            if actual_queueing_items(item_id_completed) == 0
                fprintf("ERROR HERE 5")
            end

            entering_fetching_state = fetching_state(item_id_completed, 4);
            previous_station = fetching_state(item_id_completed, 3);

            entering_PH_phase = 0;

            if entering_fetching_state == 0     % Routing
                %if verbose
                    %fprintf("Time %.6f: Picking a new station for Item %d. ", act_time, item_id_completed)
                %end
    
                next_station = pickFromDiscreteDistribution(R(previous_station+1,:,item_id_completed));
    
                if next_station <= S
                    entering_PH_phase = pickPHorCache(alpha{next_station, item_id_completed});
                    if verbose
                        %fprintf("Fetching station %d, phase %d.\n", next_station, entering_PH_phase)
                        event = sprintf("Item %d enters fetching station %d", item_id_completed, next_station);
                        row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                        fprintf(fid, fmt, row);
                    end
                %else
                    %if verbose
                        %fprintf("Cache.\n")
                    %end
                end
            else
                next_station = previous_station;
                entering_PH_phase = entering_fetching_state;
            end
    
            if isequal(next_station, S+1) || isequal(entering_PH_phase, 0)      % Enter in cache

                if actual_queueing_items(item_id_completed) == 0
                    fprintf("ERROR! \n")
                end

                if item_id_completed == 1
                    prev_cache = cell2mat(cache);
                    prev_qn = actual_queueing_items(:).';
                end

                if length(raw_waiting_times{item_id_completed}) ~= actual_queueing_items(item_id_completed)
                    disp("ERROR in computing latency")
                end
                
                mean_latency = mean_latency + sum(act_time - raw_waiting_times{item_id_completed});
                total_requests(item_id_completed) = total_requests(item_id_completed) + length(raw_waiting_times{item_id_completed});
                
                raw_waiting_times{item_id_completed} = [];

                arrived_jobs_during_fetch(item_id_completed) = arrived_jobs_during_fetch(item_id_completed) + actual_queueing_items(item_id_completed) - 1;
                actual_queueing_items(item_id_completed) = 0;
                mean_squared_fetch_time(item_id_completed) = mean_squared_fetch_time(item_id_completed) + (act_time - starting_squared_fetch_time(item_id_completed))^2;

                if rand(1) <= access_probabilities(item_id_completed, 1)

                    if isequal(replacement_policy, 'rr')
                        exit_item_idx_in_cache = randi(m(1));
                        item_that_exits = cache{1}(exit_item_idx_in_cache);
                        cache{1}(exit_item_idx_in_cache) = item_id_completed;
        
                        cache_state(item_that_exits, 1) = 0;
                        cache_state(item_id_completed, 1) = exit_item_idx_in_cache;
                        
                        if ~started_analysing
                            arrivals_at_cache(item_id_completed) = arrivals_at_cache(item_id_completed) + 1;
                        end
                    end
                    if isequal(replacement_policy, 'fifo') || isequal(replacement_policy, 'lru')
                        item_that_exits = cache{1}(m(1));
                        cache{1}(2:m(1)) = cache{1}(1:m(1)-1);
                        cache{1}(1) = item_id_completed;
        
                        cache_state(item_that_exits, 1) = 0;
                        if m(1) > 1
                            cache_state(cache{1}(2:m(1)), 1) = (2:m(1))';
                        end
                        cache_state(item_id_completed, 1) = 1;
    
                        if ~started_analysing
                            arrivals_at_cache(item_id_completed) = arrivals_at_cache(item_id_completed) + 1;
                        end
                    end
        
                    if verbose
                        %fprintf("Time %.6f: Item %d: cache level 1 from fetching station %d at PH phase %d, Item %d: outside the cache.\n", act_time, item_id_completed, previous_station, fetching_state(item_id_completed, 2), item_that_exits)
                        event = sprintf("Item %d exits from the fecthing system", item_id_completed);
                        row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                        fprintf(fid, fmt, row);
                    end
             
                end

                if isequal(station_type(previous_station), "PS")
                    PS_values(previous_station) = PS_values(previous_station) - 1;

                    active_items_in_station_mask = (fetching_state(:, 1) == previous_station) & (fetching_times_val ~= inf) & (1:n)' ~= item_id_completed;

                    if any(active_items_in_station_mask)
                        items_to_update = find(active_items_in_station_mask);
                        residual_time = fetching_times_val(items_to_update) - act_time;
                        factor = PS_values(previous_station) / (PS_values(previous_station) + 1);
                        residual_time_modified = residual_time * factor;
                        fetching_times_val(items_to_update) = residual_time_modified + act_time;
                    end
                end

                fetching_state(item_id_completed, :) = [0, 0, 0, 0];        % Reset fetching state

                if item_id_completed == 1 & store_csv
                    logFetchRow(csv_name, starting_squared_fetch_time(item_id_completed), arr_cache, arr_qn, prev_cache, prev_qn, cell2mat(cache), actual_queueing_items(:).', queued_ts, act_time);
                    queued_ts = [];
                end
    
            else        % Enter in a new fetching station or phase
                if next_station ~= previous_station
                    [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, T{next_station, item_id_completed});
                    fetching_state(item_id_completed, :) = [next_station, entering_PH_phase, next_station, next_fetching_state];
    
                    if isequal(station_type(previous_station), "PS")
                        PS_values(previous_station) = PS_values(previous_station) - 1;
                        active_items_in_prev_station_mask = (fetching_state(:, 1) == previous_station) & (fetching_times_val ~= inf) & (1:n)' ~= item_id_completed;
                        if any(active_items_in_prev_station_mask)
                            items_to_update = find(active_items_in_prev_station_mask);
                            residual_time = fetching_times_val(items_to_update) - act_time;
                            factor = PS_values(previous_station) / (PS_values(previous_station) + 1);
                            residual_time_modified = residual_time * factor;
                            fetching_times_val(items_to_update) = residual_time_modified + act_time;
                        end
                    end

                    if isequal(station_type(next_station), "PS")
                        PS_values(next_station) = PS_values(next_station) + 1;
                        active_items_in_next_station_mask = (fetching_state(:, 1) == next_station) & (fetching_times_val ~= inf) & (1:n)' ~= item_id_completed;
                        if any(active_items_in_next_station_mask)
                            items_to_update = find(active_items_in_next_station_mask);
                            residual_time = fetching_times_val(items_to_update) - act_time;
                            factor = PS_values(next_station) / (PS_values(next_station) - 1);
                            residual_time_modified = residual_time * factor;
                            fetching_times_val(items_to_update) = residual_time_modified + act_time;
                        end
                    end
                    fetching_times_val(item_id_completed) = act_time + sojourn_time * PS_values(next_station);
    
                    if verbose
                        %fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from fetching station %d.\n", act_time, item_id_completed, next_station, entering_PH_phase, previous_station)
                        event = sprintf("Item %d enters fetching station %d", item_id_completed, next_station);
                        row = [act_time, cell2mat(cache), actual_queueing_items, event];    % force row shape
                        fprintf(fid, fmt, row);
                    end
    
                else
                    [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, T{next_station, item_id_completed});
                    fetching_times_val(item_id_completed) = act_time + sojourn_time*PS_values(next_station);
    
                    %if verbose
                        %fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from phase %d.\n", act_time, item_id_completed, next_station, entering_PH_phase, fetching_state(item_id_completed, 2))
                    %end
                    fetching_state(item_id_completed, :) = [next_station, entering_PH_phase, next_station, next_fetching_state];
                end
            end
        end
    end
    
    %if verbose && exist('f_waitbar', 'var') && ishandle(f_waitbar)
    %    delete(f_waitbar);
    %end

    % Analysis
    miss_ratios = miss_time / tot_time;
    delayed_hits_ratios = cell(1,S);
    for s = 1:S
        delayed_hits_ratios{s} = delayed_hits_time{s} / tot_time;
    end
    hit_ratios = hit_time / tot_time;

    % (delayed_hits_ratios{1}(:,:,1) + delayed_hits_ratios{2}(:,:,1)) / miss_ratios(1,1) * (1 - sum(hit_ratios(:,1)))
    mean_req_fetching = queueing_times / tot_time;

    mean_fetch_time = mean_fetch_time ./ fetching_events;
    mean_squared_fetch_time = mean_squared_fetch_time ./ fetching_events;
    
    mean_number_awaiting_completition = mean_number_awaiting_completition / tot_time;

    mean_arrived_jobs_during_fetch = arrived_jobs_during_fetch ./ fetching_events;
    
    mean_latency = mean_latency / sum(total_requests);

    utilization = utilization / tot_time;
end

function logFetchRow(csvFile, t_miss, cache_arr, qn_arr, cache_end, qn_end, cache_post, qn_post, queued_ts, t_fetch_end)
%LOGFETCHROW Append one "fetch" record to a CSV log.
%   One row per fetch; fixed-length vectors stored as JSON strings too
%   (robust, easy to parse), queued_ts is variable length.

    arguments
        csvFile (1,1) string
        t_miss (1,1) double
        cache_arr {mustBeNumeric}
        qn_arr    {mustBeNumeric}
        cache_end {mustBeNumeric}
        qn_end    {mustBeNumeric}
        cache_post {mustBeNumeric}
        qn_post    {mustBeNumeric}
        queued_ts {mustBeNumeric}
        t_fetch_end (1,1) double
    end

    % Serialize vectors as JSON strings (single CSV cell each)
    cache_arr_json = string(jsonencode(cache_arr));
    qn_arr_json    = string(jsonencode(qn_arr));
    cache_end_json = string(jsonencode(cache_end));
    qn_end_json    = string(jsonencode(qn_end));
    cache_post_json = string(jsonencode(cache_post));
    qn_post_json    = string(jsonencode(qn_post));
    queued_ts_json = string(jsonencode(queued_ts));

    % Build one-row table
    T = table( ...
        t_miss, cache_arr_json, qn_arr_json, cache_end_json, qn_end_json, cache_post_json, qn_post_json, t_fetch_end, queued_ts_json,...
        'VariableNames', ["t_miss","cache_arr","qn_arr","cache_end","qn_end","cache_post","qn_post","t_fetch_end","queued_ts"]);

    % Write/append CSV with header on first write
    if isfile(csvFile)
        writetable(T, csvFile, "WriteMode","append", ...
            "WriteVariableNames", false, "QuoteStrings","all");
    else
        writetable(T, csvFile, "QuoteStrings","all");
    end
end
