function [miss_ratios, delayed_hits_ratios, hit_ratios, mean_latency, utilization] = simulate_trace_mixed_real_trace(trace_folder, file_names, last_idx, n, m, replacement_policy, S, F, station_type, access_probabilities, access_fetching, times, factors, R, recover, seed, verbose)
    
    % alpha non dipende dall'item
    % R non dipende dall'item
    % T sostituito da times e factors

    if nargin < 16
        seed = "default";
    end
    rng(seed)

    if nargin < 17
        verbose = false;
    end
    
    recovering = false;
    if recover == 0
        h = numel(m);
        mt = sum(m);
    
        init_num_arrival = mt;
        file_name = trace_folder + file_names + "_0.dat";
        data = readmatrix(file_name);
        num_events = size(data, 1);
        arrival_times = data(:,1);
        item_types = data(:,2) + 1;
    
        old_start = 0;
        trace_idx = 0;
    
        if trace_idx == last_idx
            last_file = true;
        else
            last_file = false;
        end
        
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

        mean_latency = 0.0;
        raw_waiting_times = cell(n,1);
        total_requests = zeros(1,n);
        utilization = zeros(1,S);
    
        arrivals_at_cache = zeros(1,n);
        started_analysing = false;
    
        %{
        miss_time_prev = miss_time;
        hit_time_prev = hit_time;
        dht_prev = dht;
        %}
    
        store_data = false;
    else
        recover_filename = trace_folder + "act_status_" + n + "_items_next_trace_" + recover;
        load(recover_filename);
        recovering = true;
    end
    
    % Iterations
    while ~last_file || act_arrivals < num_events
        [next_fetching_time, fetching_event_item_idx] = min(fetching_times_val);

        if act_arrivals - old_start + 1 > length(arrival_times)
            trace_idx = trace_idx + 1;
            file_name = trace_folder + file_names + "_" + trace_idx + ".dat";

            old_start = act_arrivals;

            data = readmatrix(file_name);
            num_events = act_arrivals + size(data, 1);
            arrival_times = data(:,1);
            item_types = data(:,2) + 1;

            if trace_idx == last_idx
                last_file = true;
            end

            store_data = true;
        end

        next_arrival_time = arrival_times(act_arrivals - old_start + 1);
    
        % Stats update
        if sum(arrivals_at_cache) >= init_num_arrival && ~recovering

            if ~started_analysing
                started_analysing = true;
            end

            update_time = min(next_arrival_time, next_fetching_time) - act_time;
            tot_time = tot_time + update_time;

            % ---- CACHE HITS ----------------------------------------------------
            % Find first positive cache level per row (vectorized)
            cache_mask = cache_state > 0;
            has_cache = any(cache_mask, 2);

            % First cache level per ii
            [~, cache_level] = max(cache_mask, [], 2);

            hit_idx = find(has_cache);
            lin_idx = sub2ind(size(hit_time), cache_level(hit_idx), hit_idx);
            hit_time(lin_idx) = hit_time(lin_idx) + update_time;

            % ---- NO CACHE ------------------------------------------------------
            no_cache_idx = find(~has_cache);

            if ~isempty(no_cache_idx)
                station = fetching_state(no_cache_idx, 1);
                %phase = fetching_state(no_cache_idx, 2);

                % ---- DELAYED HITS ----------------------------------------------
                delayed_mask = station > 0;
                delayed_idx  = no_cache_idx(delayed_mask);

                for k = 1:numel(delayed_idx)
                    ii = delayed_idx(k);
                    st = fetching_state(ii,1);
                    ph = fetching_state(ii,2);
                    delayed_hits_time{st}(1,ph,ii) = ...
                        delayed_hits_time{st}(1,ph,ii) + update_time;
                end

                % ---- MISSES ----------------------------------------------------
                miss_idx = no_cache_idx(~delayed_mask);
                miss_time(miss_idx) = miss_time(miss_idx) + update_time;
            end

            for ss = 1:S
                if isequal(station_type(ss), "PS")
                    if PS_values(ss) > 0
                        utilization(ss) = utilization(ss) + update_time;
                    end
                end
            end
        end

        act_time = min(next_arrival_time, next_fetching_time);
    
        if next_arrival_time < next_fetching_time           % new arrival

            if recovering
                fprintf("Starting events from file %d \n", trace_idx + 1);
                store_data = false;
                recovering = false;
            end

            if store_data
                store_filename = trace_folder + "act_status_" + n + "_items_next_trace_" + trace_idx;

                if trace_idx > 2
                    delete(trace_folder + "act_status_" + n + "_items_next_trace_" + (trace_idx - 2) + ".mat");
                end

                save(store_filename)
                disp(datetime('now'))
                fprintf("Starting events from file %d \n", trace_idx + 1);
                store_data = false;
            end

            act_arrivals = act_arrivals + 1;

            if mod(act_arrivals, 1000) == 0
                fprintf("Actual arrival: %d \n", act_arrivals)
            end
    
            item = item_types(act_arrivals - old_start);
            item_cache_level = find(cache_state(item,:) > 0, 1);
    
            if isempty(item_cache_level)        % Item not in cache
                if fetching_state(item, 1) == 0         % Not currently fetching

                    if rand(1) <= access_fetching(item)      % Can enter in fetching phase
    
                        entering_station = pickFromDiscreteDistribution(R(1,:));
                        entering_PH_phase = 0;
                        if entering_station <= S
                            entering_PH_phase = 1;
                        end
        
                        if isequal(entering_station, S+1) || isequal(entering_PH_phase, 0)    % Directly enters cache (level 1)

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
                                    fprintf("Time %.6f: Item %d: cache level 1, Item %d: outside the cache.\n", act_time, item, item_that_exits)
                                end

                            end
        
                        else        % Enters fetching phase

                            actual_queueing_items(item) = 1;
                            raw_waiting_times{item} = [act_time];

                            [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, -1/(times(item) * factors(entering_station)));
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
                                fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from outside.\n", act_time, item, entering_station, entering_PH_phase)
                            end
                        end

                    end

                else        % Item is already fetching, do nothing.
                    
                    actual_queueing_items(item) = actual_queueing_items(item) + 1;
                    raw_waiting_times{item}(end+1) = act_time;
                    
                    if verbose
                        fprintf("Time %.6f: Item %d: already in fetching, station %d, phase %d.\n", act_time, item, fetching_state(item, 1), fetching_state(item, 2))
                    end
                end
            else        % Item is already in cache

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
                            fprintf("Time %.6f: Item %d: cache level %d, Item %d: cache level %d.\n", act_time, item, target_cache_level, item_that_exits, target_cache_level-1)
                        end

                    end

                else        % Item is already in the last cache level (h)

                    if isequal(replacement_policy, 'lru')
                        original_pos_in_current_cache = cache_state(item, item_cache_level);
                        cache{h} = [item, cache{h}(1:original_pos_in_current_cache-1), cache{h}(original_pos_in_current_cache+1:end)];
                        
                        cache_state(cache{h}(:), h) = 1:m(h);
                    end

                    if verbose
                        fprintf("Time %.6f: Item %d: already at cache level %d.\n", act_time, item, h)
                    end
                end
            end
        else        % Process Fetching Completion Event
            item_id_completed = fetching_event_item_idx;
            fetching_times_val(item_id_completed) = inf;

            entering_fetching_state = fetching_state(item_id_completed, 4);
            previous_station = fetching_state(item_id_completed, 3);

            entering_PH_phase = 0;

            if entering_fetching_state == 0     % Routing
                if verbose
                    fprintf("Time %.6f: Picking a new station for Item %d. ", act_time, item_id_completed)
                end
    
                next_station = pickFromDiscreteDistribution(R(previous_station+1,:));
    
                if next_station <= S
                    entering_PH_phase = 1;
                    if verbose
                        fprintf("Fetching station %d, phase %d.\n", next_station, entering_PH_phase)
                    end
                else
                    if verbose
                        fprintf("Cache.\n")
                    end
                end
            else
                next_station = previous_station;
                entering_PH_phase = entering_fetching_state;
            end
    
            if isequal(next_station, S+1) || isequal(entering_PH_phase, 0)      % Enter in cache

                mean_latency = mean_latency + sum(act_time - raw_waiting_times{item_id_completed});
                total_requests(item_id_completed) = total_requests(item_id_completed) + length(raw_waiting_times{item_id_completed});
                raw_waiting_times{item_id_completed} = [];
                actual_queueing_items(item_id_completed) = 0;

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
                        fprintf("Time %.6f: Item %d: cache level 1 from fetching station %d at PH phase %d, Item %d: outside the cache.\n", act_time, item_id_completed, previous_station, fetching_state(item_id_completed, 2), item_that_exits)
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
    
            else        % Enter in a new fetching station or phase
                if next_station ~= previous_station
                    [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, -1/(times(item_id_completed) * factors(next_station)));
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
                        fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from fetching station %d.\n", act_time, item_id_completed, next_station, entering_PH_phase, previous_station)
                    end
    
                else
                    [sojourn_time, next_fetching_state] = samplePHStep(entering_PH_phase, -1/(times(item_id_completed) * factors(next_station)));
                    fetching_times_val(item_id_completed) = act_time + sojourn_time*PS_values(next_station);
    
                    if verbose
                        fprintf("Time %.6f: Item %d: fetching station %d, phase %d, from phase %d.\n", act_time, item_id_completed, next_station, entering_PH_phase, fetching_state(item_id_completed, 2))
                    end
                    fetching_state(item_id_completed, :) = [next_station, entering_PH_phase, next_station, next_fetching_state];
                end
            end
        end
    end

    % Analysis
    miss_ratios = miss_time / tot_time;
    delayed_hits_ratios = cell(1,S);
    for s = 1:S
        delayed_hits_ratios{s} = delayed_hits_time{s} / tot_time;
    end
    hit_ratios = hit_time / tot_time;

    mean_latency = mean_latency / sum(total_requests);

    utilization = utilization / tot_time;

    store_filename = trace_folder + "final_status_" + n + "_items";
    save(store_filename)
end