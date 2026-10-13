function [mean_number_awaiting] = estimate_mean_number_awaiting_item(n, m, S, F, p_access, alpha, T, R, station_type, access_probabilities, access_fetching, item, phi_item)

    h = numel(m);

    States = [{0, 0}];
    for ss = 1:S
        for ff = 1:F(ss)
            States = [States; {ss, ff}];
        end
    end
    
    % System without item 
    n_ = n-1;
    p_access_ = p_access; p_access_(item) = [];
    alpha_ = alpha; T_ = T; R_ = R; for ss = 1:S alpha_{ss}(:,:,item) = []; T_{ss}(:,:,item) = []; end; R_(:,:,item) = [];
    access_probabilities_ = access_probabilities; access_fetching_ = access_fetching; access_probabilities_(item,:) = []; access_fetching_(item) = [];

    gamma_ = zeros(n_,h);
    for ii = 1:n_
        for hh = 1:h
            gamma_(ii,hh) = p_access_(ii)^hh * access_fetching_(ii);
            for g = 1:hh
                gamma_(ii,hh) = gamma_(ii,hh) * access_probabilities_(ii,g);
            end
        end
    end

    constants_ = compute_constants(n_, S, F, p_access_, alpha_, T_, R_, false, access_fetching_);
    
    eta_tmp_ =[];
    for ss = 1:S
        for kk = 1:n_
            eta_tmp_(ss,kk) = sum(constants_{ss}(:,:,kk));
        end
    end
    eta_ = eta_tmp_';
    
    [phi_, xi_, pi_vals_, pi_0_] = fetching_iteration_mixed(gamma_, eta_, m, station_type);
    
    % Computation
    n_j = sum(phi_, 1);

    for ss = 1:S
        if isequal(station_type(ss), "INF")
            n_j(ss) = 0;
        end
    end
    
    Q = zeros(length(States));
    
    for ss = 1:S
        for ff = 1:F(ss)
            idx_state = find_state_in_states({ss,ff}, States);
            Q(1, idx_state) = Q(1, idx_state) + access_fetching(item,1) * p_access(item) * R(1, ss, item) * alpha{ss}(1,ff,item);
    
            for ff2 = 1:F(ss)
                if ~isequal(ff, ff2)
                    new_state = find_state_in_states({ss,ff2}, States);
                    Q(idx_state, new_state) = Q(idx_state, new_state) + T{ss}(ff,ff2,item) / (1 + n_j(ss));
                end
            end     
    
            for ss2 = 1:S
                for ff2 = 1:F(ss2)
                    if ~(isequal(ss, ss2) && isequal(ff, ff2))
                        new_state = find_state_in_states({ss2,ff2}, States);
                        Q(idx_state, new_state) = Q(idx_state, new_state) + (-sum(T{ss}(ff, :, item))) / (1 + n_j(ss)) * R(ss+1,ss2,item) * alpha{ss2}(1, ff2, item);
                    end
                end

                Q(idx_state, 1) = Q(idx_state, 1) + (-sum(T{ss}(ff, :, item))) / (1 + n_j(ss)) * R(ss+1,ss2,item) * (1 - sum(alpha{ss2}(1,:,item)));
            end
    
            Q(idx_state, 1) = Q(idx_state, 1) + (-sum(T{ss}(ff, :, item))) / (1 + n_j(ss)) * R(ss+1,S+1,item);
        end
    end
    
    for i = 1:length(States) Q(i,i) = -sum(Q(i,:)) + Q(i,i); end

    D0 = Q(2:end, 2:end);
    pi_e = [];
    for ss = 1:S
        pi_e = [pi_e, R(1,ss,item)*alpha{ss}(:,:,item)];
    end
    first_moment = sum(pi_e * (-D0)^(-1));
    second_moment = sum(2 * pi_e * (-D0)^(-2));

    mean_number_awaiting = phi_item * p_access(item) * second_moment / (2 * first_moment);

end