function constants = compute_constants(n, S, F, p_access, alpha, T, R, sym_computation, access_fetching)

    constants = [];
    alpha_prov = [];
    T_prov = [];

    if sym_computation
        for ss = 1:S
            constants = [constants, {sym(zeros(1,F(ss),n))}];
            alpha_prov = [alpha_prov, {sym(zeros(1,F(ss),1))}];
            T_prov = [T_prov, {sym(zeros(F(ss),F(ss),1))}];
            R_prov = sym(zeros(S+1, S+1, 1));
        end
    else
        for ss = 1:S
            constants = [constants, {zeros(1,F(ss),n)}];
            alpha_prov = [alpha_prov, {zeros(1,F(ss),1)}];
            T_prov = [T_prov, {zeros(F(ss),F(ss),1)}];
            R_prov = zeros(S+1, S+1, 1);
        end
    end
    %{
    parfor item = 1:n

        for ss = 1:S
            alpha_prov{ss}(:,:,1) = alpha{ss}(:,:,item);
            T_prov{ss}(:,:,1) = T{ss}(:,:,item);
        end
        R_prov(:,:,1) = R(:,:,item);

        [States_prov, Q_prov] = generate_States_QMatrix(1, [], 'const', S, F, [p_access(item)], alpha_prov, T_prov, R_prov, 'false', [], access_fetching(item), sym_computation);
        %fprintf("The reduced system has %d states.\n", size(States_prov, 1))
        pi = ctmc_solve(Q_prov);

        counter = 1;
        for ss = 1:S
            for ff = 1:F(ss)
                counter = counter + 1;
                constants{ss}(1,ff,item) = pi(counter) / pi(1);
            end
        end
    end
    %}
    constants_local = cell(S,1);

    parfor item = 1:n
        alpha_prov = cell(S,1);
        T_prov     = cell(S,1);
    
        for ss = 1:S
            alpha_prov{ss} = alpha{ss}(:,:,item);
            T_prov{ss}     = T{ss}(:,:,item);
        end
        R_prov = R(:,:,item);
    
        [States_prov, Q_prov] = generate_States_QMatrix(1, [], 'const', S, F, ...
                                    [p_access(item)], alpha_prov, T_prov, R_prov, ...
                                    'false', [], access_fetching(item), sym_computation);
        pi = ctmc_solve(Q_prov);
    
        c_temp = cell(S,1);
        counter = 1;
        for ss = 1:S
            c_temp{ss} = zeros(1,F(ss));
            for ff = 1:F(ss)
                counter = counter + 1;
                c_temp{ss}(1,ff) = pi(counter)/pi(1);
            end
        end
        constants_local{item} = c_temp;
    end
    
    for item = 1:n
        for ss = 1:S
            constants{ss}(1,1:F(ss),item) = constants_local{item}{ss};
        end
    end

end