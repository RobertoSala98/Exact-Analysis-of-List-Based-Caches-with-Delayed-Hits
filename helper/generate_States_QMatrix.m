function [States, Q] = generate_States_QMatrix(n, m, replacement_policy, S, F, p_access, alpha, T, R, ProcessSharing, access_probabilities, access_fetching, sym_computation)
    if isequal(replacement_policy, 'const')
        States = generate_states(n, m, 'rr', S, F);
    else
        States = generate_states(n, m, replacement_policy, S, F);
    end
    Q = generate_QMatrix(n, m, replacement_policy, S, F, p_access, alpha, T, R, ProcessSharing, access_probabilities, access_fetching, States, sym_computation);
end