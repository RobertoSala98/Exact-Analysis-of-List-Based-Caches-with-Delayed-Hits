function [alpha_sim, T_sim] = compute_alpha_T_simulator(n, S, alpha, T)

    alpha_sim = cell(S, n);
    T_sim = cell(S, n);
    for s = 1:S
        for ii = 1:n
            alpha_sim{s, ii} = alpha{s}(:,:,ii);
            T_sim{s, ii} = T{s}(:,:,ii);
        end
    end

end