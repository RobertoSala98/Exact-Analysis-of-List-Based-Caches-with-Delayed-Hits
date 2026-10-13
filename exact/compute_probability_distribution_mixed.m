function [pi, norm_const] = compute_probability_distribution_mixed(n, m, S, p_access, station_type, access_probabilities, access_fetching, States, sym_computation, constants)

    s = size(States,1);
    h = numel(m);

    if sym_computation
        pi = sym(zeros(s,1));
    else
        pi = zeros(s,1);
    end

    parfor i = 1:s
        state = States(i,:);

        if sym_computation
            pi_i = sym(1);
        else
            pi_i = 1;
        end

        for j = 1:h
            for k = 1:length(state{j})
                pi_i = pi_i * (p_access(state{j}(k)))^j * access_fetching(state{j}(k));

                for g = 1:j
                    pi_i = pi_i * access_probabilities(state{j}(k), g);
                end
            end
        end

        PS_values = ones(1, S);
        for ss = 1:S
            if isequal(station_type(ss), "PS")
                PS_values(ss) = sum(state{h+1} == ss);
            end
        end

        for k = 1:n
            PS_queue_item = state{h+1}(k);
            if PS_queue_item > 0
                pi_i = pi_i * constants{PS_queue_item}(1, state{h+2}(k), k);
            end
        end

        pi_i = pi_i * prod(factorial(PS_values));
        pi(i) = pi_i;
    end

    norm_const = sum(pi);
    pi = pi / norm_const;
end
