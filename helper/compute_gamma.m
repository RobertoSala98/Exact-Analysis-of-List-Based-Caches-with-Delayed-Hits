function gamma = compute_gamma(n, h, p_access, access_probabilities, access_fetching)

    gamma = zeros(n,h);
    for ii = 1:n
        for hh = 1:h
            gamma(ii,hh) = p_access(ii)^hh * access_fetching(ii);
            for g = 1:hh
                gamma(ii,hh) = gamma(ii,hh) * access_probabilities(ii,g);
            end
        end
    end

end