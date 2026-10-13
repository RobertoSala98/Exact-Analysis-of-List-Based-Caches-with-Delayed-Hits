function [alpha, T, R] = generateData(n, S, F, min_scale_factor, max_scale_factor, seed)

    if nargin > 5
        rng(seed)
    else
        rng('shuffle');
    end

    if nargin < 5
        min_scale_factor = 1000;
        max_scale_factor = 10;
    end

    alpha = [];
    T = [];

    for s = 1:S
        if nargin > 5
            alpha = [alpha, generateAlphaVector(n, F(s), seed)];
            T = [T, generaPHMatrix(n, F(s), min_scale_factor, max_scale_factor, seed)];
        else
            alpha = [alpha, generateAlphaVector(n, F(s))];
            T = [T, generaPHMatrix(n, F(s), min_scale_factor, max_scale_factor)];
        end
    end
    
    if nargin > 5
        R = generateProbabilityMatrix(n, S, seed);
    else
        R = generateProbabilityMatrix(n, S);
    end

end