function T_ = generaPHMatrix(n, F, min_scale_factor, max_scale_factor, seed)

    if nargin > 4
        rng(seed)
    else
        rng('shuffle');
    end

    T_ = cell(1,1);
    T_{1} = zeros(F,F,n);

    for k = 1:n
        Q = zeros(F,F);
        for i = 1:F
            idx = setdiff(1:F, i);
            rates = 10.^(log10(1/min_scale_factor) + (log10(max_scale_factor/F) - log10(1/min_scale_factor)) * rand(1, F-1));  % Random values in [1/scale_factor, scale_factor]
            Q(i, idx) = rates;
            margine = 10.^(log10(1/min_scale_factor) + (log10(max_scale_factor/F) - log10(1/min_scale_factor)) * rand());     % Margin in [1/scale_factor, scale_factor]
            Q(i, i) = -sum(rates) - margine;
        end
        T_{1}(:,:,k) = Q;
    end
end
