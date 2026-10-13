function alpha_ = generateAlphaVector(n, F, seed)

    if nargin > 2
        rng(seed)
    else
        rng('shuffle');
    end

    alpha_ = {zeros(1,F,n)};

    for ii = 1:n

        v = rand(1, F);
        sommaTotale = sum(v);
        if sommaTotale > 1
            v = v / sommaTotale;
            v = v * rand();
        end
    
        alpha_{1}(:,:,ii) = v;
    
    end
end