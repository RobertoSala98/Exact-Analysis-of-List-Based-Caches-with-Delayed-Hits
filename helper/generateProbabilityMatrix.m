function R = generateProbabilityMatrix(n, S, seed)

    if nargin > 2
        rng(seed)
    else
        rng('shuffle');
    end

    R = zeros(S+1, S+1, n);

    for ii = 1:n
        M = zeros(S+1, S+1);

        for i = 1:(S+1)
            row = rand(1, S+1);
            row = row / sum(row);
            M(i, :) = row;
        end

        R(:,:,ii) = M;
    end
end