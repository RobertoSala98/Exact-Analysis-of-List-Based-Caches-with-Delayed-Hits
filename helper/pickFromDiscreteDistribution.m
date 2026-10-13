function idx = pickFromDiscreteDistribution(P)
    C = cumsum(P);
    r = rand();
    idx = find(r <= C, 1, 'first');
end
