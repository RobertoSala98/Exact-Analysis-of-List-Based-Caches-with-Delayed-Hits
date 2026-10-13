function phase = pickPHorCache(alpha)
    % PICKPHORCACHE: Pick a PH phase (1...F) or 0 for direct cache entry
    % INPUT:
    %   alpha : 1xF vector of probabilities (should sum <= 1)
    % OUTPUT:
    %   phase : integer in 1...F for PH phase, or 0 for direct cache
    
    F = length(alpha);
    p_escape = min(0, 1 - sum(alpha));

    extended_probs = [alpha, p_escape];  % length F+1
    pick = find(mnrnd(1, extended_probs));  % 1 to F+1

    if pick <= F
        phase = pick;  % enter PH process
    else
        phase = 0;     % direct cache hit
    end
end