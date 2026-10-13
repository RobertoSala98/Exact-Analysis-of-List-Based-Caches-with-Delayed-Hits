function [sojourn_time, next_state] = samplePHStep(current_state, T)
    % SAMPLEPHSTEP: Simulates one jump of a PH process from the current state
    % INPUTS:
    %   current_state : current phase (1-based index)
    %   T             : sub-generator matrix (NxN)
    % OUTPUTS:
    %   sojourn_time  : time spent in current state
    %   next_state    : index of the next state, or 0 if absorbed
    
    N = size(T, 1);
    lambda = -T(current_state, current_state);      % exit rate
    sojourn_time = exprnd(1 / lambda);              % time before jump
    %sojourn_time = 1/lambda;
    %sojourn_time = unifrnd(0,1);

    % Probabilities to all possible outcomes: other states + absorption
    probs = T(current_state, :) / lambda;
    probs(current_state) = 0;                       % no self-transition
    p_absorb = 1 - sum(probs);                      % residual goes to absorption
    probs = [probs, p_absorb];                      % append absorption

    % Sample next state (1...N or N+1 == absorb)
    pick = find(mnrnd(1, probs));

    if pick <= N
        next_state = pick;
    else
        next_state = 0;  % absorbed
    end
end