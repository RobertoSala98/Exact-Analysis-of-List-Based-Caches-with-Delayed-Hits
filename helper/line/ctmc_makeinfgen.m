function Q=ctmc_makeinfgen(Q)
% This function is taken from LINE: https://github.com/imperial-qore/line-solver
A=Q-diag(diag(Q)); Q=A-diag(sum(A,2));
end