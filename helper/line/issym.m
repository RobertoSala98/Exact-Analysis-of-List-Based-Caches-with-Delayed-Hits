function b=issym(x)
% This function is taken from LINE: https://github.com/imperial-qore/line-solver
    b=strcmpi(class(x),'sym');
end
