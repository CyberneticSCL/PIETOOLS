function A = uminus(A)
% Trivial. Just flip signs of coefficients. Returns -A given quadPoly A.
A.C = -A.C;
end