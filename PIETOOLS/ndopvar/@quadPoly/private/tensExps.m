function [E, Z3] = tensExps(Z1,Z2)
% Given two monomial column vectors Z1, Z2, this
% function calculates Z3 such that E*Z3 = Z1\otimes Z2.

Zfull = kron(Z1,Z2);

[Z3,~, E] = unique(Zfull);

E = sparse(1:numel(Zfull), E, 1, numel(Zfull), numel(Z3));
end