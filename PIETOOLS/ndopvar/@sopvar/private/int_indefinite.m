function [Cout,Zint] = int_indefinite(Z)
% Given Z a tensor basis of monomials, we perform the indefinite integration
% int(Z) with respect to each variable
Zint = cellfun(@(x) x+1, Z, 'UniformOutput',false);
C = cellfun(@(x) 1./(x+1), Z, 'UniformOutput',false);
Cout = 1;
for i=1:numel(C)
    Cout = kron(Cout,C{i});
end
end