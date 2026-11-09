function P_tr = trace(P1)
%TRACE  Sum of diagonal elements.
%   TRACE(A) is the sum of the diagonal elements of A, which is
%   
% Input P1 quadPoly obj (n times n)
% Output sum_i obj(i, i)
 

if isa(P1, 'double')
    P_tr = trace(P1);
    return
end

if ~isa(P1, 'quadPoly')
    error('quadPoly:trace:badType', 'trace supports quadPoly and numeric inputs.');
end


sz = P1.dim;
if sz(1) ~= sz(2)
    error('quadPoly:trace:dimMismatch', 'quadPoly must be square.');
end
% extract monomials and coefficients.
Zt = P1.Zt;
Zs = P1.Zs;
nt = P1.nt;
ns = P1.ns;
C = P1.C;

n_Zt = 1;% length of right monomials
for idx = 1:size(Zt,2)
    n_Zt = n_Zt*size(Zt{idx}, 1);
end

n_Zs = 1;% length of left monomials
for idx = 1:size(Zs,2)
    n_Zs = n_Zs*size(Zs{idx}, 1);
end


% compute C matrix 
C_tr = C(1:n_Zs, 1:n_Zt);
for idx = 2:sz(1)
    left_idx  = n_Zs*(idx - 1)+1;
    right_idx = n_Zt*(idx - 1)+1;
    C_tr = C_tr + C(left_idx:(left_idx + n_Zs-1), right_idx:(right_idx+n_Zt-1));
end

% define trace quadpoly 
P_tr = quadPoly(C_tr, Zs, Zt, [1, 1], ns, nt);

end

