clc; clear;

  
 
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2', 's41', 's42', 's43'});
dim = [20,3]; nMons = [20, 10]; maxdeg = [5, 5]; density = 0.1;

% create a random quadpoly object 
Lbeta = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density); % size [20,3]
% create the monomial vector of the size [20]
Zs_mult    = quadPoly.randquadPoly([1, 1], [dim(1), 1], var_Zs, {'t2b'}, maxdeg, density);
Zs_mult = Zs_mult.Zs;
%The function computes CL ans Zs_new such that
% Zs_mult^T Lbeta = Zs_new^T CL
[Zs_new, ns_new, CL] = leftshiftMonomoials_AT(Lbeta, Zs_mult, var_Zs);



% for rightshiftMonomial
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2', 's41', 's42', 's43'});
dim = [2,20]; nMons = [20, 10]; maxdeg = [5, 5]; density = 0.1;
Ralpha = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density);
Zs_mult    = quadPoly.randquadPoly([1, 1], [1, dim(2)], {'t2b'}, var_Zs,  maxdeg, density);
Zs_mult = Zs_mult.Zt;

%The function computes CR ans Zs_new such that
% Ralpha*Zs_mult = CR Zs_new
[Zt_new,nt_new, CR] = rightshiftMonomials_AT(Ralpha, Zs_mult, var_Zs);


