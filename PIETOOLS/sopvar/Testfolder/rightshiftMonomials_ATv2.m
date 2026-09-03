function [Zt_new,nt_new, CR] = rightshiftMonomials_ATv2(Ralpha,Zs_mult, var_Zs, n1, n2)
% The function computes CR and Zs_new such that 
% (I_n2 otimes Ralpha) (I_n1 otimes Zs_mult) =  CR (I_n2 otimes Zs_new)
%
% INPUTS:
% Ralpha  -- n times m quadpoly object
% Zs_mult -- 1xI cell array of monomial exponents as in quadpoly objects
%            (Zs_mult{1} \otimes ... \otimes Zs_mult{I}) -- nz times 1
%            monomial vector
% var_Zs  -- 1xI cell array of varnames for Zs_mult monomials
%
%
% OUTPUTS:
% Zt_new -- cell array of monomial exponents (as in the quadpoly object)
%           this vector represents Z(s3a'| s3b, s1)
%           Z(s3a'| s3b, s1) has the size of q times 1
% nt_new -- cell array of varnames for monomials in Zt_new 
% CR     -- sparse matrix size of (n2*n times n1*q)
%
%
% This code computes (I_n2 otimes Ralpha) (I_n1 otimes Zs_mult) 
% Ralpha -- (n times m)  quadpoly 
% Zs_mult-- (nz times 1) cell array as in quadpoly
%
% [(I_n2 otimes Ralpha) (I_n1 otimes Zs_mult)]^T  = 
%  (I_n1 otimes Zs_mult)^T  (I_n2 otimes Ralpha)^T 
%
% using leftshiftMonomials_ATv2 we have
% (I_n1 otimes Zs_mult)^T  (I_n2 otimes Ralpha)^T  = (In1 otimes Zs_new)^T CL
%
% Then 
% (I_n2 otimes Ralpha) (I_n1 otimes Zs_mult) = CL^T (In1 otimes Zs_new) 

[Zt_new, nt_new, CL] = leftshiftMonomoials_ATv2(Ralpha', Zs_mult, var_Zs, n1, n2);
CR = CL';


end

