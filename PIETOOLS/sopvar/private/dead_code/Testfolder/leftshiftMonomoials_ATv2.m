function [Zs_new, ns_new, CL] = leftshiftMonomoials_ATv2(Lbeta, Zs_mult, var_Zs, n1, n2)
%The function computes CL ans Zs_new such that
% (I_n1 otimes Zs_mult)^T (I_n2 otimes Lbeta) 
%         = (I_n1 otimes Zs_new)^T CL
%
% INPUTS:
% Lbeta   -- n times m quadpoly object
% Zs_mult -- 1xI cell array of monomial exponents as in quadpoly objects
%            (Zs_mult{1} \otimes ... \otimes Zs_mult{I}) -- nz times 1
%            monomial vector 
% var_Zs  -- 1xI cell array of varnames for Zs_mult monomials
%
%
% OUTPUTS:
% Zt_new -- cell array of monomial exponents (as in the quadpoly object)
%           this vector represents Z(s4| s2a, s3a)
%           Z(s4| s2a, s3a) has the size of q times 1
% nt_new -- cell array of varnames for monomials in Zt_new 
% CL     -- sparse matrix size of (n1*q times m*n2)
%
%
% This code computes (I_n1 otimes Zs_mult)^T (I_n2 otimes Lbeta) 
% Lbeta  -- (n times m)  quadpoly 
% Zs_mult-- (nz times 1) cell array as in quadpoly
%
%
% 1st step: verify n2*n == n1*nz 
%            n1*nz == n2*n for correct MatVec/MatMat multiplications
% 2nd step: transform Lbeta = (I_n otimes ZL)^T C_Lbeta (I_m otimes ZR)^T 
%                  to (I_n otimes Z2)^T C_new 
% 3rd step: (I_n2 otimes Lbeta) = (I_{n2*n} otimes Z2)^T (I_n2 otimes C_new)
% 4th step: since n2*n = n1*nz we have 
%     (I_{n1} otimes Zs_mult)^T (I_{n2*n} otimes Z2)^T
%         =(I_{n1} otimes Zs_mult)^T (I_n1 otimes I_nz otimes Z2)^T 
%         =(I_n1 otimes (Zs_mult^T (I_nz otimes Z2)^T )
%         =(I_n1 otimes (Zs_mult^T otimes 1)(I_nz otimes Z2^T)) 
%         =(I_(n1) otimes (Zs_mult otimes Z2)^T)
% 5th step: compute Zs_new = Zs_mult otimes Z2 and remove dublicates
% The process updates (I_n2 otimes C_new). 
% Thus, CL does not have the kronecker structure. 

%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%
P_subs1 = Lbeta;
ns = Lbeta.ns;
nt = Lbeta.nt;
dim = Lbeta.dim;
length_Zmult = cellfun(@length, Zs_mult);
nz = prod(length_Zmult);

if (nz*n1 ~= n2*dim(1))
    error('sopvar:leftshiftMonomoials_ATv2:dimMismatch', 'Inner dimensions must match.');
end


%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%
% ns_RENAMED -- dummy varnames used in var_swap
nt_RENAMED = cellfun(@(name) char(strcat("<subs_>", name)), nt, 'UniformOutput', false);
    

% first we move all monomials in Lbeta to the left side  using var_swap 
for var_id = 1:length(nt)
    name_OLD = nt{var_id};
    name_NEW = nt_RENAMED{var_id};
    P_subs1 = var_swap(P_subs1, name_NEW, name_OLD);% swap with dummy varnames
end    

ns_new = P_subs1.ns;
Zs_new = P_subs1.Zs;
nt_new = P_subs1.nt;
Zt_new = P_subs1.Zt;
C_new  = P_subs1.C;

% remove <subs_> from varnames
% remove all right dummy monimials 
ns_to_renamed = cellfun(@(name) char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive')),...
    ns_new, 'UniformOutput', false); % names to rename

%%%%%%%%%%%%%%%% STEP 3-4 %%%%%%%%%%%%%%%%
% Then combine two set of monomials Zs_mult and P.Zs
Zs_full =  table2cell([cell2table(Zs_mult), cell2table(Zs_new)]);
ns_full =  table2cell([cell2table(var_Zs), cell2table(ns_to_renamed)]);


%%%%%%%%%%%%%%%% STEP 5 %%%%%%%%%%%%%%%%

% compute I_n2 otimes C_new
C_new = kron(eye(n2), C_new);
% Given Zs_full, ns_full, we sort the variable names as in combine function
P_subs1 = quadPoly(C_new, Zs_full, {[0]}, [n1, dim(2)*n2], ns_full, {'<empty_var>'}, 0)   ;
P_subs1 = combine(P_subs1); % The combine function sorts varnames, we only need to remove dublicates

% finally Zs, ns may have repeated var names
Zt_new = P_subs1.Zt;
Zs_new = P_subs1.Zs;
nt_new = P_subs1.nt;
ns_new = P_subs1.ns;
CL = P_subs1.C;
dim_new = P_subs1.dim;

% merge monomial degrees for repeated ns
is_repeated_in_ns = cellfun(@isequal, ns_new(1:end-1), ns_new(2:end));
while sum(is_repeated_in_ns) > 0
    idx_to_merge_monomials = find(is_repeated_in_ns);% find repeated var names
    idx_to_merge_monomials = idx_to_merge_monomials(1);
    ns_new{idx_to_merge_monomials+1} = []; % delete repeated var name
    Zs_new_idx=  Zs_new{idx_to_merge_monomials}; % choose repeated monoms
    Zs_new_idx_p1=  Zs_new{idx_to_merge_monomials + 1};
    Zs_new{idx_to_merge_monomials + 1} = []; % delete dublicate
    % Define new monomial degrees
    Zs_new_merged = kron(Zs_new_idx, ones(size(Zs_new_idx_p1))) ...
                    + kron(ones(size(Zs_new_idx)), Zs_new_idx_p1);
    Zs_new{idx_to_merge_monomials} = Zs_new_merged;
    Zs_new = Zs_new(~cellfun(@isempty, ns_new)); % remove empty cells
    ns_new = ns_new(~cellfun(@isempty, ns_new)); % remove empty cells
    % check if any repeated var names
    if length(ns_new) > 1
        is_repeated_in_ns = cellfun(@isequal, ns_new(1:end-1), ns_new(2:end));
    else
        is_repeated_in_ns = {0};
    end
end
% remove dublicated degrees of the monomials
P_subs1 = quadPoly(C_new, Zs_new, {[0]}, [n1, dim(2)*n2], ns_new, {'<empty_var>'}, 0);
P_subs1 = combine(P_subs1); 


Zs_new = P_subs1.Zs;
ns_new = P_subs1.ns;
CL = P_subs1.C;

end

