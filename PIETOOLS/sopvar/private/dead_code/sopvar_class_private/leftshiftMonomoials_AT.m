function [Zs_new, ns_new, CL] = leftshiftMonomoials_AT(Lbeta, Zs_mult, var_Zs)
%The function computes CL ans Zs_new such that
% Zs_mult^T Lbeta = Zs_new^T CL
%
% INPUTS:
% Lbeta   -- n times m quadpoly object
% Zs_mult -- 1xI cell array of monomial exponents as in quadpoly objects
%            (Zs_mult{1} \otimes ... \otimes Zs_mult{I}) -- n times 1
%            monomial vector
% var_Zs  -- 1xI cell array of varnames for Zs_mult monomials
%
%
% OUTPUTS:
% Zt_new -- cell array of monomial exponents (as in the quadpoly object)
%           this vector represents Z(s4| s2a, s3a)
%           Z(s4| s2a, s3a) has the size of q times 1
% nt_new -- cell array of varnames for monomials in Zt_new 
% CL     -- sparse matrix size of (q times m)
%
%


P_subs1 = Lbeta;
ns = Lbeta.ns;
nt = Lbeta.nt;
dim = Lbeta.dim;
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

% Then combine two set of monomials Zs_mult and P.Zs
Zs_full =  table2cell([cell2table(Zs_mult), cell2table(Zs_new)]);
ns_full =  table2cell([cell2table(var_Zs), cell2table(ns_to_renamed)]);

% Given Zs_full, ns_full, we sort the variable names as in combine function
P_subs1 = quadPoly(C_new, Zs_full, {[0]}, [1, dim(2)], ns_full, {'<empty_var>'}, 0)   ;
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

end

