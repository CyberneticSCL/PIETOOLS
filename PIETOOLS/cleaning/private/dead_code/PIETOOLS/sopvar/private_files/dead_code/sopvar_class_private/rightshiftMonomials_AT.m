function [Zt_new,nt_new, CR] = rightshiftMonomials_AT(Ralpha,Zs_mult, var_Zs)
%The function computes CR ans Zs_new such that
% Ralpha*Zs_mult = CR Zs_new
%
% INPUTS:
% Ralpha  -- n times m quadpoly object
% Zs_mult -- 1xI cell array of monomial exponents as in quadpoly objects
%            (Zs_mult{1} \otimes ... \otimes Zs_mult{I}) -- m times 1
%            monomial vector
% var_Zs  -- 1xI cell array of varnames for Zs_mult monomials
%
%
% OUTPUTS:
% Zt_new -- cell array of monomial exponents (as in the quadpoly object)
%           this vector represents Z(s3a'| s3b, s1)
%           Z(s3a'| s3b, s1) has the size of q times 1
% nt_new -- cell array of varnames for monomials in Zt_new 
% CR     -- sparse matrix size of (q times m)
%
%


P_subs1 = Ralpha;
ns = Ralpha.ns;
nt = Ralpha.nt;
dim = Ralpha.dim;
% ns_RENAMED -- dummy varnames used in var_swap
ns_RENAMED = cellfun(@(name) char(strcat("<subs_>", name)), ns, 'UniformOutput', false);
    
% first we move all left monomials in Ralpha to right using var_swap 
for var_id = 1:length(ns)
    name_OLD = ns{var_id};
    name_NEW = ns_RENAMED{var_id};
    P_subs1 = var_swap(P_subs1, name_OLD, name_NEW);% swap with dummy varnames
end    

ns_new = P_subs1.ns;
Zs_new = P_subs1.Zs;
nt_new = P_subs1.nt;
Zt_new = P_subs1.Zt;
C_new  = P_subs1.C;

% remove <subs_> from varnames
% remove all left dummy monimials 
nt_to_renamed = cellfun(@(name) char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive')),...
    nt_new, 'UniformOutput', false); % names to rename

% Then combine two set of monomials Zs_mult and P.Zs
Zt_full =  table2cell([cell2table(Zs_mult), cell2table(Zt_new)]);
nt_full =  table2cell([cell2table(var_Zs), cell2table(nt_to_renamed)]);

% Given Zs_full, ns_full, we sort the variable names as in combine function
P_subs1 = quadPoly(C_new, {[0]}, Zt_full,  [dim(1), 1], {'<empty_var>'}, nt_full, 0)   ;
P_subs1 = combine(P_subs1); % The combine function sorts varnames, we only need to remove dublicates

% finally Zt, nt may have repeated var names
Zt_new = P_subs1.Zt;  
nt_new = P_subs1.nt; 
CR = P_subs1.C; 

% merge monomial degrees for repeated ns
is_repeated_in_nt = cellfun(@isequal, nt_new(1:end-1), nt_new(2:end));
while sum(is_repeated_in_nt) > 0
    idx_to_merge_monomials = find(is_repeated_in_nt);% find repeated var names
    idx_to_merge_monomials = idx_to_merge_monomials(1);
    nt_new{idx_to_merge_monomials+1} = []; % delete repeated var name
    Zt_new_idx=  Zt_new{idx_to_merge_monomials}; % choose repeated monoms
    Zs_new_idx_p1=  Zt_new{idx_to_merge_monomials + 1};
    Zt_new{idx_to_merge_monomials + 1} = []; % delete dublicate
    % Define new monomial degrees
    Zs_new_merged = kron(Zt_new_idx, ones(size(Zs_new_idx_p1))) ...
                    + kron(ones(size(Zt_new_idx)), Zs_new_idx_p1);
    Zt_new{idx_to_merge_monomials} = Zs_new_merged;
    Zt_new = Zt_new(~cellfun(@isempty, nt_new)); % remove empty cells
    nt_new = nt_new(~cellfun(@isempty, nt_new)); % remove empty cells
    % check if any repeated var names
    if length(nt_new) > 1
        is_repeated_in_nt = cellfun(@isequal, nt_new(1:end-1), nt_new(2:end));
    else
        is_repeated_in_nt = {0};
    end
end


end

