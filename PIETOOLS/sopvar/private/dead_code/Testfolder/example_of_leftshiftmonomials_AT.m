clc; clear;

  
%%%%% TEST leftshiftMonomoials_AT %%%%%%
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2'});
dim = [2,2]; nMons = [5, 5]; maxdeg = [2, 2]; density = 0.1;

% create a random quadpoly object 
Lbeta = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density); % size [20,3]
% create the monomial vector of the size [20]
Zs_mult    = quadPoly.randquadPoly([1, 1], [dim(1), 1], var_Zs, {'t2b'}, maxdeg, density);
Zs_mult = Zs_mult.Zs;
%The function computes CL ans Zs_new such that
% Zs_mult^T Lbeta = Zs_new^T CL
[Zs_new, ns_new, CL] = leftshiftMonomoials_AT(Lbeta, Zs_mult, var_Zs);

P_subs1 = quadPoly(CL, Zs_new, {[0]},  [1, dim(2)], ns_new, {'<empty_var>'}, 0)   ;


%%%%% TEST rightshiftMonomials_AT %%%%%%
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2'});
dim = [2,2]; nMons = [5, 5]; maxdeg = [2, 2]; density = 0.1;
Ralpha = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density);
Zs_mult    = quadPoly.randquadPoly([1, 1], [1, dim(2)], {'t2b'}, var_Zs,  maxdeg, density);
Zs_mult = Zs_mult.Zt;

%The function computes CR ans Zs_new such that
% Ralpha*Zs_mult = CR Zs_new
[Zt_new,nt_new, CR] = rightshiftMonomials_AT(Ralpha, Zs_mult, var_Zs);

P_subs1 = quadPoly(CR, {[0]}, Zt_new , [dim(1), 1],  {'<empty_var>'},nt_new, 0)   ;



%%%%% TEST leftshiftMonomoials_ATv2 %%%%%%
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2'});
nz = 3;
n = 4;
m = 2;
n1 = 4; n2 = 3;

dim = [n ,m]; nMons = [5, 5]; maxdeg = [2, 2]; density = 0.1;

% create a random quadpoly object 
Lbeta = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density);
% create the monomial vector of the size [20]
Zs_mult = quadPoly.randquadPoly([1, 1], [nz, 1], var_Zs, {'t2b'}, maxdeg, density);
Zs_mult = Zs_mult.Zs;


CLbeta = Lbeta.C;
ZLLbeta = Lbeta.Zs; nsLbeta = Lbeta.ns;
ZRLbeta = Lbeta.Zt; ntLbeta = Lbeta.nt;

%%%%% NAIVE APPROACH FOR TESTING %%%%%
L_beta_kron = quadPoly(kron(eye(n2), CLbeta), ZLLbeta, ZRLbeta, dim*n2, nsLbeta, ntLbeta);
% P_subs1 = L_beta_kron;
% ns = P_subs1.ns;
% nt = P_subs1.nt;
% dim = P_subs1.dim;
% 
% % first we move all monomials for Ralpha to the right side
% ns_RENAMED = cellfun(@(name) char(strcat("<subs_>", name)), ns, 'UniformOutput', false);
% 
% 
% for var_id = 1:length(ns)
%     name_OLD = ns{var_id};
%     name_NEW = ns_RENAMED{var_id};
%     P_subs1 = var_swap(P_subs1, name_OLD, name_NEW);%first 7 symbols include <subs_>
% end    
% 
% nt_LBETA = P_subs1.nt; 
% C_LBETA  = P_subs1.C;
% 
% % remove <subs_> from varnames
% nt_to_renamed = cellfun(@(name) char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive')),...
%     nt_LBETA, 'UniformOutput', false); % names to rename
% 
% % Then combine two set of monomials Zs_mult and P.Zs
% Zt_full =  P_subs1.Zt;
% nt_full =  nt_to_renamed;
% 
% % Given Zs_full, ns_full, we sort the variable names as in combine function
% % P_subs1 = quadPoly(C_new, {[0]}, Zt_full,  [dim(1), dim(2)], {'<empty_var>'}, nt_full, 0);
% 
% P_subs1 = quadPoly(C_LBETA, Zs_mult, Zt_full,  [n1, dim(2)], var_Zs, nt_full, 0);


dim = [n ,m]; 
P_subs_temp = move_all_monomials_to_right(L_beta_kron);

C_temp = P_subs_temp.C; Zt_temp = P_subs_temp.Zt; nt_temp = P_subs_temp.nt; 

P_subs1 = quadPoly(C_temp, Zs_mult, Zt_temp, [n1, dim(2)*n2], var_Zs, nt_temp);

P_subs1 = move_all_monomials_to_left(P_subs1);

P_subs1 = combine(P_subs1);

P_subs1 = clean_mon(P_subs1);

P_subs1 = combine(P_subs1);





[Zs_new, ns_new, CL] = leftshiftMonomoials_ATv2(Lbeta, Zs_mult, var_Zs, n1, n2);

P_subs2 = quadPoly(CL, Zs_new, {[0]},  [n1, dim(2)*n2], ns_new, {}, 0)   ;


if P_subs2 == P_subs1
    fprintf('TEST PASSED: leftshiftMonomoials_ATv2 \n')
end
% P_subs2 = quadPoly(CL1, Zs_new1, {[0]},  [1, dim(2)], ns_new1, {'<empty_var>'}, 0)   ;



%%%%% TEST rightshiftMonomoials_ATv2 %%%%%%
var_s = {'s2a1', 's2a2'};
var_t = {'s2a3', 's2a4'};
var_Zs = union(union(var_s, var_t), {'s3a1', 's3a2'});
nz = 3;
n = 2;
m = 4;
n1 = 4; n2 = 3;

dim = [n ,m]; nMons = [5, 5]; maxdeg = [2, 2]; density = 0.1;

Ralpha = quadPoly.randquadPoly(dim, nMons, var_s, var_t, maxdeg, density);
Zs_mult    = quadPoly.randquadPoly([1, 1], [1, nz], {'t2b'}, var_Zs,  maxdeg, density);
Zs_mult = Zs_mult.Zt;

% The function computes CR ans Zs_new such that
% (I_n2 otimes Ralpha) (I_n1 otimes Zs_mult) =  CR (I_n2 otimes Zs_new) 
[Zt_new,nt_new, CR] = rightshiftMonomials_ATv2(Ralpha, Zs_mult, var_Zs, n1, n2);
 
 
P_subs2 = quadPoly(CR, {[0]}, Zt_new , [dim(1)*n2, n1],  {},nt_new, 0)   ;




CRalpha = Ralpha.C;
ZLRalpha = Ralpha.Zs; nsRalpha = Ralpha.ns;
ZRRalpha = Ralpha.Zt; ntRalpha = Ralpha.nt;

R_alpha_kron = quadPoly(kron(eye(n2), CRalpha), ZLRalpha, ZRRalpha, ...
                        dim*n2, nsRalpha, ntRalpha);
 
P_subs_temp = move_all_monomials_to_left(R_alpha_kron);

C_temp = P_subs_temp.C; Zs_temp = P_subs_temp.Zs; ns_temp = P_subs_temp.ns; 

P_subs1 = quadPoly(C_temp, Zs_temp, Zs_mult, [dim(1)*n2, n1], ns_temp, var_Zs);

P_subs1 = move_all_monomials_to_right(P_subs1);

P_subs1 = combine(P_subs1);

P_subs1 = clean_mon(P_subs1);

P_subs1 = combine(P_subs1);



if P_subs2 == P_subs1
    fprintf('TEST PASSED: rightshiftMonomials_ATv2 \n')
end

function [quadpoly_left] = move_all_monomials_to_left(quadpoly)
    P_subs1 = quadpoly;
    ns = P_subs1.ns;
    nt = P_subs1.nt;
    dim = P_subs1.dim;
    
    % ns_RENAMED -- dummy varnames used in var_swap
    nt_RENAMED = cellfun(@(name) char(strcat("<subs_>", name)), nt, 'UniformOutput', false);
        
    
    % first we move all monomials in Ralpha to the left side  using var_swap 
    for var_id = 1:length(nt)
        name_OLD = nt{var_id};
        name_NEW = nt_RENAMED{var_id};
        P_subs1 = var_swap(P_subs1, name_NEW, name_OLD);% swap with dummy varnames
    end    
    
    ns_new = P_subs1.ns; 
    C_new  = P_subs1.C;
    Zs_new = P_subs1.Zs;
    % remove <subs_> from varnames
    ns_to_renamed = cellfun(@(name) char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive')),...
        ns_new, 'UniformOutput', false); % names to rename
    
    % Then combine two set of monomials Zs_mult and P.Zs

    ns_full =  ns_to_renamed;
 
    
    quadpoly_left = quadPoly(C_new, Zs_new, {[0]},  dim, ns_full, {}, 0);
end

 



function [quadpoly_left] = move_all_monomials_to_right(quadpoly)
    P_subs1 = quadpoly;
    ns = P_subs1.ns;
    nt = P_subs1.nt;
    dim = P_subs1.dim;
    
    % first we move all monomials for Ralpha to the right side
    ns_RENAMED = cellfun(@(name) char(strcat("<subs_>", name)), ns, 'UniformOutput', false);
        
    
    for var_id = 1:length(ns)
        name_OLD = ns{var_id};
        name_NEW = ns_RENAMED{var_id};
        P_subs1 = var_swap(P_subs1, name_OLD, name_NEW);%first 7 symbols include <subs_>
    end    
    
    nt_new = P_subs1.nt; 
    C_new  = P_subs1.C;
    Zt_new = P_subs1.Zt;
    % remove <subs_> from varnames
    nt_to_renamed = cellfun(@(name) char(eraseBetween(string(name), "<", ">",'Boundaries','inclusive')),...
        nt_new, 'UniformOutput', false); % names to rename
    
    % Then combine two set of monomials Zs_mult and P.Zs

    nt_full =  nt_to_renamed;
 
    
    quadpoly_left = quadPoly(C_new, {[0]}, Zt_new, dim, {'<empty_var>'},   nt_full, 0);
end


function [quadpoly_new] = clean_mon(quadpoly)
Zs = quadpoly.Zs;
ns = quadpoly.ns; 
Zt = quadpoly.Zt;
nt = quadpoly.nt; 
C = quadpoly.C;
dim = quadpoly.dim;

if length(ns) > 1
    
    is_repeated_in_ns = cellfun(@isequal, ns(1:end-1), ns(2:end));
    while sum(is_repeated_in_ns) > 0
        idx_to_merge_monomials = find(is_repeated_in_ns);% find repeated var names
        idx_to_merge_monomials = idx_to_merge_monomials(1);
        ns{idx_to_merge_monomials+1} = []; % delete repeated var name
        Zs_new_idx=  Zs{idx_to_merge_monomials}; % choose repeated monoms
        Zs_new_idx_p1=  Zs{idx_to_merge_monomials + 1};
        Zs{idx_to_merge_monomials + 1} = []; % delete dublicate
        % Define new monomial degrees
        Zs_new_merged = kron(Zs_new_idx, ones(size(Zs_new_idx_p1))) ...
                        + kron(ones(size(Zs_new_idx)), Zs_new_idx_p1);
        Zs{idx_to_merge_monomials} = Zs_new_merged;
        Zs = Zs(~cellfun(@isempty, ns)); % remove empty cells
        ns = ns(~cellfun(@isempty, ns)); % remove empty cells
        % check if any repeated var names
        if length(ns) > 1
            is_repeated_in_ns = cellfun(@isequal, ns(1:end-1), ns(2:end));
        else
            is_repeated_in_ns = {0};
        end
    end

end

if length(nt) > 1

    % merge monomial degrees for repeated ns
    is_repeated_in_nt = cellfun(@isequal, nt(1:end-1), nt(2:end));
    while sum(is_repeated_in_nt) > 0
        idx_to_merge_monomials = find(is_repeated_in_nt);% find repeated var names
        idx_to_merge_monomials = idx_to_merge_monomials(1);
        nt{idx_to_merge_monomials+1} = []; % delete repeated var name
        Zt_new_idx=  Zt{idx_to_merge_monomials}; % choose repeated monoms
        Zs_new_idx_p1=  Zt{idx_to_merge_monomials + 1};
        Zt{idx_to_merge_monomials + 1} = []; % delete dublicate
        % Define new monomial degrees
        Zs_new_merged = kron(Zt_new_idx, ones(size(Zs_new_idx_p1))) ...
                        + kron(ones(size(Zt_new_idx)), Zs_new_idx_p1);
        Zt{idx_to_merge_monomials} = Zs_new_merged;
        Zt = Zt(~cellfun(@isempty, nt)); % remove empty cells
        nt = nt(~cellfun(@isempty, nt)); % remove empty cells
        % check if any repeated var names
        if length(nt) > 1
            is_repeated_in_nt = cellfun(@isequal, nt(1:end-1), nt(2:end));
        else
            is_repeated_in_nt = {0};
        end
    end
end

quadpoly_new = quadPoly(C, Zs, Zt,  dim, ns, nt, 0);


end


function [quadpoly_new] = clean_right_mon(quadpoly)

end
 