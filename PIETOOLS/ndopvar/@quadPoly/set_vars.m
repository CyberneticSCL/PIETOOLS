function [A_out,s_idcs,t_idcs] = set_vars(A_in,ns_new,nt_new)
% A_OUT = SET_VARS(A_IN,NS_NEW,NT_NEW) takes a 'quadPoly' object and
% expresses in terms of variables NS_NEW and NT_NEW
%
% INPUTS
% - A_in:       'quadPoly' objects representing a polynomial;
% - ns_new:     1 x M cellstr specifying the new names of the left
%               variables. Must include A_in.ns as a subset;
% - nt_new:     1 x N cellstr specifying the new names of the right
%               variables. Must include A_in.nt as a subset;
%
% OUTPUTS
% - A_out:      'quadPoly' objects representing the same polynomial as the
%               input, but now expressed in terms of the specified
%               variables, so that
%                       A_out.ns = ns_new;    A_out.nt = nt_new;
% - s_idcs,     1 x 2 cell, specifying the indices relating the variables
%    t_idcs:    in A_in to those in A_out, so that
%                   A_out.ns(s_idcs{1}) = A_in.ns(s_idcs{2})
%                   A_out.nt(t_idcs{1}) = A_in.nt(t_idcs{2})
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - set_vars
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 06/19/2026: Initial coding

% Extract the variable names and monomials
vars1_A = A_in.ns;      N1 = numel(ns_new);
vars2_A = A_in.nt;      N2 = numel(nt_new);
Zs_old = A_in.Zs;       nZs_old = cellfun(@(a) numel(a),Zs_old);
Zt_old = A_in.Zt;       nZt_old = cellfun(@(a) numel(a),A_in.Zt);

% Make sure all variables appearing in A are covered
if any(~ismember(vars1_A,ns_new)) || any(~ismember(vars2_A,nt_new))
    error("New variable names must include variable names already appearing in the function.")
end

% Determine map from old to new variables
[~,s_idcs1,s_idcs2] = intersect(ns_new,vars1_A,'stable');
[~,t_idcs1,t_idcs2] = intersect(nt_new,vars2_A,'stable');
s_idcs = {s_idcs1,s_idcs2};
t_idcs = {t_idcs1,t_idcs2};

% Augment monomial basis with constant monomial for vars not appearing 
% in the object
Zs = num2cell(zeros(1,N1));      Zs(s_idcs1) = Zs_old(s_idcs2);
Zt = num2cell(zeros(1,N2));      Zt(t_idcs1) = Zt_old(t_idcs2);

% Account for possible reordering of the variables
is_reordered_s = any(sort(s_idcs2)~=s_idcs2);
is_reordered_t = any(sort(t_idcs2)~=t_idcs2);
if ~is_reordered_s && ~is_reordered_t
    % The order of the monomials has not changed
    % --> the order of coefficients need not be changed either
    A_out = A_in;
    A_out.ns = ns_new;      A_out.Zs = Zs;
    A_out.nt = nt_new;      A_out.Zt = Zt;
    return
end

% Decompose the coefficient matrix to reorder rows/columns
[rridcs,ccidcs,vals] = find(A_in.C);

if is_reordered_s
    % Split row indices into row index of matrix-valued polynomial and Zs
    % monomial index for each variable
    ridcs = ceil(rridcs/prod(nZs_old));
    Zsidcs_arr = [ridcs(:),zeros(numel(rridcs),numel(Zs_old))];
    for i=1:numel(Zs_old)
        % Remove the contribution of variables up to i from the monomials
        rridcs = rridcs(:) - (Zsidcs_arr(:,i)-1)*prod(nZs_old(i:end));
        % Determine indices of monomial along direction i
        Zsidcs_arr(:,i+1) = ceil(rridcs/prod(nZs_old(i+1:end)));
    end

    % Account for reordering of the variables
    nZs_new = nZs_old(s_idcs2);
    rridcs = 1+(Zsidcs_arr(:,[1;1+s_idcs2(:)])-1)*cumprod([nZs_new(:);1],'reverse');
end

if is_reordered_t
    % Split column indices into column index of matrix-valued polynomial
    % and Zt monomial index for each variable
    cidcs = ceil(ccidcs/prod(nZt_old));
    Ztidcs_arr = [cidcs(:),zeros(numel(ccidcs),numel(Zt_old))];
    for i=1:numel(Zt_old)
        % Remove the contribution of variables up to i from the monomials
        ccidcs = ccidcs(:) - (Ztidcs_arr(:,i)-1)*prod(nZt_old(i:end));
        % Determine indices of monomial along direction i
        Ztidcs_arr(:,i+1) = ceil(ccidcs/prod(nZt_old(i+1:end)));
    end

    % Account for reordering of the variables
    nZt_new = nZt_old(t_idcs2);
    ccidcs = 1+(Ztidcs_arr(:,[1;1+t_idcs2(:)])-1)*cumprod([nZt_new(:);1],'reverse');
end

% Declare the new polynomial
Cnew = sparse(rridcs,ccidcs,vals,size(A_in.C,1),size(A_in.C,2));
A_out = quadPoly(Cnew,Zs,Zt,size(A_in),ns_new,nt_new,true);

end