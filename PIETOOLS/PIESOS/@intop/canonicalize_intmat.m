function [intmat,params,is_ord] = canonicalize_intmat(intmat,params,check_only)
% [INTMAT,PARAMS] = INTOP.CANONICALIZE_INTMAT(INTMAT,PARAMS) splits the
% rows of the intmat of an 'intop' object that are not ordered into sums of
% ordered rows, repeating the kernel of each row that is split, so that the
% declared functional keeps the same value, if the indicator function
% chi_{[x,y]}(z) is taken to be 1 if x<=z<=y, and 0 otherwise.
%
% [INTMAT] = INTOP.CANONICALIZE_INTMAT(INTMAT) only returns the new intmat.
%
% [~,~,IS_ORD] = INTOP.CANONICALIZE_INTMAT(INTMAT,[],TRUE) determines
% which rows are ordered, without splitting any rows.
%
% INPUTS
% - intmat:     1 x N cell with element i a q_i x 2*d_i array of integers
%               between 0 and d_i+1, specifying the intervals of
%               integration in the q_i terms of a functional along
%               direction i, as the field intmat of an 'intop' object. For
%               N=1, may also be a single array;
% - params:     (optional) q_1 x ... x q_N cell, specifying the kernel (or
%               the sparse coefficient matrix) of each term, as the first
%               input of the 'intop' constructor. Empty to only update
%               intmat;
% - check_only: (optional) logical, if true, no rows are split. Defaults to
%               false;
%
% OUTPUTS
% - intmat:     1 x N cell (or array) in which all rows are ordered. The
%               rows of the input that were ordered are retained unchanged,
%               in the same order, with the rows into which a row was split
%               in its place;
% - params:     cell of kernels matching the new intmat, in which the
%               kernel of a row that is split is repeated for each of the
%               rows obtained from it, and the kernels of rows for which
%               the region of integration is empty are removed;
% - is_ord:     1 x N cell with element i a q_i x 1 logical array, with
%               is_ord{i}(j) indicating whether row j of intmat{i} is
%               ordered (a logical array for N=1 and intmat an array);
%
% NOTES
% - This is a static method of the class 'intop': call it as
%   intop.canonicalize_intmat(...), not with an 'intop' object as input.
% - A row [l1,u1,...,ld,ud] defines the indicator functions
%   chi_{[x_lj,x_uj]}(x_j) of the variables x_1,...,x_d, with x_0=a and
%   x_{d+1}=b, and with lj=uj denoting a multiplier, x_j=x_lj. The row is
%   ordered if for every variable x_j for which both limits are variables
%   (1<=lj,uj<=d, lj~=uj, after eliminating multipliers), the inequality
%   x_lj<=x_uj is implied by the limits of the other variables. Only for
%   ordered rows does the integral with reversed limits (integrating x_j
%   from x_lj to x_uj with x_uj<x_lj, which 'apply_op' takes to be
%   negative) coincide with the integral over the region where x_lj<=x_uj
%   holds, as chi is a true indicator. E.g. [0,4,1,3,0,4] is not ordered
%   (x_2 in [x_1,x_3], but x_1 and x_3 are independent), whereas
%   [0,4,0,3,1,4] is (no variable has two variable limits).
% - The implication is determined recursively: a variable x_j with two
%   variable limits is added to the variables known to be ordered if
%   x_lj<=x_uj follows from the limits of the variables known to be
%   ordered. Using the limits of all other variables, including those of
%   variables that are not known to be ordered, would not be sound: e.g.
%   x_1,x_2 in [x_3,x_4] would then be ordered, although x_3<=x_4 is not
%   implied, and both would be negative if x_4<x_3.
% - A row that is not ordered is replaced by the cells of the ordered
%   regions: let S be the variables that are not ordered, together with the
%   variables defining their limits, and recursively the variables defining
%   the limits of these. For each ordering x_s1<=...<=x_sn of the variables
%   in S consistent with their limits, we declare a row with each x_sk
%   integrated from a up to x_s(k+1), and the last over [a,b], and all
%   other intervals unchanged. The region where the limits are reversed,
%   which has no cell, is dropped. For [0,4,1,3,0,4], we get the single row
%   [0,2,0,3,0,4], for x_1<=x_2<=x_3.
% - If no ordering of the variables in S is consistent with their limits,
%   e.g. for [3,4,1,3,0,4], where x_1>=x_3 and x_1<=x_2<=x_3, the region of
%   integration is empty (has measure zero) and the integral is zero. The
%   row is then removed, together with its kernel. If this leaves a
%   direction without rows, the functional is zero, and an error is thrown.
% - Only the integers in INTMAT are used, so that the cost is independent
%   of the number of decision variables, and is zero for rows in which no
%   variable has two variable limits. The kernels are only indexed, not
%   copied.
% - The 'intop' constructor and `isvalid` do not split rows, but reject an
%   intmat with rows that are not ordered, referring to this function.
%
% EXAMPLE
%   The row [0,4,1,3,0,4] has s2 between s1 and s3, which are
%   independent, so it is not ordered. It is replaced by the cell
%       s1<=s2<=s3:
%       intmat = {[0,4,1,3,0,4]};       params = {7};
%   [intmat,params] = intop.canonicalize_intmat(intmat,params)
%   returns intmat = {[0,2,0,3,0,4]} and params = {7}.
%
%   With s1,s2 in [s3,s4], we get two cells, s3<=s1<=s2<=s4 and
%   s3<=s2<=s1<=s4, and the kernel is repeated for both:
%       intmat = {[3,4,3,4,0,5,0,5]};
%   [intmat,params] = intop.canonicalize_intmat(intmat,{7})
%   returns intmat = {[0,2,0,4,0,1,0,5; 0,4,0,1,0,2,0,5]} and
%   params = {7;7}.
%
%   A row with an empty region is removed, together with its kernel:
%       intmat = {[0,4,0,3,1,4; 0,4,1,3,0,4; 3,4,1,3,0,4]};
%   [intmat,params] = intop.canonicalize_intmat(intmat,{'a';'b';'c'})
%   returns intmat = {[0,4,0,3,1,4; 0,2,0,3,0,4]} and
%   params = {'a';'b'}.
%
%   Only determine which rows are ordered:
%       intmat = {[0,4,0,3,1,4; 0,4,1,3,0,4]};
%   [~,~,is_ord] = intop.canonicalize_intmat(intmat,[],true)
%   returns is_ord = {[1;0]}.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - canonicalize_intmat
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
% DJ, 10/22/2026: Initial coding;

if nargin<2    
    params = [];            
end
if nargin<3
    check_only = false;
end
is_array = ~iscell(intmat);
if is_array
    intmat = {intmat};
end
N = numel(intmat);
is_ord = cell(1,N);
if ~isempty(params)
    params = reshape(params,[size(params,1:max(N,1)),1]);
end
for i=1:N
    [intmat_i,idcs,is_ord{i}] = canonicalize_dir(intmat{i},check_only);
    if check_only || all(is_ord{i})
        continue
    end
    if isempty(idcs)
        error("The integrals along direction %d are all over empty regions: the functional is zero.",i)
    end
    intmat{i} = intmat_i;
    if ~isempty(params)
        % N+1 subscripts: params remains q x 1 for N=1, also if q=1
        idx = repmat({':'},1,N+1);      idx{i} = idcs;
        params = params(idx{:});
    end
end
if is_array
    intmat = intmat{1};
    is_ord = is_ord{1};
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [intmat,idcs,is_ord] = canonicalize_dir(intmat,check_only)
% [INTMAT,IDCS,IS_ORD] = CANONICALIZE_DIR(INTMAT,CHECK_ONLY) splits the
% rows of the intmat along a single direction that are not ordered, see
% canonicalize_intmat. IDCS(k) is the row of the input from which row k of
% the output was obtained.
%
% EXAMPLE
%   intmat = [0,4,0,3,1,4; 0,4,1,3,0,4; 3,4,1,3,0,4];
%   [intmat,idcs,is_ord] = canonicalize_dir(intmat,false)
%   % returns intmat = [0,4,0,3,1,4; 0,2,0,3,0,4], idcs = [1;2] and
%   % is_ord = [1;0;0]: the first row is ordered, the second is replaced
%   % by [0,2,0,3,0,4], and the third has an empty region, and is removed.

q = size(intmat,1);
d = size(intmat,2)/2;
idcs = (1:q)';
is_ord = true(q,1);

% Only rows in which some variable has two distinct variable limits can be
% unordered
Llim = intmat(:,1:2:end);       Ulim = intmat(:,2:2:end);
cand = find(any(Llim>=1 & Llim<=d & Ulim>=1 & Ulim<=d & Llim~=Ulim,2))';
if isempty(cand)
    return
end
for r=cand
    is_ord(r) = ~any(unordered_vars(intmat(r,:)));
end
if check_only || all(is_ord)
    return
end

rows_new = num2cell(intmat,2);
for r=find(~is_ord)'
    rows_new{r} = split_row(intmat(r,:));
end
n_new = cellfun(@(a) size(a,1),rows_new);
idcs = repelem((1:q)',n_new);
intmat = vertcat(rows_new{:});
if isempty(intmat)
    intmat = zeros(0,2*d);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [bad_vars,Lr,Ur] = unordered_vars(row)
% [BAD_VARS,LR,UR] = UNORDERED_VARS(ROW) determines which variables in the
% row of an intmat have two variable limits for which the order is not
% implied by the limits of the other variables.
%
% INPUTS
% - row:        1 x 2*d array of integers specifying the limits of the
%               intervals of integration;
%
% OUTPUTS
% - bad_vars:   1 x d logical, with bad_vars(j) true if the order of the
%               limits of variable j is not implied;
% - Lr,Ur:      1 x d arrays of integers with the lower and upper limits of
%               the variables, in which each multiplier is replaced by the
%               variable (or boundary) at which it is evaluated;
%
% EXAMPLE
%   [bad_vars,Lr,Ur] = unordered_vars([0,4,1,3,0,4])
%   % returns bad_vars = [0,1,0], as x_2 is in [x_1,x_3] with independent
%   % x_1,x_3, and Lr = [0,1,0], Ur = [4,3,4]. For [3,4,0,4,2,4,0,5], x_1
%   % is in [x_3,x_4], but x_3<=x_4 follows from the upper limit of x_3,
%   % so bad_vars = [0,0,0,0]. For [0,5,4,3,0,5,1,1], x_4 is a multiplier,
%   % equal to x_1, so that Lr = [0,1,0,1], and bad_vars = [0,1,0,0].

d = numel(row)/2;
L = row(1:2:end);       U = row(2:2:end);
is_dlt = L==U;
% Eliminate multipliers: tgt(k+1) is the variable that variable k is equal
% to (following chains), tgt(1)=0 and tgt(d+2)=d+1 for the boundaries
tgt = 0:d+1;
tgt(find(is_dlt)+1) = L(is_dlt);
for k=1:d
    tgt = tgt(tgt+1);
end
Lr = tgt(L+1);          Ur = tgt(U+1);
is_varL = Lr>=1 & Lr<=d;
is_varU = Ur>=1 & Ur<=d;
% Variables of which the order of the limits must be implied
is_amb = ~is_dlt & is_varL & is_varU & Lr~=Ur;
% LE(x,y) true if x_x<=x_y is known. Start from the limits of the variables
% that need no proof, and add variables for which the order is implied
is_ok = ~is_amb & ~is_dlt;
bad_vars = is_amb;
for k=1:d
    LE = eye(d)>0;
    jL = find(is_ok & is_varL);     LE(sub2ind([d,d],Lr(jL),jL)) = true;
    jU = find(is_ok & is_varU);     LE(sub2ind([d,d],jU,Ur(jU))) = true;
    for l=1:ceil(log2(d))+1         % transitive closure
        LE = (double(LE)*double(LE))>0;
    end
    jA = find(bad_vars);
    is_new = LE(sub2ind([d,d],Lr(jA),Ur(jA)));
    if ~any(is_new)
        break
    end
    is_ok(jA(is_new)) = true;
    bad_vars(jA(is_new)) = false;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function rows = split_row(row)
% ROWS = SPLIT_ROW(ROW) splits a row of an intmat that is not ordered into
% ordered rows, see canonicalize_intmat. ROWS has no rows if the region of
% integration is empty.
%
% EXAMPLE
%   rows = split_row([3,4,3,4,0,5,0,5])
%   % returns the rows [0,2,0,4,0,1,0,5] and [0,4,0,1,0,2,0,5], for the
%   % orderings x_3<=x_1<=x_2<=x_4 and x_3<=x_2<=x_1<=x_4, whereas
%   % split_row([3,4,1,3,0,4]) returns zeros(0,6): the region is empty.

d = numel(row)/2;
[bad_vars,Lr,Ur] = unordered_vars(row);
% Variables to order: the variables that are not ordered, and the
% variables defining their limits, recursively
in_S = bad_vars;
for k=1:d
    lims = [Lr(in_S),Ur(in_S)];
    in_new = in_S;      in_new(lims(lims>=1 & lims<=d)) = true;
    if isequal(in_new,in_S)
        break
    end
    in_S = in_new;
end
var_nums = find(in_S);
n = numel(var_nums);
% Orderings of these variables consistent with their limits, ords(r,k) is
% the position of variable var_nums(k) in ordering r
ords = enumerate_valid_orders(Lr,Ur,in_S,var_nums,n);
nrows = size(ords,1);
if nrows==0
    rows = zeros(0,2*d);
    return
end
% Integrate each variable from a up to the next larger variable, the last
% over [a,b]
% var_idcs(r,:) = [k1,...,kn] implies x_k1 <= ... <= x_kn in ordering r
[~,var_idcs] = sort(ords,2);
% Ulim_tmp(r,k): upper limit of the variable at position k in ordering r
Ulim_tmp = [var_nums(var_idcs(:,2:end)),(d+1)*ones(nrows,1)];
Ulim_new = Ulim_tmp(sub2ind(size(ords),repmat((1:nrows)',1,n),ords));
rows = repmat(row,nrows,1);
rows(:,2*var_nums-1) = 0;
rows(:,2*var_nums) = Ulim_new;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %