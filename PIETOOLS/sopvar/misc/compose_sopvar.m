function Yop = compose_sopvar(Xops,vars_in,vars_tot,dom,dim_in)
% YOP = COMPOSE_SOPVAR(XOPS,VARS_IN,VARS_TOT,DOM,DIM_IN) composes operators
% X_1,...,X_nX that each act on a single variable, with the identity along
% all other variables, into one 'sopvar' operator
%       Yop = X_nX*...*X_1,
% on the function space of the variables vars_in. For example, for X_1 an
% integral operator in s1 and X_2 a multiplier in s2,
%   (Yop*u)(s1,s2) = X_2(s2)*int_{a}^{s1} K(s1,theta)*u(theta,s2) dtheta.
%
% INPUT
% - Xops:   1 x nX struct array with fields
%       op:     'sopvar' object acting in a single variable, either as
%               a 3-PI operator (the variable is in its input and output
%               space), integrating the variable out (only in its input
%               space), or introducing it (only in its output space);
%       var:    char, the name of that variable;
%           The stages are applied in the order t=1,...,nX, and each
%           variable can occur in one stage at most. The input dimension of
%           stage t must match the output dimension of stage t-1;
% - vars_in:    cellstr, the variables of the input space of Yop: those
%           that a stage acts on as a 3-PI operator or integrates out, and
%           those that no stage acts on, along which Yop is the identity;
% - vars_tot:   cellstr, the names of all variables involved;
% - dom:    N x 2 array, the intervals on which the variables vars_tot are
%           defined;
% - dim_in: integer, the dimension of the input of Yop, only used if
%           Xops is empty, in which case Yop is the identity;
%
% OUTPUT
% - Yop:    'sopvar' object, from the input dimension of X_1 to the output
%           dimension of X_nX, with input variables vars_in, and output
%           variables those of vars_in that no stage integrates out, and
%           those that a stage introduces;
%
% NOTES
% - As each variable is acted upon by one stage only, the composed kernel
%   is the matrix product of the kernels of the stages,
%   K_nX(s_nX,theta_nX)*...*K_1(s_1,theta_1), so that the parameters of
%   Yop, for each choice of the partial integral of each stage, follow from
%   the product of those of the stages, contracted over the components,
%   with the monomials tensored. This avoids lifting each stage to all
%   variables and composing them with the 'mtimes' method of 'sopvar'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - compose_sopvar
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
% DJ, 11/01/2026: Initial coding;

% No stages: the identity on the input space
if isempty(Xops)
    [~,loc] = ismember(vars_in,vars_tot);
    Yop = mat2sopvar(eye(dim_in),struct('in',{vars_in},'out',{vars_in}),struct('in',dom(loc,:),'out',dom(loc,:)));
    return
end
% Classify the variable of each stage: it is in both the input and output
% space of the stage (3, e.g. integral operator), only in the output (2,
% the stage introduces the variable, e.g. evaluation of an auxiliary value
% as a function of s), or only in the input (1, the stage integrates it
% out)
n_ops = numel(Xops);
vars_stage = {Xops.var};
type_arr = zeros(1,n_ops);                      % 3: in and out, 2: out only, 1: in only
for j=1:n_ops
    hasin = ismember(vars_stage{j},Xops(j).op.vars.in);
    hasout = ismember(vars_stage{j},Xops(j).op.vars.out);
    type_arr(j) = 2*hasout + hasin;
    % The variable is in the input space, unless the stage introduces it
    if (hasin && ~ismember(vars_stage{j},vars_in)) || (~hasin && ismember(vars_stage{j},vars_in))
        error('The variable of a stage is not in the space on which it acts.')
    end
end
% Variables of the composed operator: those not acted on by any stage pass
% through unchanged (multiplier by identity, so also in S3)
vars_pass = vars_in(~ismember(vars_in,vars_stage));     % acted on by no stage
S3 = sort([vars_stage(type_arr==3),vars_pass]);
S2 = sort(vars_stage(type_arr==2));
S1 = sort(vars_stage(type_arr==1));
vars_out_new = [S2,S3];     vars_in_new = [S3,S1];
% Bases and domains of the output (side 1) and input (side 2) variables:
% the monomial basis in a variable is that of the stage acting on it, or
% {1} if there is none. nZ_cell: number of monomials in each variable
vars_new_cell = {vars_out_new,vars_in_new};
side_cell = {'out','in'};       Zfield_cell = {'ZL','ZR'};
Z_cell = cell(1,2);             dom_cell = cell(1,2);       nZ_cell = cell(1,2);
for side_idx=1:2
    vars_side = vars_new_cell{side_idx};
    Z_cell{side_idx} = repmat({0},1,numel(vars_side));
    dom_cell{side_idx} = zeros(numel(vars_side),2);
    for k=1:numel(vars_side)
        j = find(strcmp(vars_stage,vars_side{k}),1);
        dom_cell{side_idx}(k,:) = dom(strcmp(vars_tot,vars_side{k}),:);
        if ~isempty(j)
            Z_cell{side_idx}{k} = Xops(j).op.(Zfield_cell{side_idx}){strcmp(Xops(j).op.vars.(side_cell{side_idx}),vars_side{k})};
        end
    end
    nZ_cell{side_idx} = cellfun(@numel,Z_cell{side_idx});
end
ZL = Z_cell{1};         ZR = Z_cell{2};
dom_out = dom_cell{1};  dom_in = dom_cell{2};
nZL_arr = nZ_cell{1};   nZR_arr = nZ_cell{2};
% Dimensions and monomial counts of the stages (stage t maps
% q_arr(t)-valued functions to p_arr(t)-valued ones, so q_arr(t)=p_arr(t-1)
% for a valid chain), and the empty parameters of the composition, indexed
% by the slots of S3
p_arr = zeros(1,n_ops);         q_arr = zeros(1,n_ops);
nZL_stage = ones(1,n_ops);      nZR_stage = ones(1,n_ops);
for j=1:n_ops
    p_arr(j) = Xops(j).op.dims(1);
    q_arr(j) = Xops(j).op.dims(2);
    if type_arr(j)>=2
        nZL_stage(j) = nZL_arr(strcmp(vars_out_new,vars_stage{j}));
    end
    if mod(type_arr(j),2)==1
        nZR_stage(j) = nZR_arr(strcmp(vars_in_new,vars_stage{j}));
    end
    if j>1 && q_arr(j)~=p_arr(j-1)
        error('Dimensions of the operators in the chain do not match.')
    end
end
n3 = numel(S3);
nnZL_arr = prod(nZL_arr);       nnZR_arr = prod(nZR_arr);
params_new = cell([repmat(3,1,n3),1]);
for lidx=1:numel(params_new)
    params_new{lidx} = sparse(p_arr(n_ops)*nnZL_arr,q_arr(1)*nnZR_arr);
end
% Products of the parameters of the stages, for all partial integrals.
% parts(pidx) holds a partial product A_j*...*A_1 for one choice of slots
% (field slot), as a sparse matrix with rows (output component) and columns
% (monomial in the output variables, input component, monomial in the input
% variables): each stage multiplies on the left, and extends the monomial
% indices by the factor of its own variable (Kronecker structure). Stages
% whose parameter has no entries in a slot drop the corresponding term.
parts = struct('A',{[]},'slot',{ones(1,n3)},'nZ',1,'nZp',1);
for j=1:n_ops
    X = Xops(j).op;
    parts_new = parts([]);
    for pidx=1:numel(parts)
        for lidx=1:numel(X.params)
            [rridcs,ccidcs,vals] = find(X.params{lidx});
            rridcs = rridcs(:);     ccidcs = ccidcs(:);     vals = vals(:);
            if isempty(vals)
                continue
            end
            % Split the row and column indices of the stage into component
            % and monomial index (monomial fastest)
            ridcs = floor((rridcs-1)/nZL_stage(j))+1;
            ZL_idcs = mod(rridcs-1,nZL_stage(j))+1;
            cidcs = floor((ccidcs-1)/nZR_stage(j))+1;
            ZR_idcs = mod(ccidcs-1,nZR_stage(j))+1;
            part = parts(pidx);
            if j==1
                % First stage: initialize the product with its parameter
                ccidcs_new = ((ZL_idcs-1)*q_arr(1)+(cidcs-1))*nZR_stage(j)+ZR_idcs;
                A = sparse(ridcs,ccidcs_new,vals,p_arr(1),nZL_stage(1)*q_arr(1)*nZR_stage(1));
                nZ = nZL_stage(1);      nZp = nZR_stage(1);
            else
                % Later stages: multiply the previous product from the left
                % by the stage (components Mr), then re-index the new
                % monomials (ZL_idcs_new: output, ZR_idcs_new: input) with
                % the new variable as the fastest index
                Mr = sparse(((ridcs-1)*nZL_stage(j)+(ZL_idcs-1))*nZR_stage(j)+ZR_idcs,cidcs,vals,p_arr(j)*nZL_stage(j)*nZR_stage(j),p_arr(j-1));
                [rridcs_A,ccidcs_A,vals_A] = find(Mr*part.A);
                rridcs_A = rridcs_A(:);     ccidcs_A = ccidcs_A(:);     vals_A = vals_A(:);
                ZR_idcs_t = mod(rridcs_A-1,nZR_stage(j))+1;
                tmp_idcs = floor((rridcs_A-1)/nZR_stage(j));
                ZL_idcs_t = mod(tmp_idcs,nZL_stage(j))+1;
                ridcs_new = floor(tmp_idcs/nZL_stage(j))+1;
                ZR_idcs_prev = mod(ccidcs_A-1,part.nZp)+1;
                rem_idcs = floor((ccidcs_A-1)/part.nZp);
                cidcs_in = mod(rem_idcs,q_arr(1))+1;
                ZL_idcs_prev = floor(rem_idcs/q_arr(1))+1;
                ZL_idcs_new = (ZL_idcs_prev-1)*nZL_stage(j)+ZL_idcs_t;
                ZR_idcs_new = (ZR_idcs_prev-1)*nZR_stage(j)+ZR_idcs_t;
                nZ = part.nZ*nZL_stage(j);      nZp = part.nZp*nZR_stage(j);
                ccidcs_new = ((ZL_idcs_new-1)*q_arr(1)+(cidcs_in-1))*nZp+ZR_idcs_new;
                A = sparse(ridcs_new,ccidcs_new,vals_A,p_arr(j),nZ*q_arr(1)*nZp);
            end
            if nnz(A)==0
                continue
            end
            part.A = A;     part.nZ = nZ;       part.nZp = nZp;
            % Record the slot of this stage for its variable (S3 only)
            if type_arr(j)==3
                part.slot(strcmp(S3,vars_stage{j})) = lidx;
            end
            parts_new(end+1) = part;                                        %#ok<AGROW>
        end
    end
    parts = parts_new;
end
% Convert to the canonical monomial ordering of the output and input: the
% monomial indices above are ordered by stage, whereas the sopvar orders
% the variables as vars_out_new and vars_in_new (the first variable
% slowest). Then add the product to the parameter of its slot.
for pidx=1:numel(parts)
    [ridcs_new,ccidcs_A,vals_A] = find(parts(pidx).A);
    ridcs_new = ridcs_new(:);   ccidcs_A = ccidcs_A(:);     vals_A = vals_A(:);
    nZp = parts(pidx).nZp;
    ZR_idcs_rem = mod(ccidcs_A-1,nZp);
    rem_idcs = floor((ccidcs_A-1)/nZp);
    cidcs_in = mod(rem_idcs,q_arr(1))+1;
    ZL_idcs_rem = floor(rem_idcs/q_arr(1));
    ZL_idcs_arr = zeros(numel(vals_A),n_ops);       ZR_idcs_arr = zeros(numel(vals_A),n_ops);
    for j=n_ops:-1:1
        ZL_idcs_arr(:,j) = mod(ZL_idcs_rem,nZL_stage(j));
        ZL_idcs_rem = floor(ZL_idcs_rem/nZL_stage(j));
        ZR_idcs_arr(:,j) = mod(ZR_idcs_rem,nZR_stage(j));
        ZR_idcs_rem = floor(ZR_idcs_rem/nZR_stage(j));
    end
    mono_cell = {zeros(numel(vals_A),1),zeros(numel(vals_A),1)};
    digits_cell = {ZL_idcs_arr,ZR_idcs_arr};
    for side_idx=1:2
        for k=1:numel(vars_new_cell{side_idx})
            j = find(strcmp(vars_stage,vars_new_cell{side_idx}{k}),1);
            digit_idcs = 0;
            if ~isempty(j)
                digit_idcs = digits_cell{side_idx}(:,j);
            end
            mono_cell{side_idx} = mono_cell{side_idx}*nZ_cell{side_idx}(k)+digit_idcs;
        end
    end
    ZL_idcs_can = mono_cell{1};     ZR_idcs_can = mono_cell{2};
    lidx = 1;
    if n3>0
        slot_cell = num2cell(parts(pidx).slot);
        lidx = sub2ind([3*ones(1,n3),1],slot_cell{:});
    end
    rridcs_can = (ridcs_new-1)*nnZL_arr+ZL_idcs_can+1;
    ccidcs_can = (cidcs_in-1)*nnZR_arr+ZR_idcs_can+1;
    params_new{lidx} = params_new{lidx} + sparse(rridcs_can,ccidcs_can,vals_A,p_arr(n_ops)*nnZL_arr,q_arr(1)*nnZR_arr);
end
% Declare the composed operator
Yop = sopvar(params_new,struct('in',{vars_in_new},'out',{vars_out_new}),ZL,ZR,struct('in',dom_in,'out',dom_out),[p_arr(n_ops),q_arr(1)]);

end