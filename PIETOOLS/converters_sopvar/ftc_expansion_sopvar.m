function [Pop_f2x, Pop_b2x, bc_state_tab] = ftc_expansion_sopvar(PDE)
% Construct 'copvar' operators Pop_f2x and Pop_b2x such that the PDE state x
% can be expressed in terms of the associated fundamental state xf and a
% set of core boundary values xb as
%   x = Pop_f2x * xf + Pop_b2x * xb;
% for a PDE in any number of spatial variables.
%
% INPUTS:
% - PDE:    An initialized 'pde_struct' object defining a PDE, with state
%           components ordered by the spaces they exist in (see
%           "@pde_struct/reorder_comps").
%
% OUTPUTS:
% - Pop_f2x:    'copvar' describing the map from the fundamental state xf
%               to the PDE state x defined by PDE.x_tab. Both state spaces
%               are ordered by the binary index
%                   sum_j 2^(j-1)*[component varies in PDE.vars(j,1)],
%               and only spaces containing a state component are present.
% - Pop_b2x:    'copvar' describing the map from the core boundary state xb
%               to the PDE state x. The core boundary state is ordered
%               by the same binary index of the variables it still varies
%               in.
% - bc_state_tab:   A Nx(2+2*nvars) array defining the core boundary state
%               exactly as in "PIETOOLS_FTC_Expansion_2D": for each
%               boundary state, the component index, size, the variables
%               it varies in, and the order of the derivative wrt each
%               variable.
%
% NOTES:
% - Each variable s_j in which a state component x(S) is differentiable up
%   to order d_j is expanded as
%       x = sum_{k<d_j} (s_j-a_j)^k/k! d^k x(a_j) + int_{a_j}^{s_j} (s_j-t)^{d_j-1}/(d_j-1)! d^{d_j}x(t) dt.
%   Expanding in all variables, a core boundary state is associated to each
%   non-empty set A of differentiable variables, evaluated at the lower
%   boundary, and each combination of derivative orders k_i<d_i for i in A.
%   This state varies in the variables outside A, and is differentiated
%   d_j times wrt each of these. The state with A empty is the fundamental
%   state. The states are ordered as in "PIETOOLS_FTC_Expansion_2D".
% - The construction does not depend on the number of spatial variables;
%   all operators are expressed using the PI operators of the 'sopvar'
%   class, with arbitrary numbers of variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ftc_expansion_sopvar
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
% Initial coding: Claude 10/28/2026, adapting PIETOOLS_FTC_Expansion_2D to 'copvar' and an
% arbitrary number of variables

dom = PDE.dom;
vars = PDE.vars;
x_tab = PDE.x_tab;
nvars = size(vars,1);

dep_tab = logical(x_tab(:,3:2+nvars));
diff_tab = x_tab(:,3+nvars:2+2*nvars);
ncomps = size(x_tab,1);
np_arr = x_tab(:,2);

% Sort the variables by name, as required for the 'copvar' registry
names = pvar2varname(vars(:,1));    names = names(:)';
[~,order] = sort(names);            % registry position -> variable index
inv_order = zeros(1,nvars);         inv_order(order) = 1:nvars;

% The binary index of the space of each component
pow2 = 2.^(0:nvars-1)';
if nvars>0
    xmask = dep_tab*pow2;
else
    xmask = zeros(ncomps,1);
end

% Rows of the PDE state associated to each component, within its space
x_masks = unique(xmask)';
x_dim = zeros(size(x_masks));
x_rows = cell(ncomps,1);
for ii=1:numel(x_masks)
    for comp=find(xmask==x_masks(ii))'
        x_rows{comp} = x_dim(ii) + (1:np_arr(comp));
        x_dim(ii) = x_dim(ii) + np_arr(comp);
    end
end

% % % Enumerate the core boundary states, in the order of the 2D code
bc_list = struct('comp',{},'A',{},'kmat',{},'cmask',{});
for comp=1:ncomps
    isvar = dep_tab(comp,:);
    Dval = diff_tab(comp,:);
    isdiff = isvar & Dval>0;
    dvars = find(isdiff);
    for kk=1:2^numel(dvars)-1
        % Variables that are evaluated at the lower boundary: the binary
        % digits of kk, for the differentiable variables in increasing order
        A = false(1,nvars);
        A(dvars(logical(bitget(kk,1:numel(dvars))))) = true;
        % Variable that the boundary state still varies in
        is_intrr = isvar & ~A;
        % All combinations of derivative orders, first variable fastest
        Avars = find(A);
        dA = Dval(Avars);
        sub = cell(1,numel(Avars));
        [sub{:}] = ind2sub([dA,1],(1:prod(dA))');
        kmat = zeros(prod(dA),nvars);
        for jj=1:numel(Avars)
            kmat(:,Avars(jj)) = sub{jj}-1;
        end
        bc_list(end+1) = struct('comp',comp,'A',A,'kmat',kmat, ...
                    'cmask',is_intrr*pow2); %#ok<AGROW>
    end
end
if nvars==0 || isempty(bc_list)
    b_masks = zeros(1,0);
else
    b_masks = unique([bc_list.cmask]);
end
b_dim = zeros(size(b_masks));
for ii=1:numel(bc_list)
    bidx = find(b_masks==bc_list(ii).cmask);
    nst = size(bc_list(ii).kmat,1)*np_arr(bc_list(ii).comp);
    bc_list(ii).cols = b_dim(bidx) + (1:nst);
    b_dim(bidx) = b_dim(bidx) + nst;
end

% % % Assemble the blocks of the operators
nx = numel(x_masks);        nb = numel(b_masks);
% Parameters of the blocks, as arrays of polynomials of the right size
Cf = cell(nx,nx);
Cb = cell(nx,nb);
bc_xmask = xmask([bc_list.comp]);       bc_xmask = bc_xmask(:)';
bc_cmask = [bc_list.cmask];
for ii=1:nx
    Cf{ii,ii} = new_params(sum(dep_mask(x_masks(ii),nvars)),x_dim(ii),x_dim(ii));
    for jj=1:nb
        if any(bc_cmask==b_masks(jj) & bc_xmask==x_masks(ii))
            Cb{ii,jj} = new_params(sum(dep_mask(b_masks(jj),nvars)),x_dim(ii),b_dim(jj));
        end
    end
end

for comp=1:ncomps
    ii = find(x_masks==xmask(comp));
    isvar = dep_tab(comp,:);
    Dval = diff_tab(comp,:);
    isdiff = isvar & Dval>0;
    nxc = np_arr(comp);
    rows = x_rows{comp};

    % Factors for the interior (integral) variables
    % - fundamental state: integral over all variables that are
    %   differentiable, multiplier along the rest.
    [par,fctr] = interior_factor(isvar,isdiff,Dval,vars,[],names,order);
    Cf{ii,ii}{par} = assign(Cf{ii,ii}{par},rows,rows,kron(fctr,eye(nxc)));
end

for ib=1:numel(bc_list)
    comp = bc_list(ib).comp;
    A = bc_list(ib).A;
    isvar = dep_tab(comp,:);
    Dval = diff_tab(comp,:);
    isdiff = isvar & Dval>0;
    nxc = np_arr(comp);
    rows = x_rows{comp};
    ii = find(x_masks==xmask(comp));
    jj = find(b_masks==bc_list(ib).cmask);
    % Factors in the variables evaluated at the lower boundary
    kmat = bc_list(ib).kmat;
    fctr = polynomial(ones(1,size(kmat,1)));
    for v=find(A)
        fctr = fctr.*((1./factorial(kmat(:,v)')).*(vars(v,1)-dom(v,1)).^(kmat(:,v)'));
    end
    % Factors for the interior variables
    [par,ifctr] = interior_factor(isvar,isdiff,Dval,vars,A,names,order);
    fctr = ifctr*fctr;
    Cb{ii,jj}{par} = assign(Cb{ii,jj}{par},rows,bc_list(ib).cols,kron(fctr,eye(nxc)));
end

% % % Convert the polynomial parameters to 'sopvar' blocks
reg = names(order);                     % sorted registry
regdom = dom(order,:);
space_x = false(nx,nvars);
for ii=1:nx
    sp = dep_mask(x_masks(ii),nvars);
    space_x(ii,:) = sp(order);
end
space_b = false(nb,nvars);
for jj=1:nb
    sp = dep_mask(b_masks(jj),nvars);
    space_b(jj,:) = sp(order);
end

Pf = cell(nx,nx);
Pb = cell(nx,nb);
for ii=1:nx
    sp = dep_mask(x_masks(ii),nvars);
    Pf{ii,ii} = make_block(Cf{ii,ii},sp,sp,vars,dom,names,[x_dim(ii),x_dim(ii)]);
    for jj=1:nb
        if ~isempty(Cb{ii,jj})
            sb = dep_mask(b_masks(jj),nvars);
            Pb{ii,jj} = make_block(Cb{ii,jj},sp,sb,vars,dom,names,[x_dim(ii),b_dim(jj)]);
        end
    end
end

meta_f = struct('vars',{reg},'dom',regdom,'space_out',space_x,'space_in',space_x, ...
                'dim_out',x_dim(:),'dim_in',x_dim(:));
meta_b = struct('vars',{reg},'dom',regdom,'space_out',space_x,'space_in',space_b, ...
                'dim_out',x_dim(:),'dim_in',b_dim(:));
Pop_f2x = copvar(Pf,meta_f);
Pop_b2x = copvar(Pb,meta_b);

% % % The table describing the core boundary states
bc_state_tab = zeros(0,2+2*nvars);
sorted = stable_sort_idx([bc_list.cmask]);
for ib=sorted(:)'
    comp = bc_list(ib).comp;
    kmat = bc_list(ib).kmat;
    is_intrr = dep_tab(comp,:) & ~bc_list(ib).A;
    kdeg = kmat;
    kdeg(:,is_intrr) = repmat(diff_tab(comp,is_intrr),[size(kmat,1),1]);
    rows_tab = [repmat([x_tab(comp,1),np_arr(comp)],[size(kmat,1),1]), ...
                repmat(is_intrr,[size(kmat,1),1]), kdeg];
    bc_state_tab = [bc_state_tab; rows_tab]; %#ok<AGROW>
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function idx = stable_sort_idx(v)
% Indices sorting the vector v in ascending order, preserving the order of
% equal elements.
[~,idx] = sortrows([v(:),(1:numel(v))']);
end

function in = dep_mask(m,nvars)
% Logical indices of the variables included in the space with binary index m.
in = logical(bitget(m,1:nvars));
end

function C = new_params(n3,q,p)
% Cell array of zero polynomial matrices, parameters of a sopvar with n3
% shared variables.
if n3==0,       sz = [1,1];
elseif n3==1,   sz = [3,1];
else,           sz = 3*ones(1,n3);
end
C = cell(sz);
for k=1:numel(C)
    C{k} = polynomial(zeros(q,p));
end
end

function M = assign(M,rows,cols,val)
M(rows,cols) = val;
end

function [par,fctr] = interior_factor(isvar,isdiff,Dval,vars,A,names,order)
% Index of the parameter in a block with the variables outside A as shared
% variables, and the corresponding kernel factor.
% The kernel is the product over variables in the interior of
% (s-t)^(D-1)/(D-1)! for differentiable variables, and 1 otherwise.
nvars = size(vars,1);
if isempty(A),  A = false(1,nvars);     end
is_intrr = isvar & ~A;
fctr = polynomial(1);
for v=find(is_intrr & isdiff)
    fctr = fctr*((vars(v,1)-vars(v,2))^(Dval(v)-1)/factorial(Dval(v)-1));
end
% Shared variables sorted by name; index 2 if integrated, 1 otherwise
shared = order(ismember(order,find(is_intrr)));
n3 = numel(shared);
sub = ones(1,max(n3,1));
for k=1:n3
    sub(k) = 1 + double(isdiff(shared(k)));
end
if n3<=1
    par = sub(1);
else
    sub_c = num2cell(sub);
    par = sub2ind(3*ones(1,n3),sub_c{:});
end
end

function P = make_block(C,sp_out,sp_in,vars,dom,names,dims)
% Build a sopvar from the parameters C of a block between the spaces defined
% by the logical indices sp_out and sp_in.
nvars = numel(sp_out);
[~,order] = sort(names);
S3 = order(sp_out(order) & sp_in(order));
S2 = order(sp_out(order) & ~sp_in(order));
S1 = order(~sp_out(order) & sp_in(order));
n3 = numel(S3);
if n3==0,       sz = [1,1];
elseif n3==1,   sz = [3,1];
else,           sz = 3*ones(1,n3);
end
K = reshape(C,sz);
% The variables in S1 appear in the kernels as dummy variables, those in
% S2 as primary variables; for S3, kernels already use both.
for k=1:numel(K)
    for v=S1
        K{k} = subs(K{k},vars(v,1),vars(v,2));
    end
end
P = sopvar_from_kernels(K,S2,S3,S1,vars,dom,dims);
end
