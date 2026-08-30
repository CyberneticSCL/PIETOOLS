function [prog,Pop,Qcell,alpha_list] = possopvar(prog,dim,vars,dom,deg,options)
% [prog,Pop,Qcell,alpha_list] = POSSOPVAR(prog,dim,vars,dom,deg,options)
% declares a positive semidefinite, self-adjoint 'sdopvar' decision operator
%
%       Pop: L_2^m[S3] -> L_2^m[S3],       Pop = Pop* >= 0
%
% following Sec. 9 of the sopvar document. The operator is parameterized as
%
%       Pop = sum_{i,j} (Z_{alpha_i})* Q_{ij} (Z_{alpha_j}),     Q >= 0
%
% where each basis operator carries a single multi-index alpha and acts as
%
%       (Z_alpha x)(s) = int_theta I_alpha(s-theta) ...
%                           (I_m kron Z^alpha(theta,s)) x(theta) dtheta,
%
%       I_alpha(s) = prod_k I_{alpha_k}(s_k),
%
%       I_0(r) = delta(r),   I_1(r) = 1 for r>=0,   I_-1(r) = 1 for r<=0.
%
% Composing the adjoint pair and eliminating the integration variable theta
% by the semiseparable integral of Sec. 6.1 gives the sopvar kernel form
%
%       (Pop x)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds'
%
%       K_gamma(s,s') = (I_m kron ZL(s))' C_gamma(d) (I_m kron ZR(s'))
%
% with C_gamma(d) = unvec(A_gamma + B_gamma'*d) affine in the entries d of Q.
%
% For a single spatial variable this reproduces the L_2 -> L_2 block of
% 'poslpivar': alpha=1 is the multiplier block Z1(s), alpha=2 is the lower
% integral block Z2(eta,s), and alpha=3 is the upper integral block
% Z3(eta,s).
%
% INPUT
% - prog:   'struct' specifying an LPI/SOS program to modify;
% - dim:    scalar m, or 1x2 array [m,m], specifying the dimension of the
%           function space L_2^m. Since Pop is self-adjoint the input and
%           output dimensions must agree;
% - vars:   1 x n3 'cellstr' object specifying the names of the spatial
%           variables, or a scalar n3, in which case the variables are named
%           's1',...,'sn3'. These are the variables common to the input and
%           output space, i.e. the set S3 of the sopvar document;
% - dom:    n3 x 2 array of type 'double', with dom(k,:) = [ak,bk] the
%           domain of variable k. May be given as a single 1x2 row, in which
%           case the same domain is used for every variable;
% - deg:    degrees of the monomial bases Z^alpha(theta,s). May be given as
%             * a scalar d, in which case every basis has degree at most d
%               in each theta_k and each s_k;
%             * a 'struct' with any of the fields
%                   deg.int:    scalar or 1 x n3, maximal degree in each
%                               integration variable theta_k. Defaults to 1;
%                   deg.mult:   scalar or 1 x n3, maximal degree in each
%                               output variable s_k. Defaults to deg.int;
%                   deg.joint:  scalar, maximal total degree of theta and s
%                               combined. Defaults to no further restriction;
%               For one spatial variable, (deg.int,deg.mult,deg.joint)
%               correspond to the elements of 'poslpivar's d{2};
%             * a cell array with one such scalar or struct per included
%               basis operator, in the order of the rows of alpha_list;
%           Note that Z^alpha is taken to have degree 0 in s_k whenever
%           alpha_k = 1, since the factor delta(s_k-theta_k) identifies the
%           two variables and any s_k dependence would be redundant;
% - options: (optional) 'struct' specifying other options, with fields
%   options.psatz     binary value, set to 1 to enforce positivity of the
%                     operator only on the domain, by including the factor
%                     g(theta) = prod_k (theta_k-ak)*(bk-theta_k) in the
%                     construction. Defaults to 0;
%   options.include   specification of which basis operators to include.
%                     Either an N x n3 array whose rows are the desired
%                     multi-indices alpha (entries in {1,2,3}, with 1 the
%                     multiplier, 2 the lower integral and 3 the upper
%                     integral), or a logical array of length 3^n3 masking
%                     the multi-indices in their standard linear order, or a
%                     numeric vector of linear indices into that order.
%                     Defaults to all 3^n3 multi-indices;
%
% OUTPUT
% - prog:       'struct' specifying the same program as the input, but now
%               including the decision variables defining Pop and a
%               constraint enforcing Q>=0;
% - Pop:        m x m 'sdopvar' object representing a positive semidefinite,
%               self-adjoint PI operator decision variable;
% - Qcell:      N x N cell array of 'cellstr' objects naming the decision
%               variables in each block Q_{ij} of the matrix Q>=0, as
%               returned by 'sosquadvar';
% - alpha_list: N x n3 array whose rows are the multi-indices alpha of the
%               included basis operators, in the order used for Qcell;
%
% NOTES
% This implements only the L_2 -> L_2 part of the construction, i.e. the
% reduction of Sec. 9.1 in which the created and lost variables S1 and S2
% are ignored. A positive operator on R^n x L_2^m[S3] additionally requires
% the R^n -> R^n, L_2 -> R^n and R^n -> L_2 blocks, which are formed from
% the same Q>=0 but need a container object to hold them.
%
% Following Sec. 9.3, the products (Z_{alpha_i})* Q_{ij} (Z_{alpha_j}) are
% formed as 'dpvar' objects via 'sosquadvar' and then converted to the
% kernel format with 'dpvar2sdvar'. This reuses existing code at the cost of
% carrying three groups of variables through the 'dpvar' arithmetic, and is
% the temporary measure the document describes; the coefficient matrices
% B_gamma could instead be built directly from vec(Q) by index arithmetic.
%
% See also POSLPIVAR, SOSQUADVAR, INT_SEMISEP, DPVAR2SDVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - possopvar
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
% MP, 08/22/2026: Initial coding


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Matrix dimension
if isscalar(dim)
    m = dim;
elseif numel(dim)==2
    if dim(1)~=dim(2)
        error("A self-adjoint operator must have equal input and output dimensions.")
    end
    m = dim(1);
else
    error("Matrix dimension of the operator should be specified as a scalar or 1x2 array.")
end
if m<1 || m~=round(m)
    error("Matrix dimension of the operator should be a positive integer.")
end

% % % Spatial variables
if isnumeric(vars) && isscalar(vars)
    n3 = vars;
    vars = cell(1,n3);
    for k=1:n3
        vars{k} = ['s',num2str(k)];
    end
elseif isa(vars,'polynomial')
    vars = vars.varname(:)';
elseif ischar(vars)
    vars = {vars};
end
if ~iscellstr(vars)
    error("Spatial variables should be specified as a 1 x n3 'cellstr' object.")
end
vars = vars(:)';        n3 = numel(vars);
if numel(unique(vars))~=n3
    error("Spatial variable names should be unique.")
end

% Declare names for the integration variable theta and the input dummy
% variable s'. Neither appears in the returned object: the dummy variables
% of an sdopvar are implicit in vars.in.
% These suffixes are reserved outright, rather than only where they would
% actually collide, so that the set of admissible names does not depend on
% which other variables happen to be present.
vars_int = strcat(vars,'_int');
vars_dum = strcat(vars,'_dum');
is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end

% % % Domain
if isempty(dom) && n3==0
    dom = zeros(0,2);
end
if size(dom,2)~=2
    error("Domains should be specified as an n3 x 2 array.")
end
if size(dom,1)==1 && n3~=1
    dom = repmat(dom,n3,1);
elseif size(dom,1)~=n3
    error("Domains should be specified as an n3 x 2 array for n3 spatial variables.")
end
if any(dom(:,2)<=dom(:,1))
    error("Each domain should satisfy dom(k,1) < dom(k,2).")
end

% % % Options
if nargin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if ~isfield(options,'psatz') || isempty(options.psatz)
    options.psatz = 0;
end
if isfield(options,'sep') && ~isempty(options.sep) && options.sep
    error("The 'sep' option is not supported by 'possopvar'.")
end

% % % Basis operators to include
%   1 <-> multiplier (delta),  2 <-> lower integral,  3 <-> upper integral
if n3==0
    alpha_all = zeros(1,0);
else
    alpha_all = fliplr(dec2base(0:3^n3-1,3,n3)-'0')+1;
end
if ~isfield(options,'include') || isempty(options.include)
    alpha_list = alpha_all;
else
    incl = options.include;
    if islogical(incl)
        if numel(incl)~=size(alpha_all,1)
            error("A logical 'include' should have 3^n3 elements.")
        end
        alpha_list = alpha_all(incl(:),:);
    elseif n3>0 && size(incl,2)==n3
        if any(incl(:)<1) || any(incl(:)>3) || any(incl(:)~=round(incl(:)))
            error("Multi-indices in 'include' should have entries 1, 2 or 3.")
        end
        alpha_list = incl;
    else
        if any(incl(:)<1) || any(incl(:)>size(alpha_all,1)) ...
                || any(incl(:)~=round(incl(:)))
            error("Linear indices in 'include' should lie between 1 and 3^n3.")
        end
        alpha_list = alpha_all(incl(:),:);
    end
end
% Note that for n3=0 there is exactly one basis operator, indexed by the
% empty multi-index, so the row count rather than 'isempty' must be tested.
if size(alpha_list,1)==0
    error("At least one basis operator must be included.")
end
if size(unique(alpha_list,'rows'),1)~=size(alpha_list,1)
    error("Multi-indices in 'include' should be distinct.")
end
nblk = size(alpha_list,1);

% % % Degrees
deg_list = process_degrees(deg,nblk,n3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build the monomial bases Z^alpha(theta,s) and Z^alpha(theta,s')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The exponents of Z^alpha are stored over the variable list
% [theta_1,...,theta_n3, s_1,...,s_n3], so that Z1c{i} and Z2c{i} differ
% only in the name of the second group of variables.
Z1c = cell(1,nblk);
Z2c = cell(1,nblk);
mdim = m*ones(nblk,1);
ndim = mdim;
for i=1:nblk
    caps_int = deg_list{i}.int;
    caps_mult = deg_list{i}.mult;
    % A multiplier in variable k identifies s_k with theta_k, so any s_k
    % dependence of the basis would be redundant.
    caps_mult(alpha_list(i,:)==1) = 0;

    Ei = build_exponent_grid([caps_int,caps_mult],deg_list{i}.joint);
    Ti = size(Ei,1);

    Z1c{i} = polynomial(speye(Ti),Ei,[vars_int(:);vars(:)],[Ti,1]);
    Z2c{i} = polynomial(speye(Ti),Ei,[vars_int(:);vars_dum(:)],[Ti,1]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declare Q>=0 and form the products N{i,j} = Z^alpha_i' Q_ij Z^alpha_j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[prog,N,Qcell] = sosquadvar(prog,Z1c,Z2c,mdim,ndim,'pos');

% Collect the full list of decision variables, so that every block is
% expressed over a common decision variable basis.
Zd = {};
for k=1:numel(Qcell)
    Zd = [Zd; reshape(cellstr(string(Qcell{k})),[],1)];                     %#ok<AGROW>
end
Zd = unique(Zd);
ndec = numel(Zd);
% Index lookup for the decision variables, built once. Each block below
% needs the positions of its own names in Zd; doing that with ismember
% per block is O(ndec) per block, which is the dominant cost once there
% are many decision variables.
dmap = containers.Map(Zd,num2cell(1:ndec));                                % MMP, 08/29/2026

% Multiplier used to restrict positivity to the domain. Note that it is
% evaluated at the integration variable, as in 'poslpivar'.
if options.psatz
    gfun = polynomial(1);
    for k=1:n3
        thk = polynomial(vars_int(k));
        gfun = gfun*(thk-dom(k,1))*(dom(k,2)-thk);
    end
else
    gfun = polynomial(1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eliminate the integration variable block by block
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variable split used to read the products as
%   N{i,j} = (I_m kron Zs(s) kron Zth(theta))' C (I_m kron Zsp(s'))
% Placing theta last in vars.out makes the theta monomials the fastest row
% index, which is the layout 'int_semisep' expects of its input.
vars_conv = struct();
vars_conv.out = [vars,vars_int];
vars_conv.in = vars_dum;

vars_io = struct();
vars_io.in = vars;
vars_io.out = vars;
dom_io = struct();
dom_io.in = dom;
dom_io.out = dom;

% Adjoining a basis operator flips its lower and upper integrals.
negmap = [1,3,2];

Pop = [];
for i=1:nblk
    for j=1:nblk
        % % % Coefficients of the product, over the split (s,theta | s')
        Rij = gfun*N{i,j};
        [Pblk,ZLb,ZRb] = dpvar2sdvar(Rij,vars_conv);

        Zs  = ZLb(1:n3);
        Zth = ZLb(n3+1:2*n3);
        Zsp = ZRb;
        Ns  = prod(cellfun(@numel,Zs));
        NG  = prod(cellfun(@numel,Zth));
        Nsp = prod(cellfun(@numel,Zsp));
        g1 = m*Ns;      g2 = m*Nsp;
        if Pblk.m~=g1*NG || Pblk.n~=g2
            error("Internal error: unexpected coefficient dimensions.")
        end

        % Express the affine coefficient A + B'*d over the common decision
        % variable basis Zd.
        Bblk = remap_dvars(Pblk.B,cellstr(string(Pblk.dvarname)),dmap,ndec);   % MMP, 08/29/2026

        % % % Eliminate theta
        % Both A and each row of B are subjected to the same linear map, so
        % they are stacked as additional column blocks of the input and
        % sliced apart afterwards.
        G = struct();
        G.C = pack_sheets([Pblk.A.';Bblk],Pblk.m,Pblk.n);
        G.Z = Zth;
        [Cgam,ZLnew,ZRnew] = int_semisep(G,negmap(alpha_list(i,:)),alpha_list(j,:),dom);
        NL = prod(cellfun(@numel,ZLnew));
        NR = prod(cellfun(@numel,ZRnew));

        % % % Merge the pre-existing and newly generated monomials
        [ZLf,ML] = merge_monomial_product(Zs,ZLnew);
        [ZRf,MR] = merge_monomial_product(Zsp,ZRnew);
        Lmat = kron(speye(m),ML).';
        Rmat = kron(speye(m),MR);

        params = struct();
        params.A = cell(numel(Cgam),1);
        params.B = cell(numel(Cgam),1);
        for k=1:numel(Cgam)
            [Ak,Bk] = unpack_sheets(Cgam{k},g1*NL,g2*NR,ndec);
            [params.A{k},params.B{k}] = lr_multiply(Lmat,Ak,Bk,Rmat);
        end
        params.A = reshape(params.A,[3*ones(1,n3),1,1]);
        params.B = reshape(params.B,[3*ones(1,n3),1,1]);

        Pij = sdopvar(params,vars_io,Zd,ZLf,ZRf,dom_io,[m,m]);
        if isempty(Pop)
            Pop = Pij;
        else
            Pop = Pop+Pij;
        end
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg_list = process_degrees(deg,nblk,n3)
% Expand the degree specification into one struct per basis operator, each
% with fields 'int', 'mult' (1 x n3 arrays) and 'joint' (scalar).

if iscell(deg)
    if numel(deg)~=nblk
        error("A cell 'deg' should have one entry per included basis operator.")
    end
    deg_list = cell(1,nblk);
    for i=1:nblk
        deg_list{i} = process_degrees_one(deg{i},n3);
    end
else
    deg_list = repmat({process_degrees_one(deg,n3)},1,nblk);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = process_degrees_one(deg,n3)
% Expand a single degree specification.

if isnumeric(deg) && isscalar(deg)
    deg = struct('int',deg);
elseif ~isa(deg,'struct')
    error("Monomial degrees should be specified as a scalar or 'struct' object.")
end

if ~isfield(deg,'int') || isempty(deg.int)
    deg.int = 1;
end
if ~isfield(deg,'mult') || isempty(deg.mult)
    deg.mult = deg.int;
end

d = struct();
d.int = expand_caps(deg.int,n3,'int');
d.mult = expand_caps(deg.mult,n3,'mult');

if ~isfield(deg,'joint') || isempty(deg.joint)
    d.joint = sum(d.int)+sum(d.mult);
elseif ~isscalar(deg.joint)
    error("The joint degree should be specified as a scalar.")
else
    d.joint = deg.joint;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function caps = expand_caps(caps,n3,name)
% Expand a scalar or 1 x n3 degree bound into a 1 x n3 array.

caps = reshape(caps,1,[]);
if isscalar(caps)
    caps = repmat(caps,1,n3);
elseif numel(caps)~=n3
    error("The '"+string(name)+"' degree should be a scalar or have one "...
          +"element per spatial variable.")
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Bnew = remap_dvars(B,dvars_old,dmap,ndec)
% Move the rows of B from the ordering dvars_old to the global ordering that
% dmap indexes, inserting zero rows for decision variables absent from B.
% The lookup is a hash per name rather than a search over the whole global
% list, so the cost is in the number of names B actually carries.

loc = zeros(numel(dvars_old),1);
for i = 1:numel(dvars_old)
    if ~isKey(dmap,dvars_old{i})
        error("Internal error: unrecognized decision variable.")
    end
    loc(i) = dmap(dvars_old{i});
end

Bnew = sparse(ndec,size(B,2));
if ~isempty(loc)
    Bnew(loc,:) = B;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CG = pack_sheets(V,mrow,ncol)
% Reshape each row of V into an mrow x ncol matrix and lay the results out
% side by side, so that sheet sg occupies columns (sg-1)*ncol+1 : sg*ncol.

nsheet = size(V,1);
[sg,lin,val] = find(V);
irow = mod(lin-1,mrow)+1;
icol = floor((lin-1)/mrow)+1;
CG = sparse(irow,(sg-1)*ncol+icol,val,mrow,nsheet*ncol);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B] = unpack_sheets(C,mrow,ncol,nsheet)
% Inverse of 'pack_sheets' applied to the output of 'int_semisep': the
% column index of the result is (sheet,ncol,NR) with the semiseparable index
% fastest, so each sheet still occupies a contiguous block of ncol columns.
% Sheet 1 is returned as the vectorized A, the remaining nsheet sheets as
% the rows of B.
%
% The whole array is scanned once and then partitioned by column, rather
% than sliced first. Slicing a sparse matrix whose column count carries the
% decision variable dimension, as C(:,ncol+1:end) does, copies that block
% before anything useful happens, and the copy dominates once there are many
% decision variables.

if size(C,1)~=mrow || size(C,2)~=(1+nsheet)*ncol
    error("Internal error: unexpected coefficient dimensions.")
end

[irow,icol,val] = find(C);
is_A = icol<=ncol;

A = sparse((icol(is_A)-1)*mrow+irow(is_A),1,val(is_A),mrow*ncol,1);

jcol = icol(~is_A)-ncol;
sg = floor((jcol-1)/ncol)+1;
cc = jcol - (sg-1)*ncol;
B = sparse(sg,(cc-1)*mrow+irow(~is_A),val(~is_A),nsheet,mrow*ncol);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aout,Bout] = lr_multiply(L,A,B,R)
% Given vec(C(d)) = A + B'*d, return the coefficients of vec(L*C(d)*R), using
% vec(L*C*R) = (R' kron L)*vec(C).
%
% The constant term is handled by reshaping, which avoids the Kronecker
% product entirely. For B the product is formed explicitly, as in
% @sdopvar/plus: here L and R are monomial selection matrices with a single
% nonzero per row, so the Kronecker product has only nnz(L)*nnz(R) nonzeros
% and applying it to all decision variables at once is far cheaper than
% looping over the rows of B, of which there may be many thousands.

X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
Aout = Y(:);

Bout = B*kron(R.',L).';

end
