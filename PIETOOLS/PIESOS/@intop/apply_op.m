function c = apply_op(Kop,v)
% C = APPLY_OP(KOP,V) applies the functional KOP to the function V,
% returning the result C
%
% INPUTS
% - Kop:    m x n 'intop' object acting on
%           L2^{n}[[a1,b1]^d1 x ... x [aN,bN]^dN];
% - v:      n x p 'polynomial' (or 'double') array, representing a
%           polynomial function in L2^{n}[[a1,b1]^d1 x ... x [aN,bN]^dN];
%
% OUTPUTS
% - c:      m x p 'polynomial' array (m x p 'dpvar' array if Kop depends on
%           decision variables), representing the output of Kop upon acting
%           on v;
%
% NOTES
% - The function v must be defined in terms of the spatial variables of the
%   functional, Kop.vars.in. If v also depends on other variables, so does
%   c.
% - v may also be given as a 1x1 cell {v}, the form in which a single
%   factor of a degree-one monomial is passed on by 'apply_op' of
%   'tensopvar' and 'tensopmat'.
% - With the kernels expressed in terms of the monomials z(s) defined by
%   Kop.Z.in, the integral over the region r of the kernel times v is
%       P_r * W_r,     W_r = int_r kron(I_n,z(s))*v(s) ds,
%   where P_r = Kop.params{r}. The integrals W_r do not depend on the
%   decision variables: only the final product of the sparse P_r with the
%   dense W_r involves them. As the limits of integration along a
%   variable only involve that variable, the integral of a monomial
%   z(s)*w(s) over the region r is the product of the integrals along the
%   different variables. Hence, W_r follows from the integrals along each
%   variable of the products of the monomials of z with those of v, which
%   are contracted with the coefficients of v one variable at a time.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - apply_op
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
% DJ, 10/15/2026: Initial coding;

% % Process the inputs
% Check that the functional is properly specified
if ~isa(Kop,'intop')
    error("First argument must be of type 'intop'.")
elseif ~Kop.chkval && ~isvalid(Kop)
    error("The 'intop' object is not properly declared.")
end
% Check that the function is properly specified
if iscell(v) && isscalar(v)
    v = v{1};
end
if isa(v,'double')
    v = polynomial(v);
elseif ~isa(v,'polynomial')
    error("Second argument must be of type 'polynomial'.")
end
% Check that the dimensions match
[m,n] = size(Kop);
if n~=size(v,1)
    error("Row dimension of function must match column dimension of functional.")
end
p = size(v,2);
nd = numel(Kop.Zd);

% % Define Kop and v in terms of the same variables
% Extract the variables of integration
varnames = Kop.vars.in;
N = numel(varnames);
varnames_in = Kop.varnames.in_dum;
varnames_out = Kop.varnames.in;
pvars = cellfun(@(a) polynomial(a(:)),varnames_in,'UniformOutput',false);

% Replace variables in v by dummy variables for integration
v_varname = v.varname;
for i=1:N
    vnames_old = varnames_out{i};   % spatial variables
    vnames_new = varnames_in{i};    % dummy variables

    [~,ia,ib] = intersect(v_varname,vnames_old,'stable');
    v_varname(ia) = vnames_new(ib);
end
% Also account for degree 1 case, e.g. s1 instead of s1_1
varnames_in1 = cellfun(@(a) a{1},varnames_in,'UniformOutput',false);
[~,ia,ib] = intersect(v_varname,Kop.vars.in,'stable');
v_varname(ia) = varnames_in1(ib);
v.varname = v_varname;

% % Apply the operator
if N==0
    % No spatial variables: Kop is a plain constant multiplier, which multiplies v
    % (v can still depend on variables other than the spatial variables)
    Kc = kernels(Kop);
    c = Kc{1}*v;
    return
end

% Separate the monomials of v: the variables of integration, direction by direction, and
% the other variables, on which the result depends
[nmons,~,nZ] = Zdims(Kop.Z.in);
dpos = [0,cumsum(Kop.depmat.in)];               % dummy variables of direction i: dpos(i)+1:dpos(i+1)
dumnames = [cell(1,0),varnames_in{:}];
np = n*p;
Cv = full(v.C);                                 % monomial in each row, entry (j,l) of v in column (l-1)*n+j
nmv = size(Cv,1);
if nmv==0 || ~any(Cv(:))
    % v is zero
    c = declare_result(zeros(m*(nd+1),p),Kop.Zd,m,p);
    return
end
vdeg = full(v.degmat);
[is_dum,loc] = ismember(v.varname(:),dumnames(:));
degD = zeros(nmv,numel(dumnames));
degD(:,loc(is_dum)) = vdeg(:,is_dum);
namesF = v.varname(~is_dum);
if isempty(namesF)
    Fdeg = zeros(1,0);      iF = ones(nmv,1);
else
    [Fdeg,~,iF] = unique(vdeg(:,~is_dum),'rows');
end
nF = size(Fdeg,1);
% Distinct monomials of v in the dummy variables of each direction
Mdeg = cell(1,N);       iM = zeros(nmv,N);      nmv_i = zeros(1,N);
for i=1:N
    [Mdeg{i},~,iM(:,i)] = unique(degD(:,dpos(i)+1:dpos(i+1)),'rows');
    nmv_i(i) = size(Mdeg{i},1);
end
% Coefficients of v in a tensor, with a dimension for the monomials in each direction, and one for
% the monomials in the other variables, the entries (f,j,l) in the last dimension, f fastest
szT = [nmv_i,nF];
idx = num2cell([iM,iF],1);
lin = sub2ind(szT,idx{:});
T = zeros(prod(szT),np);
for cc=1:np
    T(:,cc) = accumarray(lin,Cv(:,cc),[prod(szT),1]);
end
dimsT = [nmv_i,nF*np];
T = reshape(T,dimsT);

% Integral of the product of each monomial in the kernels and each monomial in v, along each
% direction, for each of the terms
int_order = Kop.int_order;
A = cell(1,N);
for i=1:N
    Zi = Kop.Z.in{i};       Mi = Mdeg{i};
    nzi = size(Zi,1);       nmi = size(Mi,1);
    % The exponents of all products, with the monomial in Zi changing fastest
    S = Zi(repmat((1:nzi)',nmi,1),:) + Mi(repelem((1:nmi)',nzi,1),:);
    [Eu,~,iE] = unique(S,'rows');
    nE = size(Eu,1);
    pE = polynomial(speye(nE),Eu,varnames_in{i}(:),[nE,1]);
    ntrms_i = size(Kop.intmat{i},1);
    A{i} = cell(ntrms_i,1);
    for j=1:ntrms_i
        F = integrate_term(pE,Kop,int_order,pvars,i,j);
        A{i}{j} = reshape(F(iE),nzi,nmi);
    end
end

% Contract with the integrals along each direction, and combine with the coefficients
Cres = zeros(m*(nd+1),nF*p);
sz_params = [size(Kop.params,1:N),1];
Cres = apply_regions(Cres,T,dimsT,1,cell(1,0),Kop,A,sz_params,nmons,nZ,n,p,nF);

% Declare the result
if isempty(namesF)
    c = declare_result(Cres,Kop.Zd,m,p);
elseif nd==0
    % Polynomial in the other variables: monomials as rows, entry (i,l) in column (l-1)*m+i
    [ii,jj,vv] = find(Cres);
    ii = ii(:);     jj = jj(:);     vv = vv(:);
    f = mod(jj-1,nF)+1;     l = ceil(jj/nF);
    c = polynomial(sparse(f,(l-1)*m+ii,vv,nF,m*p),Fdeg,namesF(:),[m,p]);
else
    c = dpvar(sparse(Cres),Fdeg,namesF(:),Kop.Zd,[m,p]);
end

end


function Cres = apply_regions(Cres,T,dimsT,i,jidx,Kop,A,sz_params,nmons,nZ,n,p,nF)
% Contract the tensor T, with dimensions dimsT =
% [nmv_i,...,nmv_N,nF*n*p,nmons_1,...,nmons_{i-1}], along the first
% dimension with the integrals A{i}{j} for each term j along direction i,
% and add the product of the coefficients of each of the regions with the
% integrals W_r, once all directions are integrated
N = numel(Kop.vars.in);
if i>N
    % T has dimensions [nF*n*p,nmons_1,...,nmons_N]. Order the monomials with the first variable
    % changing slowest, and move the matrix columns j to the rows, W_r has rows (j-1)*nZ+z
    Tz = permute(T,[N+1:-1:2,1]);
    Tz = reshape(Tz,[nZ,nF,n,p]);
    W = reshape(permute(Tz,[1,3,2,4]),[nZ*n,nF*p]);
    k = sub2ind(sz_params,jidx{:});
    Cres = Cres + Kop.params{k}*W;
    return
end
nterms_i = numel(A{i});
Tm = reshape(T,dimsT(1),[]);
for j=1:nterms_i
    % Integrate along direction i for term j, and move the monomials of the kernel to the last dimension
    Tj = A{i}{j}*Tm;
    dims_new = [dimsT(2:end),nmons(i)];
    Tj = permute(reshape(Tj,[nmons(i),dimsT(2:end)]),[2:numel(dims_new),1]);
    Cres = apply_regions(Cres,Tj,dims_new,i+1,[jidx,{j}],Kop,A,sz_params,nmons,nZ,n,p,nF);
end
end


function F = integrate_term(pE,Kop,int_order,pvars,i,j)
% Integrate the polynomial vector pE in the dummy variables of direction i
% over the domain of term j along this direction, returning the vector of
% values
intmat_i = Kop.intmat{i};
di = size(int_order{i},2);
A = intmat_i(:,1:2:end);
B = intmat_i(:,2:2:end);
% Determine the variables
pvars_i = [Kop.dom.in(i,1); pvars{i}; Kop.dom.in(i,2)];
% Determine the order of integration of the variables in term j
[~,order_j] = sort(int_order{i}(j,:));
Zj = pE;
for k=1:di
    % Determine the interval of integration for variable si_k
    idx = order_j(k);
    var = pvars_i(idx+1);
    L = pvars_i(A(j,idx)+1);
    U = pvars_i(B(j,idx)+1);
    % Perform the integration
    if isequal(A(j,idx),B(j,idx)) && ~ismember(A(j,idx),[0;di+1])   % currently, prohibit boundary evaluation
        % Multiplier operator (incl. evaluation at either boundary,
        % per chi_{[x,y]}(z):=delta(x,z) for x==y)
        Zj = subs(Zj,var,U);
    else
        % Integral operator
        Zj = int_simple(Zj,var,L,U);
    end
end
F = reshape(double_poly(Zj),[],1);
end


function W = double_poly(Zv)
% Convert a polynomial without variables to a double array
if any(any(Zv.degmat~=0))
    error("The integrated polynomial still depends on spatial variables.")
end
sz = size(Zv);
W = reshape(full(sum(Zv.C,1)),sz(1),sz(2));
end


function c = declare_result(Cres,Zd,m,p)
% Declare the result as a polynomial or dpvar
nd = numel(Zd);
if nd==0
    c = polynomial(Cres);
else
    c = dpvar(sparse(Cres),zeros(1,0),cell(0,1),Zd,[m,p]);
end
end
