function [x,xvarname,dims_in,vars,dom,depmat_in] = rand_polyopvar(p,dims_in,vars,dom,depmat_in)
% [X,XVARNAME,DIMS_IN,VARS,DOM,DEPMAT_IN] = RAND_POLYOPVAR(P,DIMS_IN,VARS,
% DOM,DEPMAT_IN) generates a 'polyopvar' object representing a column
% vector of P (possibly vector-valued) distributed state variables
%   x = [x1;x2;...;xp];
% with x_i of size DIMS_IN(i) (generalizes the scalar-only constructor
% polyopvar(name,vars,dom) to DIMS_IN(i)>1).
%
% INPUTS
% - p:          (optional) scalar integer, number of states. Defaults to
%               randi([1,4]);
% - dims_in:    (optional) p x 1 array of positive integers, vector size
%               of each state. A scalar "m" draws each size as
%               randi([1,m],p,1). Defaults to m=3. Use ones(p,1) for the
%               all-scalar case;
% - vars:       (optional) 1 x N cellstr of spatial variable names, shared
%               by every state. A scalar "N" generates 's1',...,'sN'.
%               Defaults to N=randi([1,3]);
% - dom:        (optional) N x 2 array, domain of each variable (a 1 x 2
%               array is repeated for every variable). Defaults to
%               repmat([0,1],N,1);
% - depmat_in:  (optional) p x N logical array, element (i,j) true if state
%               i depends on variable j. Defaults to true(p,N) (see
%               NOTES);
%
% OUTPUTS
% - x:          sum(DIMS_IN) x 1 'polyopvar' representing the stacked
%               states, each an identity-mapped (degree-1, Z=(1:p)') state
%               of the specified size;
% - xvarname:   p x 1 cellstr, names of the states ('x1',...,'xp');
% - dims_in:    p x 1 array, vector size of each state, as used;
% - vars:       1 x N cellstr, spatial variable names, as used;
% - dom:        N x 2 array, domain of each variable, as used;
% - depmat_in:  p x N logical array, dependency of each state on each
%               variable, as used. The outputs dims_in, vars, dom and
%               depmat_in have the shorthand inputs expanded, for building
%               matching numeric test values, e.g.
%                   xvals{i} = rand_poly([dims_in(i),1],...
%                                       vars(depmat_in(i,:)),deg);
%
% NOTES
% - 'plus.m'/'times.m' do not support combining states with different
%   output (spatial-dependency) spaces (see Test_polyopvar_arith.m), so
%   DEPMAT_IN defaults to all-true, which is safe for '+' and '.*'. Pass a
%   non-default DEPMAT_IN only for vertcat/common_basis-style tests.

% DJ (CC), 10/15/2026: Initial coding;

if nargin<1 || isempty(p)
    p = randi([1,4]);
end
if nargin<2 || isempty(dims_in)
    dims_in = randi([1,3],p,1);
elseif isscalar(dims_in)
    dims_in = randi([1,dims_in],p,1);
elseif numel(dims_in)~=p
    error("Number of state sizes 'dims_in' must match number of states 'p'.")
end
dims_in = dims_in(:);

if nargin<3 || isempty(vars)
    N = randi([1,3]);
    vars = cell(1,N);
    for i=1:N
        vars{i} = ['s',num2str(i)];
    end
elseif isnumeric(vars) && isscalar(vars)
    N = vars;
    vars = cell(1,N);
    for i=1:N
        vars{i} = ['s',num2str(i)];
    end
elseif ~iscellstr(vars)
    error("Spatial variables 'vars' must be specified as 1xN 'cellstr', or scalar N.")
end
vars = reshape(vars,1,[]);
N = numel(vars);

if nargin<4 || isempty(dom)
    dom = repmat([0,1],N,1);
elseif all(size(dom)==[1,2])
    dom = repmat(dom,N,1);
elseif ~all(size(dom)==[N,2])
    error("Domain 'dom' must be specified as Nx2 array for N spatial variables.")
end

if nargin<5 || isempty(depmat_in)
    depmat_in = true(p,N);
elseif ~islogical(depmat_in) || ~all(size(depmat_in)==[p,N])
    error("'depmat_in' must be specified as a px N 'logical' array.")
end

xvarname = cell(p,1);
for i=1:p
    xvarname{i} = ['x',num2str(i)];
end

% Block-diagonal identity operator, as in polyopvar(name,vars,dom) but eye(dims_in(i))
Ccell = cell(p,p);
for i=1:p
    vars_i = struct('in',{vars(depmat_in(i,:))},'out',{vars(depmat_in(i,:))});
    dom_i = struct('in',dom(depmat_in(i,:),:),'out',dom(depmat_in(i,:),:));
    Ci = mat2sopvar(eye(dims_in(i)),vars_i,dom_i);
    Ccell{i,i} = tensopvar(Ci);
end
Cop = tensopmat(Ccell);

x = polyopvar();
x.C = Cop;
x.Z = sparse((1:p)');
x.varname = xvarname;
x.dims.in = dims_in;        x.dims.out = dims_in;
x.vars.in = vars;           x.vars.out = vars;
x.dom.in = dom;             x.dom.out = dom;
x.depmat.in = depmat_in;    x.depmat.out = depmat_in;
x.chkval = true;

end
