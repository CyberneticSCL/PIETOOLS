function varargout = pie_test_utils(cmd,varargin)
% PIE_TEST_UTILS Helpers for the exact, solver-free conversion tests of 1D
% PDEs with periodic and other boundary conditions, on the 'copvar' PIEs
% of the converter.
%
% The PIE state xf is a cell with a block for each component of the PIE,
% as the operators T, A, K are 'copvar' objects with a block row/column for
% each of these components. The helpers evaluate these operators exactly on
% polynomials, block by block.
%
%   xf = pie_test_utils('xf',PIE,PDE,x)       blocks of the PIE state of a
%                                             PDE state x (polynomial column)
%   y  = pie_test_utils('apply',P,v)          P*v, for the cell v of blocks
%   [a,m] = pie_test_utils('aux',PIE)         indices of the finite-dimensional
%                                             (auxiliary) blocks, and of the others
%   xp = pie_test_utils('pde_state',PDE,x,k)  rows of x of the PDE state k
%
% DB, 11/07/2026: Initial coding; replaces the 'opvar' interface of the
%                 tests (PIE.F, PIE.misc.periodic) that the converter no
%                 longer has;

switch cmd
    case 'xf',          [varargout{1:nargout}] = fund_state(varargin{:});
    case 'apply',       [varargout{1:nargout}] = apply_cop(varargin{:});
    case 'aux',         [varargout{1:nargout}] = aux_idcs(varargin{:});
    case 'pde_state',   [varargout{1:nargout}] = pde_state(varargin{:});
    otherwise,          error('Unknown command "%s".',cmd);
end
end

function xp = pde_state(PDE,x,k)
offsets = [0;cumsum(PDE.x_tab(:,2))];
xp = x(offsets(k)+1:offsets(k+1));
end

function [aux,main] = aux_idcs(PIE)
dep = PIE.space.x.depmat;
fin = ~any(dep,2);
aux = find(fin)';       main = find(~fin)';
end

function xf = fund_state(PIE,PDE,x)
% Component b of xf is C_b*D^{diff_b} x_k, for the PDE state k=idx(b)
p2p = PIE.pde2pie.x;
nb = size(p2p.idx,1);
s = PDE.vars(1,1);
xf = cell(nb,1);
x = polynomial(x);
for b=1:nb
    k = p2p.idx(b,1);
    xk = pde_state(PDE,x,k);
    for j=1:p2p.diff(b,1)
        xk = diff(xk,s);
    end
    xf{b} = apply_block(PIE.pde2pie.x.C.C{b,b},xk);
end
end

function y = apply_cop(P,v)
% y{r} = sum_c P.C{r,c}*v{c}
C = P.C;
y = cell(size(C,1),1);
for r=1:size(C,1)
    acc = [];
    for c=1:size(C,2)
        if isempty(C{r,c}) || isempty(v{c})
            continue
        end
        term = apply_block(C{r,c},v{c});
        if isempty(acc)
            acc = term;
        else
            acc = acc+term;
        end
    end
    y{r} = acc;
end
end

function out = apply_block(B,v)
% Evaluate a 'sopvar' block on a polynomial column v, exactly
v = polynomial(reshape(v,numel(v),1));
P = sopvar2opvar(B);
q = polynomial(zeros(P.dim(1,2),1));
f = polynomial(zeros(P.dim(2,2),1));
if P.dim(1,2)>0, q = v; end
if P.dim(2,2)>0, f = v; end
s = P.var1;     sd = P.var2;
if P.dim(1,1)>0
    out = polynomial(zeros(P.dim(1,1),1));
    if P.dim(1,2)>0, out = out+P.P*q; end
    if P.dim(2,2)>0, out = out+int(P.Q1*f,s,P.I(1),P.I(2)); end
else
    out = polynomial(zeros(P.dim(2,1),1));
    if P.dim(1,2)>0, out = out+P.Q2*q; end
    if P.dim(2,2)>0
        out = out+P.R.R0*f+int(P.R.R1*subs(f,s,sd),sd,P.I(1),s) ...
                +int(P.R.R2*subs(f,s,sd),sd,s,P.I(2));
    end
end
end
