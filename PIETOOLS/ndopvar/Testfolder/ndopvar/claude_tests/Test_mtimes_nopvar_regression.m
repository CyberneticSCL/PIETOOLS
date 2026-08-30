function Test_mtimes_nopvar_regression()
% TEST_MTIMES_NOPVAR_REGRESSION checks composition of 'nopvar' objects
% against an independent implementation of the definition of the operator,
% provided by APPLY_NOPVAR_QUAD and APPLY_NOPVAR_POLY.
%
% For 'nopvar' objects Qop and Rop and a polynomial function x, the test
% verifies that
%       ((Qop*Rop)*x)(s) = (Qop*(Rop*x))(s)
% at randomly sampled points s, where the left-hand side is evaluated from
% the coefficients returned by MTIMES, and the right-hand side is obtained
% by applying Rop and Qop separately.
%
% Beyond the nominal cases, the following regressions are covered:
% - composition of operators of unequal monomial degree, including degrees
%   that differ only along a trailing spatial direction (MTIMES previously
%   built its permutation matrices from, and stamped the result with, the
%   degree of Qop as supplied rather than the common degree);
% - CHANGE_DEGREE for operators with fewer monomials than spatial
%   directions;
% - composition of operators whose parameters are partly left empty;
% - detection of mismatched spatial domains, including a mismatch confined
%   to the first spatial direction.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_mtimes_nopvar_regression
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
% MMP, 08/23/2026: Initial coding

rng(20260823);
tol = 1e-8;
S = struct('npass',0,'nfail',0);

fprintf('\n=== composition, equal monomial degrees ===\n');
S = compose(S,tol,'1D  d=2,       [0,1]',      1,[2,3],[3,2],2,2,[0 1]);
S = compose(S,tol,'1D  d=3,       [-3,-1]',    1,[2,3],[3,2],3,3,[-3 -1]);
S = compose(S,tol,'1D  d=0,       [0.5,2.5]',  1,[2,3],[3,2],0,0,[0.5 2.5]);
S = compose(S,tol,'1D  d=2,       [-2,0]',     1,[2,3],[3,2],2,2,[-2 0]);
S = compose(S,tol,'2D  d=[2;2]',               2,[2,3],[3,2],[2;2],[2;2],[0 1;-1 1]);
S = compose(S,tol,'2D  d=[0;2],   a1=1',       2,[2,3],[3,2],[0;2],[0;2],[1 3;-2 4]);
S = compose(S,tol,'3D  d=[1;2;1]',             3,[2,3],[3,2],[1;2;1],[1;2;1],[1 3;-2 4;0.5 2]);
S = compose(S,tol,'3D  d=[0;0;0]',             3,[2,2],[2,2],[0;0;0],[0;0;0],[0 1;-1 2;1 4]);

fprintf('\n=== composition, unequal monomial degrees ===\n');
S = compose(S,tol,'1D  d=2 x d=3',             1,[2,3],[3,2],2,3,[0 1]);
S = compose(S,tol,'1D  d=3 x d=2',             1,[2,3],[3,2],3,2,[0 1]);
S = compose(S,tol,'2D  [2;1] x [1;2]',         2,[2,3],[3,2],[2;1],[1;2],[0 1;-1 1]);
S = compose(S,tol,'2D  [1;1] x [2;2]',         2,[2,3],[3,2],[1;1],[2;2],[0 1;-1 1]);
S = compose(S,tol,'2D  [2;1] x [2;3] (trail)', 2,[2,3],[3,2],[2;1],[2;3],[0 1;-1 2]);
S = compose(S,tol,'3D  [1;1;1] x [1;1;2]',     3,[2,2],[2,2],[1;1;1],[1;1;2],[0 1;-1 2;1 3]);
S = compose(S,tol,'3D  [0;0;0] x [1;0;0]',     3,[2,2],[2,2],[0;0;0],[1;0;0],[0 1;0 1;0 1]);
S = compose(S,tol,'4D  [0;0;0;0] x [1;0;0;0]', 4,[2,2],[2,2],[0;0;0;0],[1;0;0;0],[0 1;0 1;0 1;0 1]);

fprintf('\n=== composition with empty (zero) parameters ===\n');
S = compose(S,tol,'1D  Qop.C{2} empty',        1,[2,3],[3,2],2,2,[0 1],{'Q',2});
S = compose(S,tol,'2D  Rop.C{1,3} empty',      2,[2,3],[3,2],[2;1],[2;1],[0 1;-1 1],{'R',7});

fprintf('\n=== randomised sweep, 1D and 2D ===\n');
S = sweep(S,tol,30);

fprintf('\n=== associativity ===\n');
S = assoc(S,tol,'1D',1,[0 1]);
S = assoc(S,tol,'2D',2,[0 1;-1 2]);

fprintf('\n=== error handling ===\n');
S = expect_error(S,'domain mismatch, direction 1',    @() badprod([0 1;0 1],[0 2;0 1]));
S = expect_error(S,'domain mismatch, direction 2',    @() badprod([0 1;0 1],[0 1;0 2]));
S = expect_error(S,'different number of directions',  @() badprod([0 1;0 1],[0 1]));
S = expect_error(S,'double x nopvar, inner dim',      @() badmat('left'));
S = expect_error(S,'nopvar x double, inner dim',      @() badmat('right'));

fprintf('\n%d passed, %d failed\n\n',S.npass,S.nfail);
if S.nfail>0
    error("Test_mtimes_nopvar_regression: %d test(s) failed.",S.nfail)
end

end



%% Compose two random operators and compare against the definition
function S = compose(S,tol,name,N,dimQ,dimR,dQ,dR,dom,mt)

fprintf('  %-32s : ',name);
try
    [Qop,Rop] = rand_pair(N,dimQ,dimR,dQ,dR,dom);
    if nargin>=10
        % Blank out one of the parameters, which should be read as zero
        if strcmp(mt{1},'Q')
            Qop.C{mt{2}} = [];
        else
            Rop.C{mt{2}} = [];
        end
    end
    Pop = Qop*Rop;

    % Check that the composition is a properly formed operator
    if any(isnan(Pop.dim))
        error("Composition returned an operator of inconsistent dimensions.")
    end
    if ~isequal(size(Pop.dom),size(dom)) || any(Pop.dom(:)~=dom(:))
        error("Composition returned an operator on the wrong domain.")
    end

    err = compare(Qop,Rop,Pop,dom,N,dimR(2));
    if err<tol
        fprintf('pass  (err %.2e)\n',err);      S.npass = S.npass+1;
    else
        fprintf('FAIL  (err %.2e)\n',err);      S.nfail = S.nfail+1;
    end
catch ME
    fprintf('FAIL  (%s)\n',ME.message);         S.nfail = S.nfail+1;
end

end



%% Randomised sweep over dimensions, degrees and domains
function S = sweep(S,tol,ntrial)

endpts = [0 1; 1 3; -2 0; -3 -1; 0.5 2.5; -1 1];
nbad = 0;
for it=1:ntrial
    N = 1 + mod(it,2);
    m = 1+randi(2);     p = 1+randi(2);     n = 1+randi(2);
    dQ = randi([0 3],N,1);      dR = randi([0 3],N,1);
    dom = endpts(randi(size(endpts,1),N,1),:);
    try
        [Qop,Rop] = rand_pair(N,[m,p],[p,n],dQ,dR,dom);
        Pop = Qop*Rop;
        err = compare(Qop,Rop,Pop,dom,N,n);
        if err>=tol
            error("mismatch, err = %.2e",err)
        end
    catch ME
        nbad = nbad+1;
        fprintf('  FAIL  N=%d dim=[%d %d %d] dQ=[%s] dR=[%s] dom=[%s] : %s\n', ...
            N,m,p,n,num2str(dQ'),num2str(dR'),num2str(dom(:)'),ME.message);
    end
end
fprintf('  %-32s : %d of %d passed\n','sweep',ntrial-nbad,ntrial);
S.npass = S.npass + (ntrial-nbad);
S.nfail = S.nfail + nbad;

end



%% (Qop*Rop)*Sop should agree with Qop*(Rop*Sop)
function S = assoc(S,tol,name,N,dom)

fprintf('  %-32s : ',name);
try
    [Qop,Rop] = rand_pair(N,[2,2],[2,2],1,1,dom);
    [Sop,~]   = rand_pair(N,[2,2],[2,2],1,1,dom);
    Lop = (Qop*Rop)*Sop;
    Mop = Qop*(Rop*Sop);
    dx = zeros(N,1);        Xc = randn(1,2);
    err = 0;    nrm = 1;
    for trial=1:3
        s = dom(:,1)+(dom(:,2)-dom(:,1)).*rand(N,1);
        yL = apply_nopvar_poly(Lop,Xc,dx,s);
        yM = apply_nopvar_poly(Mop,Xc,dx,s);
        err = max(err,max(abs(yL-yM)));
        nrm = max(nrm,max(abs(yL)));
    end
    err = err/nrm;
    if err<tol
        fprintf('pass  (err %.2e)\n',err);      S.npass = S.npass+1;
    else
        fprintf('FAIL  (err %.2e)\n',err);      S.nfail = S.nfail+1;
    end
catch ME
    fprintf('FAIL  (%s)\n',ME.message);         S.nfail = S.nfail+1;
end

end



%% Check that an invalid composition is rejected
function S = expect_error(S,name,fh)

fprintf('  %-32s : ',name);
try
    fh();
    fprintf('FAIL  (no error raised)\n');       S.nfail = S.nfail+1;
catch
    fprintf('pass  (error raised)\n');          S.npass = S.npass+1;
end

end



%% Products that should be rejected
function badprod(dom1,dom2)

N1 = size(dom1,1);      N2 = size(dom2,1);
[v1,v2] = decl_vars(max(N1,N2));
Qop = rand_nopvar([2,2],ones(N1,1),dom1,v1(1:N1),v2(1:N1));
Rop = rand_nopvar([2,2],ones(N2,1),dom2,v1(1:N2),v2(1:N2));
Qop*Rop;                                                    %#ok<VUNUS>

end

function badmat(side)

[v1,v2] = decl_vars(1);
Rop = rand_nopvar([3,2],1,[0 1],v1,v2);
if strcmp(side,'left')
    randn(2,7)*Rop;                                         %#ok<VUNUS>
else
    Rop*randn(7,2);                                         %#ok<VUNUS>
end

end



%% Compare ((Qop*Rop)*x)(s) against (Qop*(Rop*x))(s)
function relerr = compare(Qop,Rop,Pop,dom,N,nx)

dx = ones(N,1);
Xc = randn(prod(dx+1),nx);
xfun = @(t) apply_nopvar_poly(Rop,Xc,dx,t);
% (Rop*x)(t) has degree at most 2*deg(Rop)+deg(x)+1 along each direction,
% so the integrands passed to APPLY_NOPVAR_QUAD have degree at most
% deg(Qop)+2*deg(Rop)+deg(x)+1, and Gauss--Legendre is exact for those
% provided 2*nq-1 is at least that large
dg = max(Qop.deg) + 2*max(Rop.deg) + max(dx) + 1;
nq = ceil((dg+1)/2) + 2;

err = 0;    nrm = 1;
for trial=1:3
    s = dom(:,1) + (dom(:,2)-dom(:,1)).*rand(N,1);
    y1 = apply_nopvar_quad(Qop,xfun,s,nq);
    y2 = apply_nopvar_poly(Pop,Xc,dx,s);
    err = max(err,max(abs(y1-y2)));
    nrm = max(nrm,max(abs(y1)));
end
relerr = err/nrm;

end



%% Declare a random pair of composable operators
function [Qop,Rop] = rand_pair(N,dimQ,dimR,dQ,dR,dom)

[var1,var2] = decl_vars(N);
Qop = rand_nopvar(dimQ,dQ,dom,var1,var2);
Rop = rand_nopvar(dimR,dR,dom,var1,var2);

end



%% Declare a random 'nopvar' object with dense coefficient matrices
function Pop = rand_nopvar(dim,deg,dom,var1,var2)

N = size(dom,1);
if isscalar(deg)
    deg = deg*ones(N,1);
end
Pop = nopvar();
Pop.deg = deg;
Pop.dom = dom;
Pop.vars = [var1,var2];
Pop.C = cell([3*ones(1,N),1]);

nZ = prod(deg+1);
sz_C = [size(Pop.C),1];
idcs = cell(1,N);
for ii=1:numel(Pop.C)
    [idcs{:}] = ind2sub(sz_C,ii);
    is_int = logical(cell2mat(idcs)-1);
    Pop.C{ii} = sparse(randn(dim(1)*nZ,dim(2)*prod(deg(is_int)+1)));
end

end



%% Declare N primary and N dummy spatial variables
function [var1,var2] = decl_vars(N)

var1_name = [repmat('s',[N,1]),num2str((1:N)')];
var2_name = [var1_name,repmat('_dum',[N,1])];
var1 = polynomial(mat2cell(var1_name,ones(N,1),size(var1_name,2)));
var2 = polynomial(mat2cell(var2_name,ones(N,1),size(var2_name,2)));

end
