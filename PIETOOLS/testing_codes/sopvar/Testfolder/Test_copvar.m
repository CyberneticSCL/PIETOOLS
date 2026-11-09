function Test_copvar()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR exercises the 'copvar' container: construction and validation,
% the rejection cases, and the semantics of plus, ctranspose and mtimes.
%
% VERIFICATION APPROACH
% No operation is tested against its own inverse. Each is tested against the
% action of the operator, evaluated independently by 'apply_sopvar', which
% builds the kernels from the sopvar definition and integrates them with the
% 'polynomial' class:
%
%   plus        apply((A+B)_{ij},x) == apply(A_{ij},x) + apply(B_{ij},x)
%   ctranspose  <y, A_{ij} x> == <(A*)_{ji} y, x>, both sides integrated
%   mtimes      apply((A*B)_{ij},x) == sum_k apply(A_{ik}, apply(B_{kj},x))
%
% The adjoint check is the one that can catch a swapped block index, and the
% composition check is the one that can catch a transposed inner loop, so
% the grids below are deliberately NON-SQUARE (2x3 times 3x2) with distinct
% spaces and distinct component counts in every position. A square grid with
% matched dimensions cannot detect either error.
%
% SPACE STRUCTURE
% Spaces are R^n, L_2[s1] and L_2[s1,s2], which is the shape of a PIE: a
% finite-dimensional component alongside L2 components over nested variable
% sets. This exercises n3 = 0, 1 and 2 and the n_i = 0 case the class header
% admits.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_copvar
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
% Initial coding MMP, 09/17/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  Test_mopvar -> Test_copvar. File was 'Test_mopvar.m'.

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
tol = 1e-9;
np = 0;     nf = 0;
rng(0);

% Spaces: S{1} = R, S{2} = L2[s1], S{3} = L2[s1,s2]; all on [0,1].
S = {cell(1,0),{'s1'},{'s1','s2'}};
dm = [0 1];

fprintf('\n===== Test_copvar =====\n');

%% ---------------------------------------------------------------- 1
fprintf('\n-- 1. construction and verify --\n');
% A: 2x3, output spaces (S1,S3), input spaces (S1,S2,S3), asymmetric dims.
out_A = [1 3];          din_A = [1 2 3];
qA = [2 3];             pA = [2 3 4];
CA = mkgrid(S,dm,out_A,din_A,qA,pA,1,true(2,3));
A = copvar(CA);
[np,nf] = chk(np,nf,'size(A) is the block grid',isequal(size(A),[2,3]));
[np,nf] = chk(np,nf,'registry is the sorted union',...
    isequal(A.vars,{'s1','s2'}));
[np,nf] = chk(np,nf,'domains held once',isequal(A.dom,[0 1;0 1]));
[np,nf] = chk(np,nf,'space_out masks',...
    isequal(A.space_out,[false false; true true]));
[np,nf] = chk(np,nf,'space_in masks',...
    isequal(A.space_in,[false false; true false; true true]));
[np,nf] = chk(np,nf,'dim_out',isequal(A.dim_out(:),qA(:)));
[np,nf] = chk(np,nf,'dim_in',isequal(A.dim_in(:),pA(:)));
info = verify(A);
[np,nf] = chk(np,nf,'verify(A) passes',info.true==1,strjoin(info.flags,' | '));

% Zero blocks: metadata must still come out identical.
occ = true(2,3);    occ(1,2) = false;   occ(2,3) = false;
Az = copvar(mkgrid(S,dm,out_A,din_A,qA,pA,1,occ));
info = verify(Az);
[np,nf] = chk(np,nf,'verify passes with zero blocks',info.true==1,...
    strjoin(info.flags,' | '));
[np,nf] = chk(np,nf,'zero-block metadata matches the full grid',...
    isequal(Az.space_out,A.space_out) && isequal(Az.space_in,A.space_in) ...
    && isequal(Az.dim_out,A.dim_out) && isequal(Az.dim_in,A.dim_in));

%% ---------------------------------------------------------------- 2
fprintf('\n-- 2. rejection cases --\n');
[np,nf] = chkerr(np,nf,'conflicting domain for a shared variable',...
    'copvar:domConflict',@() copvar(bad_dom(S,dm)));
[np,nf] = chkerr(np,nf,'row maps into two different spaces',...
    'copvar:rowSpaceMismatch',@() copvar(bad_rowspace(S,dm)));
[np,nf] = chkerr(np,nf,'row has two different output dimensions',...
    'copvar:rowDimMismatch',@() copvar(bad_rowdim(S,dm)));
[np,nf] = chkerr(np,nf,'column maps out of two different spaces',...
    'copvar:colSpaceMismatch',@() copvar(bad_colspace(S,dm)));
[np,nf] = chkerr(np,nf,'entirely empty row',...
    'copvar:emptyRow',@() copvar({mkblk(S,dm,2,2,3,3,1),[];[],[]}));
[np,nf] = chkerr(np,nf,'non-operator block',...
    'copvar:badBlockClass',@() copvar({mkblk(S,dm,2,2,3,3,1),42}));
[np,nf] = chkerr(np,nf,'block carrying a dims matrix',...
    'copvar:blockDimsNotScalar',@() copvar({blkdims_matrix()}));
% A legitimate container whose row orders vars.out differently per column
% must NOT be rejected: this is the canonical-order case that 'isequal' on
% the name lists would have failed.
[np,nf] = chk(np,nf,'row with per-column vars.out order is accepted',...
    verify(copvar(mixed_order(S,dm))).true==1);

%% ---------------------------------------------------------------- 3
fprintf('\n-- 3. plus, against apply_sopvar --\n');
B = copvar(mkgrid(S,dm,out_A,din_A,qA,pA,2,true(2,3)));
AB = A + B;
info = verify(AB);
[np,nf] = chk(np,nf,'verify(A+B) passes',info.true==1,strjoin(info.flags,' | '));
[np,nf] = chk(np,nf,'A+B keeps the grid and spaces',...
    isequal(size(AB),[2,3]) && isequal(AB.space_out,A.space_out) ...
    && isequal(AB.space_in,A.space_in));
err = 0;
for i = 1:2
    for j = 1:3
        x = testfun(S{din_A(j)},pA(j));
        lhs = apply_sopvar(AB.C{i,j},x);
        rhs = apply_sopvar(A.C{i,j},x) + apply_sopvar(B.C{i,j},x);
        err = max(err,polyerr(lhs,rhs));
    end
end
[np,nf] = chk(np,nf,sprintf('apply(A+B) == apply(A)+apply(B), err %.2e',err),err<tol);

% Zero blocks pass through untouched.
AZ = Az + B;
err = 0;
for i = 1:2
    for j = 1:3
        x = testfun(S{din_A(j)},pA(j));
        lhs = apply_sopvar(AZ.C{i,j},x);
        rhs = apply_sopvar(B.C{i,j},x);
        if ~isempty(Az.C{i,j})
            rhs = rhs + apply_sopvar(Az.C{i,j},x);
        end
        err = max(err,polyerr(lhs,rhs));
    end
end
[np,nf] = chk(np,nf,sprintf('zero block acts as zero in plus, err %.2e',err),err<tol);

%% ---------------------------------------------------------------- 4
fprintf('\n-- 4. ctranspose, against the adjoint property --\n');
At = A';
info = verify(At);
[np,nf] = chk(np,nf,'verify(A'') passes',info.true==1,strjoin(info.flags,' | '));
[np,nf] = chk(np,nf,'grid and metadata transpose',...
    isequal(size(At),[3,2]) && isequal(At.space_out,A.space_in) ...
    && isequal(At.space_in,A.space_out) ...
    && isequal(At.dim_out(:),A.dim_in(:)) && isequal(At.dim_in(:),A.dim_out(:)));
% <y_i, A_{ij} x_j> over s^i  ==  <(A*)_{ji} y_i, x_j> over s^j.
err = 0;
for i = 1:2
    for j = 1:3
        x = testfun(S{din_A(j)},pA(j));
        y = testfun(S{out_A(i)},qA(i));
        l = ipair(y, apply_sopvar(A.C{i,j},x), S{out_A(i)}, dm);
        r = ipair(apply_sopvar(At.C{j,i},y), x, S{din_A(j)}, dm);
        err = max(err,abs(l-r));
    end
end
[np,nf] = chk(np,nf,sprintf('<y,Ax> == <A*y,x> on every block, err %.2e',err),err<tol);

%% ---------------------------------------------------------------- 5
fprintf('\n-- 5. mtimes, against composition of actions --\n');
% G: 3x2, input spaces (S1,S2), output spaces = the input spaces of A.
out_G = din_A;          din_G = [1 2];
qG = pA;                pG = [2 3];
G = copvar(mkgrid(S,dm,out_G,din_G,qG,pG,3,true(3,2)));
AG = A*G;
info = verify(AG);
[np,nf] = chk(np,nf,'verify(A*G) passes',info.true==1,strjoin(info.flags,' | '));
[np,nf] = chk(np,nf,'A*G has A''s output and G''s input structure',...
    isequal(size(AG),[2,2]) && isequal(AG.space_out,A.space_out) ...
    && isequal(AG.space_in,G.space_in) ...
    && isequal(AG.dim_out(:),A.dim_out(:)) && isequal(AG.dim_in(:),G.dim_in(:)));
err = 0;
for i = 1:2
    for j = 1:2
        x = testfun(S{din_G(j)},pG(j));
        lhs = apply_sopvar(AG.C{i,j},x);
        rhs = [];
        for k = 1:3
            t = apply_sopvar(A.C{i,k}, apply_sopvar(G.C{k,j},x));
            if isempty(rhs),    rhs = t;    else,   rhs = rhs + t;  end
        end
        err = max(err,polyerr(lhs,rhs));
    end
end
[np,nf] = chk(np,nf,sprintf('apply(A*G) == sum_k apply(A_ik,apply(G_kj)), err %.2e',err),...
    err<tol);
% Dimension guard: A*A must be rejected, the inner spaces do not match.
[np,nf] = chkerr(np,nf,'inner grid mismatch rejected',...
    'mtimes:gridMismatch',@() A*A);
% Scalar factor.
A2 = 2*A;
err = 0;
for i = 1:2
    for j = 1:3
        x = testfun(S{din_A(j)},pA(j));
        err = max(err,polyerr(apply_sopvar(A2.C{i,j},x),2*apply_sopvar(A.C{i,j},x)));
    end
end
[np,nf] = chk(np,nf,sprintf('2*A scales every block, err %.2e',err),err<tol);

% Structural sparsity: a zero block must annihilate its product.
GZ = G;     GZ.C{2,1} = [];
GZ = copvar(GZ.C);
AGZ = A*GZ;
err = 0;
for i = 1:2
    x = testfun(S{din_G(1)},pG(1));
    lhs = apply_sopvar(AGZ.C{i,1},x);
    rhs = [];
    for k = [1 3]
        t = apply_sopvar(A.C{i,k}, apply_sopvar(GZ.C{k,1},x));
        if isempty(rhs),    rhs = t;    else,   rhs = rhs + t;  end
    end
    err = max(err,polyerr(lhs,rhs));
end
[np,nf] = chk(np,nf,sprintf('zero block annihilates its product, err %.2e',err),err<tol);

%% ----------------------------------------------------------------
fprintf('\n===== %d passed, %d FAILED =====\n\n',np,nf);
if nf>0
    error('Test_copvar:failures','%d check(s) failed.',nf);
end
end

%% ======================= helpers =======================
function [np,nf] = chk(np,nf,name,tf,extra)
if nargin<5,    extra = '';     end
if tf
    np = np+1;      fprintf('   ok    %s\n',name);
else
    nf = nf+1;      fprintf('   FAIL  %s  %s\n',name,extra);
end
end

function [np,nf] = chkerr(np,nf,name,id,fh)
% Require the specific identifier, not merely that something threw: a
% different error would mean the guard under test never ran.
try
    fh();
    nf = nf+1;      fprintf('   FAIL  %s  (no error raised)\n',name);
catch ME
    if strcmp(ME.identifier,id)
        np = np+1;  fprintf('   ok    %s\n',name);
    else
        nf = nf+1;
        fprintf('   FAIL  %s  (got %s: %s)\n',name,ME.identifier,ME.message);
    end
end
end

function C = mkgrid(S,dm,out_idx,din_idx,q,p,seed,occ)
% M x N cell of random sopvar blocks over the given output/input spaces.
rng(seed);
M = numel(out_idx);     N = numel(din_idx);
C = cell(M,N);
for i = 1:M
    for j = 1:N
        if ~occ(i,j),   continue,   end
        C{i,j} = mkblk(S,dm,out_idx(i),din_idx(j),q(i),p(j),seed+10*i+j);
    end
end
end

function b = mkblk(S,dm,io,ji,q,p,seed)
rng(seed);
vs = struct('out',{S{io}},'in',{S{ji}});
ds = struct('out',repmat(dm,numel(vs.out),1),'in',repmat(dm,numel(vs.in),1));
gs = struct('out',ones(1,numel(vs.out)),'in',ones(1,numel(vs.in)));
b = rand_sopvar([q,p],vs,ds,gs,0.6);
end

function x = testfun(vars,n)
% n x 1 polynomial test function in the given variables, degree 2 each.
x = polynomial(zeros(n,1));
for k = 1:n
    t = polynomial(k);
    for v = 1:numel(vars)
        s = polynomial(1,1,{vars{v}},[1,1]);
        t = t + (k+v)*s + (v-0.5)*s^2;
    end
    x(k) = t;
end
end

function v = ipair(a,b,vars,dm)
% <a,b> = int a'*b over the domain of the listed variables.
e = a'*b;
for v_ = 1:numel(vars)
    e = int(e,vars{v_},dm(1),dm(2));
end
v = double(e);
end

function e = polyerr(a,b)
d = a-b;
if isa(d,'polynomial')
    if isempty(d.coefficient),  e = 0;  return,     end
    e = full(max(abs(d.coefficient(:))));
else
    e = full(max(abs(d(:))));
end
if isempty(e),   e = 0;     end
end

function C = bad_dom(S,dm)
% Same variable, two domains, in two blocks of one container.
b1 = mkblk(S,dm,2,2,2,2,1);
b2 = mkblk(S,[0 2],2,2,2,2,2);
C = {b1;b2};
end

function C = bad_rowspace(S,dm)
% Row 1 maps into L2[s1] in one column and into L2[s1,s2] in the other.
C = {mkblk(S,dm,2,2,2,2,1), mkblk(S,dm,3,2,2,2,2)};
end

function C = bad_rowdim(S,dm)
C = {mkblk(S,dm,2,2,2,2,1), mkblk(S,dm,2,3,3,2,2)};
end

function C = bad_colspace(S,dm)
% Column 1 maps out of L2[s1] in one row and out of L2[s1,s2] in the other.
C = {mkblk(S,dm,2,2,2,2,1); mkblk(S,dm,2,3,2,2,2)};
end

function C = mixed_order(S,dm)
% Legitimate: row 1 maps into L2[s1,s2] from L2[s1] and from L2[s1,s2]. The
% first block has S2 = {s2}, S3 = {s1}, so vars.out = {'s2','s1'}; the
% second has S2 = {}, S3 = {s1,s2}, so vars.out = {'s1','s2'}. The row is
% consistent as a SET but not as a list.
C = {mkblk(S,dm,3,2,2,2,7), mkblk(S,dm,3,3,2,3,8)};
end

function b = blkdims_matrix()
% A block carrying sopvar's block-dims matrix instead of a 1x2.
v = struct('out',{{'s1'}},'in',{{'s1'}});
d = struct('out',[0 1],'in',[0 1]);
b = sopvar({sparse(10,10),sparse(10,10),sparse(10,10)},v,{[0;1]},{[0;1]},d,[2 2;3 3]);
end
