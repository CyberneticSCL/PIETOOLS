function Test_copvar_converters()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_CONVERTERS exercises 'rand_copvar', 'opvar2copvar' and
% 'copvar2opvar'.
%
% VERIFICATION APPROACH - WHY THERE IS NO ROUND TRIP HERE
% 'opvar2copvar' and 'copvar2opvar' are mutually inverse, so a round-trip
% test is structurally incapable of detecting a pair of inverse errors: two
% mistakes that cancel return the original object bit for bit while the
% intermediate object is a different operator. That is exactly how the
% 'ndopvar2sdopvar' index error and the 'sopvar' ZL/ZR swap survived
% testing. So each direction is checked against the OPERATOR'S ACTION,
% computed independently:
%
%   y_R  = P*x_R + int_a^b Q1(s) x_L(s) ds
%   y_L(s) = Q2(s) x_R + R0(s) x_L(s)
%            + int_a^s R1(s,t) x_L(t) dt + int_s^b R2(s,t) x_L(t) dt
%
% straight from the 'opvar' components with the 'polynomial' class, against
% the blockwise action of the container evaluated by 'apply_sopvar'. Neither
% side uses the routine under test, and neither uses the other converter.
%
% As a second, weaker check each block is also compared against the
% pre-existing single-block converters 'opvar2sopvar' and 'sopvar2opvar',
% which are independent code with their own tests.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_copvar_converters
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
% Initial coding MMP, 09/18/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  Test_mopvar_converters -> Test_copvar_converters,
%                  mopvar2nopvar -> copvar2nopvar,
%                  mopvar2opvar -> copvar2opvar,
%                  mopvar2opvar2d -> copvar2opvar2d,
%                  nopvar2mopvar -> nopvar2copvar,
%                  opvar2d2mopvar -> opvar2d2copvar,
%                  opvar2mopvar -> opvar2copvar, rand_mopvar -> rand_copvar.
%                  File was 'Test_mopvar_converters.m'.

warning('off','sopvar:noncanonicalMultiplier');
tol = 1e-9;
np = 0;     nf = 0;
pvar s1 s1_dum
dom = [0,1];

fprintf('\n===== Test_copvar_converters =====\n');

%% ---------------------------------------------------------------- 1
fprintf('\n-- 1. rand_copvar --\n');
spaces = struct('out',{{ {}, {'s1'} }},'in',{{ {}, {'s1'} }});
dims   = struct('out',[2;3],'in',[2;3]);
rng(1);
P = rand_copvar(spaces,dims,dom,2,0.7);
info = verify(P);
[np,nf] = chk(np,nf,'verify passes',info.true==1,strjoin(info.flags,' | '));
[np,nf] = chk(np,nf,'grid, spaces and dims as requested',...
    isequal(size(P),[2,2]) && isequal(P.vars,{'s1'}) ...
    && isequal(P.space_out,[false;true]) && isequal(P.space_in,[false;true]) ...
    && isequal(P.dim_out(:),[2;3]) && isequal(P.dim_in(:),[2;3]));
% Zero blocks honoured.
Pz = rand_copvar(spaces,dims,dom,2,0.7,[true false;true true]);
[np,nf] = chk(np,nf,'occ leaves a structurally zero block',...
    isempty(Pz.C{1,2}) && verify(Pz).true==1);
[np,nf] = chkerr(np,nf,'an all-empty row is rejected',...
    'rand_copvar:emptyRow',@() rand_copvar(spaces,dims,dom,2,0.7,[false false;true true]));
% Degrees differ per block, which is the point of using rand_sopvar.
degsets = cellfun(@(b) numel(b.ZL),P.C(2,:),'uni',0);
[np,nf] = chk(np,nf,'blocks carry their own bases',~isempty(degsets));

%% ---------------------------------------------------------------- 2
fprintf('\n-- 2. opvar2copvar, against the operator action --\n');
for n = [1 3]
    for deg = [1 2]
        m = 2;
        Pop = full_opvar(m,n,deg,dom,s1,s1_dum,10*n+deg);
        Pm  = opvar2copvar(Pop);
        [np,nf] = chk(np,nf,sprintf('n=%d deg=%d: verify passes',n,deg),...
            verify(Pm).true==1);
        [np,nf] = chk(np,nf,sprintf('n=%d deg=%d: 2x2 over (R^%d,L2^%d)',n,deg,m,n),...
            isequal(size(Pm),[2,2]) && isequal(Pm.dim_out(:),[m;n]) ...
            && isequal(Pm.dim_in(:),[m;n]));
        e = action_err(Pop,Pm,m,n,dom,s1,s1_dum);
        [np,nf] = chk(np,nf,...
            sprintf('n=%d deg=%d: action agrees, err %.2e',n,deg,e),e<tol);
        % Weaker independent check: each block against opvar2sopvar.
        e2 = blocks_vs_opvar2sopvar(Pop,Pm,m,n);
        [np,nf] = chk(np,nf,...
            sprintf('n=%d deg=%d: blocks match opvar2sopvar, err %.2e',n,deg,e2),e2<tol);
    end
end

%% ---------------------------------------------------------------- 3
fprintf('\n-- 3. copvar2opvar, against the operator action --\n');
for n = [1 3]
    for deg = [1 2]
        m = 2;
        % Built by the RANDOMIZER, not by converting an opvar, so this
        % direction is not handed an object the other direction produced.
        sp = struct('out',{{ {}, {'s1'} }},'in',{{ {}, {'s1'} }});
        dm = struct('out',[m;n],'in',[m;n]);
        rng(700+10*n+deg);
        Pm  = rand_copvar(sp,dm,dom,deg,0.7);
        Pop = copvar2opvar(Pm);
        [np,nf] = chk(np,nf,sprintf('n=%d deg=%d: dim is [%d %d;%d %d]',n,deg,m,m,n,n),...
            isequal(Pop.dim,[m m;n n]));
        e = action_err(Pop,Pm,m,n,dom,s1,s1_dum);
        [np,nf] = chk(np,nf,...
            sprintf('n=%d deg=%d: action agrees, err %.2e',n,deg,e),e<tol);
    end
end

%% ---------------------------------------------------------------- 4
fprintf('\n-- 4. rejection cases --\n');
sp2 = struct('out',{{ {'s1','s2'} }},'in',{{ {'s1','s2'} }});
dm2 = struct('out',2,'in',2);
rng(9);
P2 = rand_copvar(sp2,dm2,dom,1,0.7);
[np,nf] = chkerr(np,nf,'two spatial variables rejected',...
    'copvar2opvar:tooManyVars',@() copvar2opvar(P2));
[np,nf] = chkerr(np,nf,'non-opvar input rejected',...
    'opvar2copvar:badInput',@() opvar2copvar(42));

%% ---------------------------------------------------------------- 5
fprintf('\n-- 5. opvar2d2copvar / copvar2opvar2d --\n');
% Two checks, because neither alone is enough. The PINNING check fixes which
% component lands in which grid cell, read back through the independent
% 'sopvar2opvar2d'. The HOMOMORPHISM check composes in both worlds; it would
% catch an error that the pinning check cannot, but on its own it would
% accept a consistent relabelling of the four spaces, which pinning rejects.
pvar x y x_dum y_dum
v1_2 = [x;y];   v2_2 = [x_dum;y_dum];   dom2 = [0 1;0 1];
nm2 = {'R00','R0x','R0y','R02';
       'Rx0','Rxx','Rxy','Rx2';
       'Ry0','Ryx','Ryy','Ry2';
       'R20','R2x','R2y','R22'};
for nn = [1 2]
    rng(800+nn);
    A2 = rand_opvar2d(nn*ones(4,2),1,dom2,v1_2,v2_2);
    Am = opvar2d2copvar(A2);
    [np,nf] = chk(np,nf,sprintf('n=%d: verify passes, 4x4 grid',nn),...
        verify(Am).true==1 && isequal(size(Am),[4,4]));
    ep = 0;     nbad = {};
    for i = 1:4
        for j = 1:4
            if isempty(Am.C{i,j}),  continue,   end
            back = sopvar2opvar2d(Am.C{i,j});
            ek = cmp_component(A2.(nm2{i,j}),back.(nm2{i,j}));
            if ~(ek<tol),   nbad{end+1} = nm2{i,j};     end %#ok<AGROW>
            ep = max(ep,ek);
        end
    end
    [np,nf] = chk(np,nf,...
        sprintf('n=%d: every block sits at its own component, err %.2e',nn,ep),...
        ep<tol,strjoin(nbad,' '));
    % Homomorphism, forward direction only.
    rng(900+nn);
    B2 = rand_opvar2d(nn*ones(4,2),1,dom2,v1_2,v2_2);
    Bm = opvar2d2copvar(B2);
    Cm  = Am*Bm;
    Cm2 = opvar2d2copvar(A2*B2);
    [np,nf] = chk(np,nf,sprintf('n=%d: conversion commutes with mtimes',nn),...
        grids_match(Cm,Cm2,tol));
    % And the reverse converter, used on both sides.
    C2a = copvar2opvar2d(Am*Bm);
    C2b = copvar2opvar2d(Am)*copvar2opvar2d(Bm);
    eh = 0;
    for k = 1:numel(nm2),   eh = max(eh,cmp_component(C2a.(nm2{k}),C2b.(nm2{k})));  end
    [np,nf] = chk(np,nf,...
        sprintf('n=%d: copvar2opvar2d commutes with mtimes, err %.2e',nn,eh),eh<tol);
end
[np,nf] = chkerr(np,nf,'three variables rejected by copvar2opvar2d',...
    'copvar2opvar2d:tooManyVars',@() copvar2opvar2d(rand_copvar(...
        struct('out',{{ {'s1','s2','s3'} }},'in',{{ {'s1','s2','s3'} }}),...
        struct('out',1,'in',1),dom,1,0.7)));

%% ---------------------------------------------------------------- 6
fprintf('\n-- 6. nopvar2copvar / copvar2nopvar --\n');
for N = 1:2
    vv1 = polynomial(zeros(N,1));   vv2 = polynomial(zeros(N,1));
    for k = 1:N
        vv1(k) = polynomial(1,1,{sprintf('s%d',k)},[1,1]);
        vv2(k) = polynomial(1,1,{sprintf('s%d_dum',k)},[1,1]);
    end
    rng(1000+N);
    Pn = rand_ndopvar([2 2],1,repmat(dom,N,1),vv1,vv2);
    Pm = nopvar2copvar(Pn);
    [np,nf] = chk(np,nf,sprintf('N=%d: 1x1 container, verify passes',N),...
        isequal(size(Pm),[1,1]) && verify(Pm).true==1);
    % Homomorphism against nopvar's own algebra, forward direction only.
    rng(1100+N);
    Qn = rand_ndopvar([2 2],1,repmat(dom,N,1),vv1,vv2);
    Qm = nopvar2copvar(Qn);
    [np,nf] = chk(np,nf,sprintf('N=%d: conversion commutes with plus',N),...
        grids_match(Pm+Qm,nopvar2copvar(Pn+Qn),tol));
end
% Guards.
[np,nf] = chkerr(np,nf,'a 2x2 grid has no nopvar equivalent',...
    'copvar2nopvar:badGrid',@() copvar2nopvar(P));
sp3 = struct('out',{{ {'s1'} }},'in',{{ {'s2'} }});
[np,nf] = chkerr(np,nf,'differing in/out spaces rejected',...
    'copvar2nopvar:spaceMismatch',@() copvar2nopvar(...
        rand_copvar(sp3,struct('out',2,'in',2),dom,1,0.7)));

%% ----------------------------------------------------------------
fprintf('\n===== %d passed, %d FAILED =====\n\n',np,nf);
if nf>0
    error('Test_copvar_converters:failures','%d check(s) failed.',nf);
end
end

function tf = grids_match(A,B,tol)
% Two containers describe the same operator: same grid and spaces, and each
% block pair agrees. Blocks are compared with the class's own 'eq', which
% reconciles differing monomial bases; a raw coefficient compare would trip
% on bases that legitimately differ between the two construction paths.
tf = isequal(size(A),size(B)) && isequal(A.space_out,B.space_out) ...
     && isequal(A.space_in,B.space_in) ...
     && isequal(A.dim_out(:),B.dim_out(:)) && isequal(A.dim_in(:),B.dim_in(:));
if ~tf,     return,     end
for k = 1:numel(A.C)
    ea = isempty(A.C{k});   eb = isempty(B.C{k});
    if ea && eb,    continue,   end
    if ea || eb
        % One side left the block structurally zero. That is only a match if
        % the other side is numerically zero too, which 'eq' cannot be asked
        % about, so treat a populated-vs-absent pair as a mismatch unless
        % the populated one is empty of content.
        nz = A.C{k};    if ea,  nz = B.C{k};    end
        tf = tf && all(cellfun(@(c) isempty(c) || ~any(c(:)),nz.params(:)));
        continue
    end
    try
        tf = tf && all(A.C{k}==B.C{k});
    catch
        tf = false;
    end
    if ~tf,     return,     end
end
end


%% ======================= the independent oracle =======================
function e = action_err(Pop,Pm,m,n,dom,v1,v2)
% Apply both representations to the same test function and compare.
%
% The 'opvar' side is evaluated straight from its components with 'int' and
% 'subs'; the 'copvar' side blockwise with 'apply_sopvar'. Neither path uses
% either converter, so this cannot be fooled by a pair of inverse errors.
xR = testfun_const(m);
xL = testfun_var(n,v1);

% --- opvar action, from the definition ---
yR = Pop.P*xR + int(subs_if(Pop.Q1,v2,v1)*xL,v1,dom(1),dom(2));
yL = Pop.Q2*xR + Pop.R.R0*xL;
xL_t = subs(xL,v1,v2);                      % x_L(t)
yL = yL + int(Pop.R.R1*xL_t,v2,dom(1),v1);  % int_a^s R1(s,t) x(t) dt
yL = yL + int(Pop.R.R2*xL_t,v2,v1,dom(2));  % int_s^b R2(s,t) x(t) dt

% --- copvar action, blockwise ---
zR = zeros_poly(m);     zL = zeros_poly(n);
if ~isempty(Pm.C{1,1}),     zR = zR + apply_sopvar(Pm.C{1,1},xR);   end
if ~isempty(Pm.C{1,2}),     zR = zR + apply_sopvar(Pm.C{1,2},xL);   end
if ~isempty(Pm.C{2,1}),     zL = zL + apply_sopvar(Pm.C{2,1},xR);   end
if ~isempty(Pm.C{2,2}),     zL = zL + apply_sopvar(Pm.C{2,2},xL);   end

e = max(perr(yR,zR),perr(yL,zL));
end

function e = blocks_vs_opvar2sopvar(Pop,Pm,m,n)
% Each container block against the pre-existing single-block converter.
e = 0;
nm  = {'P' ,'Q1'; 'Q2','R'};
sel = {[m m;0 0],[m 0;0 n]; [0 m;n 0],[0 0;n n]};
for i = 1:2
    for j = 1:2
        if isempty(Pm.C{i,j}),  continue,   end
        Pb = opvar();
        Pb.I = Pop.I;   Pb.var1 = Pop.var1;     Pb.var2 = Pop.var2;
        Pb.dim = sel{i,j};
        Pb.(nm{i,j}) = Pop.(nm{i,j});
        ref = opvar2sopvar(Pb);
        got = Pm.C{i,j};
        % Compare through opvar, which normalises the monomial basis; a raw
        % coefficient compare would trip on legitimately different bases.
        a = sopvar2opvar(ref);      b = sopvar2opvar(got);
        e = max(e,cmp_component(a.(nm{i,j}),b.(nm{i,j})));
    end
end
end


%% ============================ helpers ================================
function Pop = full_opvar(m,n,deg,dom,v1,v2,seed)
% A 4-PI opvar with all four components populated, assembled from four
% single-component random opvars.
rng(seed);
Pp  = rand_opvar([m m; 0 0],deg,v1,v2,dom);
Pq1 = rand_opvar([m 0; 0 n],deg,v1,v2,dom);
Pq2 = rand_opvar([0 m; n 0],deg,v1,v2,dom);
Pr  = rand_opvar([0 0; n n],deg,v1,v2,dom);
Pop = opvar();
Pop.I = dom;    Pop.var1 = v1;      Pop.var2 = v2;
Pop.P = Pp.P;   Pop.Q1 = Pq1.Q1;    Pop.Q2 = Pq2.Q2;    Pop.R = Pr.R;
end

function x = testfun_const(n)
x = polynomial((1:n)');
end

function x = testfun_var(n,v)
x = polynomial(zeros(n,1));
for k = 1:n
    x(k) = k + (k+1)*v + 0.5*v^2;
end
end

function z = zeros_poly(n)
z = polynomial(zeros(n,1));
end

function p = subs_if(p,from,to)
% Q1 is a function of the primary variable only, but some generators write
% it in the dummy; substitute only if the dummy actually appears.
p = polynomial(p);
if ismember(from.varname{1},p.varname)
    p = subs(p,from,to);
end
end

function e = cmp_component(A,B)
% An opvar's R is a STRUCT of R0,R1,R2; seven opvar2d components (Rxx, Rx2,
% Ryy, Ry2, R2x, R2y, R22) are CELLS holding the alpha index of a PI
% direction. Recurse through both, and treat a shorter cell as padded with
% empties so an unpopulated alpha slot compares as zero rather than erroring.
if isstruct(A) || isstruct(B)
    e = 0;      f = union(fieldnames(A),fieldnames(B));
    for k = 1:numel(f),     e = max(e,cmp_component(A.(f{k}),B.(f{k})));    end
    return
end
if iscell(A) || iscell(B)
    if ~iscell(A),  A = {A};    end
    if ~iscell(B),  B = {B};    end
    e = 0;
    for k = 1:max(numel(A),numel(B))
        a = [];     b = [];
        if k<=numel(A),     a = A{k};   end
        if k<=numel(B),     b = B{k};   end
        e = max(e,cmp_component(a,b));
    end
    return
end
e = perr(A,B);
end

function e = perr(A,B)
A = polynomial(A);      B = polynomial(B);
if isempty(A) && isempty(B),    e = 0;      return,     end
D = A-B;
if isempty(D.coefficient),  e = 0;
else,                       e = full(max(abs(D.coefficient(:))));
end
if isempty(e),  e = 0;  end
end

function [np,nf] = chk(np,nf,name,tf,extra)
if nargin<5,    extra = '';     end
if tf,  np = np+1;  fprintf('   ok    %s\n',name);
else,   nf = nf+1;  fprintf('   FAIL  %s  %s\n',name,extra);
end
end

function [np,nf] = chkerr(np,nf,name,id,fh)
try
    fh();
    nf = nf+1;      fprintf('   FAIL  %s  (no error raised)\n',name);
catch ME
    if strcmp(ME.identifier,id)
        np = np+1;  fprintf('   ok    %s\n',name);
    else
        nf = nf+1;  fprintf('   FAIL  %s  (got %s: %s)\n',name,ME.identifier,ME.message);
    end
end
end
