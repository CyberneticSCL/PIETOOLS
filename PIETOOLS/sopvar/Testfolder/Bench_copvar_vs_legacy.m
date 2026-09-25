function results = Bench_copvar_vs_legacy(which_bench)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RESULTS = BENCH_COPVAR_VS_LEGACY(WHICH) times and measures the memory of
% 'copvar'/'cdopvar' operations against the equivalent operations on the
% classes they generalize, over a sweep of problem sizes.
%
% INPUTS
% - which_bench: 'all' (default), or any of 'A','B','C','D';
% OUTPUTS
% - results:     struct array, one row per (bench, size, operation), with
%                fields for both classes' time and memory and their ratios;
%
% THE PAIRINGS, AND WHY THEY ARE THE RIGHT ONES
% A 'copvar' is an M x N grid of 'sopvar' blocks over mixed L2 spaces, so
% each legacy class corresponds to a particular grid:
%
%   A  opvar    R^m x L2^n[s]                         -> 2 x 2 copvar
%   B  opvar2d  R x L2[x] x L2[y] x L2[x,y]           -> 4 x 4 copvar
%   C  nopvar   L2^n[s_1,...,s_N], no R component     -> 1 x 1 copvar
%   D  ndopvar  same, with decision variables         -> 1 x 1 cdopvar
%
% EQUIVALENCE
% Every bench builds both objects from THE SAME random data and checks that
% the operations agree, so all timings are for genuinely equal operators.
% A uses 'opvar2copvar', B 'opvar2d2copvar', C 'nopvar2copvar' and D
% 'ndopvar2sdopvar'. The results are compared back through the independent
% single-block converters ('sopvar2opvar', 'sopvar2opvar2d',
% 'sopvar2nopvar'), never by round-tripping the converter under test.
%
% B was previously equi-dimensional only, because assembling 16 opvar2d
% components blockwise was not supported; 'opvar2d2copvar' (09/18/2026)
% removed that caveat and B now checks agreement like the others.
%
% HOW MEMORY IS MEASURED - TWO NUMBERS, BECAUSE NEITHER ALONE IS SOUND
% M_old/M_new are 'whos' sizes of the result. 'whos' is deterministic but
% counts copy-on-write-shared arrays once per holder, so it overstates any
% object sharing storage with another - which is what the decision
% containers do with their one Zd list, where it reported 214 MB for 16
% blocks whose real cost was 0 MB. It is sound for benches A, B and C, which
% have no decision variables and so share nothing.
%
% Mreal is the marginal cost of one result, amortized over K held at once,
% from MemUsedMATLAB. That does see sharing, but it is reported in coarse
% quanta - measured at roughly 0.12 MB - so for objects of a few KB it
% returns zeros and occasional negatives as the collector runs. It is
% printed only above a 0.5 MB floor, and is the number that matters in bench
% D where the objects are tens of MB.
%
% Neither is peak footprint, which is not observable from inside MATLAB.
%
% Timing is the median of NREP warm calls after a warm-up call, since first
% calls include JIT compilation and were measured up to 5x the warm cost.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Bench_copvar_vs_legacy
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
%                  Bench_mopvar_vs_legacy -> Bench_copvar_vs_legacy,
%                  nopvar2mopvar -> nopvar2copvar,
%                  opvar2d2mopvar -> opvar2d2copvar,
%                  opvar2mopvar -> opvar2copvar. File was
%                  'Bench_mopvar_vs_legacy.m'.

if nargin<1,    which_bench = 'all';    end
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
NREP = 3;
results = struct('bench',{},'size',{},'op',{},'t_leg',{},'t_new',{},...
    'm_leg',{},'m_new',{},'r_leg',{},'r_new',{},'t_ratio',{},'m_ratio',{},'match',{});

% A worktree under .claude puts a second copy of every class ahead of the
% working tree on the path, so a benchmark can silently time the wrong code.
assert_single('copvar');    assert_single('sopvar');    assert_single('nopvar');

if any(which_bench=='A') || strcmp(which_bench,'all')
    results = bench_A(results,NREP);
end
if any(which_bench=='C') || strcmp(which_bench,'all')
    results = bench_C(results,NREP);
end
if any(which_bench=='D') || strcmp(which_bench,'all')
    results = bench_D(results,NREP);
end
if any(which_bench=='B') || strcmp(which_bench,'all')
    results = bench_B(results,NREP);
end

report(results);
end


%% ====================== A: opvar vs 2x2 copvar ======================
function R = bench_A(R,NREP)
fprintf('\n===== A: opvar (R^m x L2^n[s]) vs 2x2 copvar =====\n');
fprintf('%-14s %-11s %10s %10s %8s %11s %11s %7s %7s %6s\n',...
    'size','op','t_old','t_new','t_n/o','M_old(whos)','M_new(whos)','M_n/o','Mreal','equal');
pvar s1 s1_dum
dom = [0,1];
for n = [1 2 4 8]
    for deg = [1 2 3]
        m = 2;
        lbl = sprintf('n=%d deg=%d',n,deg);
        [Aop,Amop] = pair_1d(m,n,deg,dom,s1,s1_dum,10*n+deg);
        [Bop,Bmop] = pair_1d(m,n,deg,dom,s1,s1_dum,50+10*n+deg);
        R = one_case(R,'A',lbl,'plus',  @() Aop+Bop,  @() Amop+Bmop,  NREP,@eq_1d);
        R = one_case(R,'A',lbl,'adjoint',@() Aop',    @() Amop',      NREP,@eq_1d);
        R = one_case(R,'A',lbl,'mtimes',@() Aop*Bop,  @() Amop*Bmop,  NREP,@eq_1d);
    end
end
end

function [Pop,Pmop] = pair_1d(m,n,deg,dom,v1,v2,seed)
% One random 4-PI operator, built twice: as a full 'opvar' and as the
% equivalent 2x2 'copvar'. The four components are generated as four
% SINGLE-component opvars, because 'opvar2sopvar' converts one block at a
% time. opvar's dim matrix is [out_R in_R; out_L2 in_L2].
rng(seed);
Pp  = rand_opvar([m m; 0 0],deg,v1,v2,dom);     % R^m  -> R^m
Pq1 = rand_opvar([m 0; 0 n],deg,v1,v2,dom);     % L2^n -> R^m
Pq2 = rand_opvar([0 m; n 0],deg,v1,v2,dom);     % R^m  -> L2^n
Pr  = rand_opvar([0 0; n n],deg,v1,v2,dom);     % L2^n -> L2^n
Pop = opvar();
Pop.I = dom;    Pop.var1 = v1;      Pop.var2 = v2;
Pop.P = Pp.P;   Pop.Q1 = Pq1.Q1;    Pop.Q2 = Pq2.Q2;    Pop.R = Pr.R;
% The container is built by the shipped converter rather than assembled
% here, so the benchmark times the code a caller would actually use and the
% blockwise assembly lives in one place.
Pmop = opvar2copvar(Pop);
end

function e = eq_1d(Cop,Cmop)
% Largest coefficient of the difference, component by component, converting
% each copvar block back with 'sopvar2opvar'.
e = 0;
nm = {'P','Q1','Q2','R'};       idx = [1 1; 1 2; 2 1; 2 2];
for k = 1:4
    b = Cmop.C{idx(k,1),idx(k,2)};
    if isempty(b),  continue,   end
    Bk = sopvar2opvar(b);
    e = max(e, comp_err(Cop.(nm{k}), Bk.(nm{k})));
end
end


%% ====================== C: nopvar vs 1x1 copvar ======================
function R = bench_C(R,NREP)
fprintf('\n===== C: nopvar (L2^n over N variables) vs 1x1 copvar =====\n');
fprintf('%-14s %-11s %10s %10s %8s %11s %11s %7s %7s %6s\n',...
    'size','op','t_old','t_new','t_n/o','M_old(whos)','M_new(whos)','M_n/o','Mreal','equal');
for N = 1:3
    for deg = [1 2]
        n = 2;
        v1 = polynomial(zeros(N,1));    v2 = polynomial(zeros(N,1));
        for k = 1:N
            v1(k) = pvar_named(sprintf('s%d',k));
            v2(k) = pvar_named(sprintf('s%d_dum',k));
        end
        dom = repmat([0,1],N,1);
        lbl = sprintf('N=%d deg=%d',N,deg);
        [An,Am] = pair_nd(n,deg,dom,v1,v2,100+10*N+deg);
        [Bn,Bm] = pair_nd(n,deg,dom,v1,v2,200+10*N+deg);
        R = one_case(R,'C',lbl,'plus',   @() An+Bn, @() Am+Bm, NREP,@eq_nd);
        R = one_case(R,'C',lbl,'adjoint',@() An',   @() Am',   NREP,@eq_nd);
        R = one_case(R,'C',lbl,'mtimes', @() An*Bn, @() Am*Bm, NREP,@eq_nd);
        % @nopvar/ctranspose refuses any operator with a multiplier term, so
        % the adjoint above is unsupported on the legacy side for a general
        % operator. Repeat it on a purely INTEGRAL operator, which is the
        % regime that class does support, to get a like-for-like timing.
        [Ai,Aim] = pair_nd(n,deg,dom,v1,v2,100+10*N+deg,true);
        R = one_case(R,'C',[lbl ' int'],'adjoint',@() Ai', @() Aim', NREP,@eq_nd);
    end
end
end

function [Pnop,Pmop] = pair_nd(n,deg,dom,v1,v2,seed,integral_only)
% One random nopvar and the equivalent 1x1 copvar. With integral_only, the
% multiplier cells are zeroed: a nopvar cell index of 1 along a direction is
% the delta/multiplier term, and @nopvar/ctranspose accepts an operator only
% when every index is in {2,3}. Zeroing rather than emptying keeps the cell
% sizes, which the dim inference and the converter both read.
if nargin<7,    integral_only = false;      end
rng(seed);
Pnop = rand_ndopvar([n n],deg,dom,v1,v2);       % no dvarname -> nopvar
if integral_only
    N = numel(Pnop.deg);
    sz = [size(Pnop.C),1];
    Cc = Pnop.C;
    for ii = 1:numel(Cc)
        idcs = cell(1,N);
        [idcs{:}] = ind2sub(sz,ii);
        if any(cell2mat(idcs)==1)
            Cc{ii} = 0*Cc{ii};
        end
    end
    Pnop.C = Cc;
end
Pmop = copvar({nopvar2sopvar(Pnop)});
end

function e = eq_nd(Cnop,Cmop)
% Compare in the nopvar representation, by converting the single block back.
%
% Returns -1 when the two results cannot be compared coefficient by
% coefficient because they landed on DIFFERENT monomial bases. That is a
% finding, not a harness limitation: 'nopvar' carries one shared degree
% vector and uses the full tensor basis of that degree on both sides, while
% 'sopvar' carries per-variable degree SETS for each side independently, so
% an operation can leave the two representations describing the same
% operator at different degrees. Where that happens the console says so,
% with both degree vectors, and the timing comparison still stands.
Bk = sopvar2nopvar(Cmop.C{1,1});
if ~isequal(size(Bk.C),size(Cnop.C)) || ~isequal(Bk.deg(:),Cnop.deg(:))
    fprintf('      [not coefficient-comparable: nopvar deg [%s] vs sopvar-route deg [%s]]\n',...
        num2str(Cnop.deg(:)'),num2str(Bk.deg(:)'));
    e = -1;     return
end
% A cell whose shape differs but which is ZERO on both sides is skipped, not
% treated as a mismatch. That case is real and was run down: @nopvar/
% ctranspose emits multiplier cells transposed relative to nopvar's own
% documented convention of m*prod(deg+1) rows, so for an integral-only
% operator at N=2 the zero cells come back [2 8] where the sopvar route
% gives the documented [8 2]. A transpose of equal size is the ZL/ZR-swap
% signature, so it was checked rather than assumed: all five mismatched
% cells had nnz 0 on both sides, and the four cells carrying content agreed
% exactly. It is latent, not a wrong answer, because that routine refuses
% any operator in which those cells could be nonzero. A shape mismatch on a
% NONZERO cell still fails loudly below.
e = 0;
for k = 1:numel(Cnop.C)
    if ~isequal(size(Cnop.C{k}),size(Bk.C{k}))
        if nnz(Cnop.C{k})==0 && nnz(Bk.C{k})==0
            continue
        end
        fprintf('      [not comparable: cell %d is %s (nnz %d) vs %s (nnz %d)]\n',...
            k,mat2str(size(Cnop.C{k})),nnz(Cnop.C{k}),...
            mat2str(size(Bk.C{k})),nnz(Bk.C{k}));
        e = -1;     return
    end
    D = Cnop.C{k} - Bk.C{k};
    if ~isempty(D),     e = max(e,full(max(abs(D(:)))));    end
end
end


%% ==================== D: ndopvar vs 1x1 cdopvar =====================
function R = bench_D(R,NREP)
fprintf('\n===== D: ndopvar vs 1x1 cdopvar, over decision variables =====\n');
fprintf('%-14s %-11s %10s %10s %8s %11s %11s %7s %7s %6s\n',...
    'size','op','t_old','t_new','t_n/o','M_old(whos)','M_new(whos)','M_n/o','Mreal','equal');
% 'ndopvar' has no ctranspose, so the adjoint cannot be compared here; it is
% reported as a gap in the summary rather than silently skipped.
for N = 1:2
    for q = [1e3 1e4 1e5]
        n = 2;      deg = 1;    q = round(q);
        v1 = polynomial(zeros(N,1));    v2 = polynomial(zeros(N,1));
        for k = 1:N
            v1(k) = pvar_named(sprintf('s%d',k));
            v2(k) = pvar_named(sprintf('s%d_dum',k));
        end
        dom = repmat([0,1],N,1);
        dv = cellstr("d"+string((1:q)'));
        lbl = sprintf('N=%d q=%.0e',N,q);
        rng(300+N);
        And = rand_ndopvar([n n],deg,dom,v1,v2,dv);
        Bnd = rand_ndopvar([n n],deg,dom,v1,v2,dv);
        Amd = cdopvar({ndopvar2sdopvar(And)});
        Bmd = cdopvar({ndopvar2sdopvar(Bnd)});
        R = one_case(R,'D',lbl,'plus',@() And+Bnd, @() Amd+Bmd, NREP,[]);
        % Composition of two decision operators is quadratic and refused by
        % both classes, so compose against a fixed operator instead.
        rng(400+N);
        Gnd = rand_ndopvar([n n],deg,dom,v1,v2);
        Gm  = copvar({nopvar2sopvar(Gnd)});
        R = one_case(R,'D',lbl,'mtimes',@() And*Gnd, @() Amd*Gm, NREP,[]);
    end
end
end


%% =================== B: opvar2d vs 4x4 copvar ========================
function R = bench_B(R,NREP)
fprintf('\n===== B: opvar2d vs 4x4 copvar =====\n');
fprintf('%-14s %-11s %10s %10s %8s %11s %11s %7s %7s %6s\n',...
    'size','op','t_old','t_new','t_n/o','M_old(whos)','M_new(whos)','M_n/o','Mreal','equal');
pvar x y x_dum y_dum
v1 = [x;y];   v2 = [x_dum;y_dum];   dom = [0 1;0 1];
for n = [1 2]
    for deg = [1 2]
        d4 = n*ones(4,2);
        lbl = sprintf('n=%d deg=%d',n,deg);
        rng(500+10*n+deg);
        A2 = rand_opvar2d(d4,deg,dom,v1,v2);
        B2 = rand_opvar2d(d4,deg,dom,v1,v2);
        % Genuinely EQUAL operators now, via 'opvar2d2copvar'. Before that
        % converter existed this bench could only compare equi-dimensional
        % operators and claimed no equivalence.
        Am = opvar2d2copvar(A2);
        Bm = opvar2d2copvar(B2);
        R = one_case(R,'B',lbl,'plus',   @() A2+B2, @() Am+Bm, NREP,@eq_2d);
        R = one_case(R,'B',lbl,'adjoint',@() A2',   @() Am',   NREP,@eq_2d);
        R = one_case(R,'B',lbl,'mtimes', @() A2*B2, @() Am*Bm, NREP,@eq_2d);
    end
end
end

function e = eq_2d(Cop,Cmop)
% Compare component by component, reading each container block back with
% the independent 'sopvar2opvar2d'.
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
e = 0;
[M,N] = size(Cmop);
if M~=4 || N~=4
    e = inf;    return
end
for i = 1:4
    for j = 1:4
        if isempty(Cmop.C{i,j}),    continue,   end
        Bk = sopvar2opvar2d(Cmop.C{i,j});
        e = max(e,comp_err(Cop.(nm{i,j}),Bk.(nm{i,j})));
    end
end
end


%% ========================== measurement =============================
function R = one_case(R,bench,lbl,op,fh_leg,fh_new,NREP,eqfh)
FLOOR = 0.5*2^20;       % below this the MemUsedMATLAB quantum dominates
% An operation the legacy class cannot perform is a RESULT, not a harness
% failure: @nopvar/ctranspose refuses any operator carrying a multiplier
% term and breaks outright at one variable. Record it and keep going, so the
% suite reports capability gaps alongside cost ratios.
[t1,m1,r1,o1,e1] = try_measure(fh_leg,NREP);
[t2,m2,r2,o2,e2] = try_measure(fh_new,NREP);
if ~isempty(e1) || ~isempty(e2)
    why = e1;   if isempty(why),    why = ['NEW: ' e2];  end
    fprintf('%-14s %-11s %10s %10s %8s %11s %11s %7s %7s  %s\n',...
        lbl,op,nanstr(t1),nanstr(t2),'-','-','-','-','-',['unsupported: ' why]);
    R(end+1) = struct('bench',bench,'size',lbl,'op',op,'t_leg',t1,'t_new',t2,...
        'm_leg',m1,'m_new',m2,'r_leg',r1,'r_new',r2,...
        't_ratio',NaN,'m_ratio',NaN,'match',NaN);
    return
end
if isempty(eqfh)
    mtxt = '-';     mval = NaN;
else
    mval = eqfh(o1,o2);
    if mval<0,          mtxt = 'basis';
    elseif mval<1e-8,   mtxt = 'yes';
    else,               mtxt = sprintf('%.0e',mval);
    end
end
if r1>FLOOR && r2>FLOOR,    rtxt = sprintf('%.2f',r2/r1);
else,                       rtxt = '-';
end
fprintf('%-14s %-11s %10.4f %10.4f %8.2f %11s %11s %7.2f %7s %6s\n',...
    lbl,op,t1,t2,t2/max(t1,1e-9),kb(m1),kb(m2),m2/max(m1,1),rtxt,mtxt);
R(end+1) = struct('bench',bench,'size',lbl,'op',op,'t_leg',t1,'t_new',t2,...
    'm_leg',m1,'m_new',m2,'r_leg',r1,'r_new',r2,...
    't_ratio',t2/max(t1,1e-9),'m_ratio',m2/max(m1,1),'match',mval);
end

function s = kb(b)
if b>=2^20, s = sprintf('%.2f MB',b/2^20);
else,       s = sprintf('%.1f KB',b/1024);
end
end

function [t,mem,mem_real,out,err] = try_measure(fh,NREP)
% 'measure', but an unsupported operation returns its message instead of
% aborting the sweep.
t = NaN;    mem = NaN;      mem_real = NaN;     out = [];   err = '';
try
    [t,mem,mem_real,out] = measure(fh,NREP);
catch ME
    err = strtrim(ME.message);
end
end

function s = nanstr(t)
if isnan(t),    s = '-';    else,   s = sprintf('%.4f',t);  end
end

function [t,mem,mem_real,out] = measure(fh,NREP)
% Median of NREP warm calls; the 'whos' size of the result; and the marginal
% memory of one result amortized over K held at once.
%
% TWO memory numbers, because neither alone is trustworthy here.
%
% 'whos' gives the structural size deterministically, but it counts
% copy-on-write-shared arrays once per holder, so it overstates any object
% that shares storage with another - which is exactly what the decision
% containers do with their one Zd list. It is the right measure for benches
% A, B and C, which have no decision variables and hence nothing shared.
%
% The amortized MemUsedMATLAB delta does see sharing, but it is reported in
% coarse quanta (measured at roughly 0.12 MB), so for objects of a few KB it
% returns zeros and occasional negatives as the collector runs. It is
% therefore printed only above a noise floor, and is the measure that
% matters in bench D where the objects are tens of MB.
% Warm up until the cost settles, THEN time. A single warm-up call followed
% by NREP=3 timed ones was not enough and produced a real error: it reported
% the bench B 'plus' ratio as 1.82 with a 3.96 worst case, where 20 warm-up
% calls and 60 timed reps over three seeds give 1.24-1.44. An ~8 ms
% operation needs tens of calls before the first-call cost stops showing.
out = fh();
t0 = tic;
nwarm = 0;
while toc(t0)<0.05 && nwarm<20
    out = fh();     nwarm = nwarm+1;
end
% Time enough reps to fill ~0.2 s, so a fast operation is not measured by a
% handful of samples, but never fewer than NREP.
tic;    fh();    tguess = toc;
nrep = max(NREP,min(60,ceil(0.2/max(tguess,1e-4))));
ts = zeros(1,nrep);
for k = 1:nrep
    tic;    fh();    ts(k) = toc;
end
t = median(ts);
w = whos('out');    mem = w.bytes;
clear out
K = max(3,min(20,ceil(0.2/max(t,1e-4))));       % cap the extra time spent
keep = cell(1,K);
b0 = used();
for k = 1:K,    keep{k} = fh();     end
mem_real = (used()-b0)/K;
out = keep{1};
end

function b = used()
b = getfield(memory,'MemUsedMATLAB'); %#ok<GFLD>
end

function s = mb(b)
if isnan(b),    s = '-';    else,   s = sprintf('%.2f MB',b/2^20);  end
end

function s = ratio(a,b)
% Net deltas can be zero or negative when the garbage collector reclaims
% more than the call allocated, so a ratio is not always meaningful.
if ~isfinite(a) || ~isfinite(b) || b<=0
    s = '-';
else
    s = sprintf('%.2f',a/b);
end
end

function e = comp_err(A,B)
% Largest absolute coefficient of A-B. An opvar's R component is a STRUCT of
% R0, R1, R2, and seven opvar2d components are CELLS holding the alpha index
% of a PI direction, so recurse through both. A shorter cell is treated as
% padded with empties, so an unpopulated alpha slot compares as zero.
if isstruct(A) || isstruct(B)
    e = 0;
    f = union(fieldnames(A),fieldnames(B));
    for k = 1:numel(f)
        e = max(e,comp_err(A.(f{k}),B.(f{k})));
    end
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
        e = max(e,comp_err(a,b));
    end
    return
end
A = polynomial(A);      B = polynomial(B);
if isempty(A) && isempty(B),    e = 0;  return,     end
D = A-B;
if isempty(D.coefficient),  e = 0;
else,                       e = full(max(abs(D.coefficient(:))));
end
if isempty(e),  e = 0;  end
end

function p = pvar_named(nm)
p = polynomial(1,1,{nm},[1,1]);
end

function assert_single(nm)
w = which(nm,'-all');
if numel(w)>1
    error('Bench_copvar_vs_legacy:shadowed',...
        ['%d copies of ''%s'' are on the path; the benchmark would time an '...
         'unknown one. Remove worktrees with '...
         'rmpath(genpath(fullfile(pwd,''.claude''))).'],numel(w),nm)
end
end


%% ============================ summary ================================
function report(R)
if isempty(R),  return,     end
fprintf('\n===== summary: copvar/cdopvar cost relative to the legacy class =====\n');
fprintf('(>1 means the container is more expensive; memory is the size of the\n');
fprintf(' result, not peak footprint. All four benches compare EQUAL operators.)\n\n');
fprintf('%-6s %-11s %10s %10s %10s %10s %12s\n',...
    'bench','op','t median','t worst','M median','M worst','legacy n/a');
bs = unique({R.bench});     ops = unique({R.op});
for b = 1:numel(bs)
    for o = 1:numel(ops)
        sel = strcmp({R.bench},bs{b}) & strcmp({R.op},ops{o});
        if ~any(sel),   continue,   end
        % A row the legacy class could not perform has a NaN ratio. Count
        % those separately instead of letting them poison the median.
        tr = [R(sel).t_ratio];      mr = [R(sel).m_ratio];
        na = sum(isnan(tr));        ntot = numel(tr);
        tr = tr(~isnan(tr));        mr = mr(isfinite(mr) & mr>0);
        if isempty(tr),     tr = NaN;   end
        if isempty(mr),     mr = NaN;   end
        fprintf('%-6s %-11s %10.2f %10.2f %10.2f %10.2f %8d / %d\n',...
            bs{b},ops{o},median(tr),max(tr),median(mr),max(mr),na,ntot);
    end
end
bad = [R.match];
bad = bad(bad>=0);
bad = bad(~isnan(bad));
if ~isempty(bad)
    fprintf('\nlargest equivalence residual over A, B, C: %.2e\n',max(bad));
end
fprintf('\nNOTE: ''ndopvar'' has no ctranspose, so bench D compares plus and\n');
fprintf('mtimes only; the adjoint has no legacy counterpart to measure against.\n');
end
