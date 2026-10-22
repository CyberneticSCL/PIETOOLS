function Test_monopmat()
% TEST_MONOPMAT tests monopmat, using assertions, for opts.otype(1)==1
% (pointwise tensor-PI product, Parts A-D) and opts.otype(1)==2
% (Kronecker/tensor product, Parts E-F). The case otype(2)==1 is not
% exercised, as monopmat returns an error for it.
%
% Oracle, independent of the index arithmetic in monopmat:
%  - Single-state (q=1, dj=1) cases: enumerate_channels lists every basis
%    "channel" (multiplier/lower/upper integral active along each 3-PI
%    direction, crossed with every output and input monomial) and builds
%    its closed-form response to x=1 from the PI-operator definition with
%    polynomial/int arithmetic. Channels are matched bijectively (exact
%    count, no duplicates, no leftover rows) against apply_sopvar(Pop,1).
%    This is insensitive to row order (immaterial for a basis), but a
%    wrong selection or value fails to match any channel or collides with
%    a claimed one;
%  - Multi-factor (dj>1) and multi-state (q>1) cases check the composition
%    (pointwise product via apply_op, blkdiag pattern) against the same
%    single-state ground truth applied to each piece;
%  - Parts E-F (otype(1)==2): each factor k is enumerated by
%    enumerate_channels on a per-factor spec, its output variables are
%    suffixed '_k' (the apply_op convention for otype(1)==2; all factors
%    of a state share vars_tmp.out, so the suffix is the 1-based factor
%    index), and the factors are combined by cartesian product (matching
%    kron()).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_monopmat
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
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: 3 x 1 'exclude' with several variables is OR, not
%                        AND (oracle changed, N3=2 sweep cases and J1-J2 added);

% NOTES
% - Needs apply_sopvar.m (sopvar/Testfolder) on the path; it is a test
%   utility, not part of monopmat;
% - Channel-level ground truth is derived at nj=1 only. For nj>1 and q>1
%   blocks the test only checks equality with a standalone build at the
%   same nj;
% - Under otype=[2,2] a dj-factor state scales its rows and columns as
%   nj^dj;
% - Domains avoid 0 as an endpoint, so a lower integral of a constant
%   cannot coincide with a monomial at the evaluation point x=1 (false-pass
%   risk).

nfail = 0;
ntotal = 0;

%% ---- Part A: regression cases (ZL/ZR swap, N=0 input variables) ----
[ok,msg] = run_case('A1: ZL/ZR swap (asymmetric degs.in vs degs.out)', @case_A1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('A2: shipped example (N=0 input vars) no longer crashes', @case_A2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part B: single-state (q=1,dj=1) structural sweep ----
% Sweeps N (input var count), N3 (shared/3-PI var count), relative variable
% ordering, exclude patterns, and separability, all at otype(1)=1.
cases = build_sweep_cases();
for i=1:numel(cases)
    cs = cases{i};
    [ok,msg] = run_case(sprintf('B%d: %s',i,cs.label), @() verify_single_state(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part C: multi-factor (dj>1, otype(1)=1 pointwise product) ----
[ok,msg] = run_case('C1: dj=2 pointwise-product full bijection', @case_C1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part D: multi-state (q=2) blkdiag structure ----
[ok,msg] = run_case('D1: q=2 heterogeneous states: blkdiag + diagonal-block match', @case_D1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part E: single-state (q=1,dj>1) Kronecker/tensor product sweep ----
% opts.otype=[2,2]: tensor product of dj>1 factors, mirroring Part B's breadth
% (dj, N3, an extra non-shared input var, exclude/sep, asymmetric degrees).
[ok,msg] = run_case('E0: dj=1, otype(1)=1 vs otype(1)=2 agree (boundary case)', @case_E0);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

tensor_cases = build_tensor_sweep_cases();
for i=1:numel(tensor_cases)
    cs = tensor_cases{i};
    [ok,msg] = run_case(sprintf('E%d: %s',i,cs.label), @() verify_tensor_state(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part F: multi-state (q=2) heterogeneous dj, otype(1)=2 ----
[ok,msg] = run_case('F1: q=2, one dj=2 Kronecker state + one dj=1 state', @case_F1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('F2: q=2, both states dj>1, heterogeneous structure', @case_F2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part G: states mapped to a common output space ----
% Output-only variables for a state, also if shared with other states.
[ok,msg] = run_case('G1: mixed dependence, all states mapped to the same output space', @case_G1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part H: a state depending on the 2nd variable only ----
[ok,msg] = run_case('H1: sep=true, state in s2 only (reported call, hand rows)', @case_H1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);
[ok,msg] = run_case('H2: exclude + sep, state in s2 only (array and cell forms)', @case_H2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);
[ok,msg] = run_case('H3: 3 variables, all dependence sets, sep and exclude', @case_H3);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);
[ok,msg] = run_case('H4: q=2 (s1 | s2), sep=true, own and common output space', @case_H4);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part J: 3 x 1 exclude over several variables (OR) ----
[ok,msg] = run_case('J1: 3 x 1 exclude, states in 1 and 2 of 2 variables (hand rows, array/cell/full forms)', @case_J1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);
[ok,msg] = run_case('J2: 3 x 1 exclude in 2D, kept terms by hand', @case_J2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part I: omitted inputs ----
[ok,msg] = run_case('I1: opts omitted (5 inputs), degs and opts omitted (4)', @case_I1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

fprintf('\n================ Test_monopmat: %d / %d cases passed ================\n', ntotal-nfail, ntotal);
if nfail>0
    error('Test_monopmat: %d / %d cases FAILED.', nfail, ntotal);
else
    fprintf('Test_monopmat: ALL CASES PASSED (otype(1) in {1,2}).\n');
end

end


%% ========================================================================
function [ok,msg] = run_case(name,fn)
% Run one test case, reporting PASS/FAIL without halting the whole suite.
try
    fn();
    ok = true; msg = 'OK';
catch ME
    ok = false; msg = ME.message;
end
fprintf('%-62s %s\n', name, tern(ok,'PASS',['FAIL: ' msg]));
end

function [nfail,ntotal] = tally(nfail,ntotal,ok,~)
ntotal = ntotal+1;
if ~ok, nfail = nfail+1; end
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end


%% ========================================================================
function case_A1()
% ZL/ZR swap regression: equal cardinality (so no dimension-check can catch
% a swap) but DIFFERENT degs.in vs degs.out, so a swap is visible in the
% achievable output degree, not just the row count.
cs = struct();
cs.vars.in={'s1'}; cs.vars.out={'s1'};
cs.dom.in=[0,2]; cs.dom.out=[0,2];
cs.depmat.in=1; cs.depmat.out=1;
cs.ZR={[0;1]}; cs.ZL={[0;3]};
cs.exclude=false(3,1); cs.sep=false;
[matched,total] = verify_single_state(cs);
assert(matched==total && total==10, 'A1: expected 10/10, got %d/%d', matched, total);
end

function case_A2()
% No-crash regression: N=0 input variables (rtn_in=[]), M=1 output variable
% (a ZL/ZR swap fails on numel(ZL)=0 vs numel(vout)=1).
rtn_in = []; rtn_out = 2;
vars_in = {'s1','s2'}; vars_out = {'s1','s2'};
vars = struct(); vars.in = vars_in(rtn_in); vars.out = vars_out(rtn_out);
dom_in = [0,1;-1,0]; dom_out = [0,1;-1,0];
dom = struct(); dom.in = dom_in(rtn_in,:); dom.out = dom_out(rtn_out,:);
N3_idcs = intersect(rtn_in,rtn_out);
depmat_in = [1,1]; depmat_out = [1,1];
depmat = struct(); depmat.in = depmat_in(:,rtn_in); depmat.out = depmat_out(:,rtn_out);
degs_in = {{(0:1)',(0:1)'}};
for i=1:numel(degs_in), degs_in{i} = degs_in{i}(rtn_in); end
degs_out = [2,3]; degs_out = degs_out(rtn_out);
degs = struct(); degs.in = degs_in; degs.out = degs_out;
opts.sep = [true,false]; opts.sep = opts.sep(N3_idcs);
opts.exclude = logical([0,0,1;0,1,0;1,0,0]);
if ~ismember(1,N3_idcs), opts.exclude = all(opts.exclude,2); end
if ~ismember(2,N3_idcs), opts.exclude = all(opts.exclude,1); end
opts.otype = 1;
Zop = monopmat(2,vars,dom,depmat,degs,opts); %#ok<NASGU>
end


%% ========================================================================
function cases = build_sweep_cases()
% Structural sweep for the single-state (q=1,dj=1) verifier. Domains avoid 0 as
% an endpoint: a nonzero lower bound distinguishes a lower integral of a
% constant from a pure monomial at x=1, which could otherwise coincide.
cases = {};

c = struct(); c.label='N=0 (no input vars)';
c.vars.in = {}; c.vars.out = {'y1'};
c.dom.in = zeros(0,2); c.dom.out = [0.5,2.5];
c.depmat.in = zeros(1,0); c.depmat.out = 1;
c.ZR = {}; c.ZL = {[0;1;2]};
c.exclude = false; c.sep = false(1,0);
cases{end+1} = c;

c = struct(); c.label='N=1,N3=0 (full integral)';
c.vars.in = {'x1'}; c.vars.out = {'y1'};
c.dom.in = [0.2,1.3]; c.dom.out = [-0.7,0.9];
c.depmat.in = 1; c.depmat.out = 1;
c.ZR = {[0;1;2]}; c.ZL = {[0;1]};
c.exclude = false; c.sep = false(1,0);
cases{end+1} = c;

c = struct(); c.label='N=1,N3=1, all terms included';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.4,2.1]; c.dom.out = [0.4,2.1];
c.depmat.in = 1; c.depmat.out = 1;
c.ZR = {[0;1]}; c.ZL = {[0;1;2]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

c = struct(); c.label='N=1,N3=1, multiplier excluded';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.4,2.1]; c.dom.out = [0.4,2.1];
c.depmat.in = 1; c.depmat.out = 1;
c.ZR = {[0;1]}; c.ZL = {[0;1;2]};
c.exclude = [true;false;false]; c.sep = false;
cases{end+1} = c;

c = struct(); c.label='N=1,N3=1, separable (upper==lower)';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.4,2.1]; c.dom.out = [0.4,2.1];
c.depmat.in = 1; c.depmat.out = 1;
c.ZR = {[0;1]}; c.ZL = {[0;1;2]};
c.exclude = false(3,1); c.sep = true;
cases{end+1} = c;

c = struct(); c.label='N=2,N3=1,N2=1, 3-PI var listed first';
c.vars.in = {'s1','x1'}; c.vars.out = {'s1'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1];
c.depmat.in = [1,1]; c.depmat.out = 1;
c.ZR = {[0;1],[0;1;2]}; c.ZL = {[0;3]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

c = struct(); c.label='N=2,N3=1,N2=1, 3-PI var listed last';
c.vars.in = {'x1','s1'}; c.vars.out = {'s1'};
c.dom.in = [-1.3,-0.2; 0.4,2.1]; c.dom.out = [0.4,2.1];
c.depmat.in = [1,1]; c.depmat.out = 1;
c.ZR = {[0;1;2],[0;1]}; c.ZL = {[0;3]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

c = struct(); c.label='N=2,N3=2, same relative order';
c.vars.in = {'s1','s2'}; c.vars.out = {'s1','s2'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1; -1.3,-0.2];
c.depmat.in = [1,1]; c.depmat.out = [1,1];
c.ZR = {[0;1],[0;1]}; c.ZL = {[0;2],[0;3]};
c.exclude = false([3,3,1]); c.sep = [false,false];
cases{end+1} = c;

c = struct(); c.label='N=2,N3=2, crossed in/out order, partial exclude';
c.vars.in = {'s1','s2'}; c.vars.out = {'s2','s1'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [-1.3,-0.2; 0.4,2.1];
c.depmat.in = [1,1]; c.depmat.out = [1,1];
c.ZR = {[0;1],[0;1;2]}; c.ZL = {[0;1;2;3],[0;5]};
excl = false([3,3,1]); excl(1,1,1)=true; excl(3,3,1)=true;
c.exclude = excl; c.sep = [false,false];
cases{end+1} = c;

% 3 x 1 exclude with several variables: OR over the variables (10/21/2026)
for e = {[false;false;true],[false;true;false],[true;false;false],[false;true;true]}
    c = struct(); c.label=sprintf('N=2,N3=2, 3 x 1 exclude %s (OR over variables)',mat2str(double(e{1}')));
    c.vars.in = {'s1','s2'}; c.vars.out = {'s1','s2'};
    c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1; -1.3,-0.2];
    c.depmat.in = [1,1]; c.depmat.out = [1,1];
    c.ZR = {[0;1],[0;1]}; c.ZL = {[0;2],[0;3]};
    c.exclude = e{1}; c.sep = [false,false];
    cases{end+1} = c;
end

c = struct(); c.label='N1=1,N2=1,N3=1 mixed';
c.vars.in = {'s1','x1'}; c.vars.out = {'s1','y1'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1; 0.1,0.6];
c.depmat.in = [1,1]; c.depmat.out = [1,1];
c.ZR = {[0;1],[0;1;2]}; c.ZL = {[0;2],[0;1]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

end

function [matched,total] = verify_single_state(cs)
% Check monopmat's basis for this single-state case is bijective with the
% independent channel enumeration below.
vars = struct(); vars.in = cs.vars.in; vars.out = cs.vars.out;
dom = struct(); dom.in = cs.dom.in; dom.out = cs.dom.out;
depmat = struct(); depmat.in = cs.depmat.in; depmat.out = cs.depmat.out;
degs = struct(); degs.in = {cs.ZR}; degs.out = {cs.ZL};
opts = struct(); opts.otype = 1; opts.exclude = cs.exclude; opts.sep = cs.sep;

Zop = monopmat(1,vars,dom,depmat,degs,opts);
Pop = Zop.ops{1,1}.ops{1,1};

[chans,texpect] = enumerate_channels(cs);
nr = size(Pop,1);
assert(nr==texpect, 'row count mismatch: monopmat=%d, independent=%d', nr, texpect);

Px = cleanpoly(apply_sopvar(fill_empty_params(Pop),polynomial(1)),1e-12);
row_used = false(nr,1);
matched = 0;
for c=1:numel(chans)
    found = 0;
    for r=1:nr
        if ~row_used(r)
            d = cleanpoly(Px(r)-chans{c},1e-10);
            if nnz(full(d.coefficient))==0
                found = r; break
            end
        end
    end
    if found>0
        row_used(found) = true;
        matched = matched+1;
    end
end
total = numel(chans);
assert(matched==total, 'only %d / %d channels matched a unique row', matched, total);
assert(sum(row_used)==nr, 'only %d / %d rows accounted for by a channel', sum(row_used), nr);
end

function [chans,total] = enumerate_channels(cs)
% Independent (no monopmat/sopvar internals) enumeration of every basis
% channel's closed-form response to input x=1, from the PI-operator definition:
%   multiplier term:   y(s) += s^zl * x(s)             [x=1 -> just s^zl]
%   lower integral:    y(s) += s^zl * int_a^s t^zr dt
%   upper integral:    y(s) += s^zl * int_s^b t^zr dt
%   N2 (input-only) variable: always fully, definitely integrated.
vars_in = cs.vars.in; vars_out = cs.vars.out;
M = numel(vars_out);
shared = intersect(vars_in,vars_out,'stable');
N3 = numel(shared);
is3_out = ismember(vars_out,shared);
is3_in  = ismember(vars_in,shared);
n2_idx = find(~is3_in);    % input-only (non-shared) positions

if M>0
    pv_out = polynomial(vars_out(:));
else
    pv_out = polynomial(zeros(0,1));
end

shared_out_pos = zeros(1,N3); shared_in_pos = zeros(1,N3);
for d=1:N3
    shared_out_pos(d) = find(strcmp(vars_out,shared{d}));
    shared_in_pos(d)  = find(strcmp(vars_in, shared{d}));
end

% Effective exclude array ([3,...,3,1]), after folding separability as monopmat
% documents: sep(d) sets exclude(..,3,..) true and ORs it into exclude(..,2,..).
if N3==0
    excl = false;
else
    excl = logical(cs.exclude);
    if numel(excl)==3
        % 3 x 1 array e: a term is excluded if its type is excluded along
        % ANY variable (type k_d in variable d: e(k_d) for some d)
        e = excl(:);
        tmp = false([3*ones(1,N3),1]);
        for lin=1:3^N3
            ks = cell(1,N3+1);
            [ks{:}] = ind2sub([3*ones(1,N3),1],lin);
            tmp(lin) = any(e(cell2mat(ks(1:N3))));
        end
        excl = tmp;
    end
    excl = reshape(excl,[3*ones(1,N3),1]);
    sep = cs.sep;
    for d=find(sep)
        idx2 = repmat({1:3},1,N3); idx2{d}=2; idx2=[idx2,{1}];
        idx3 = repmat({1:3},1,N3); idx3{d}=3; idx3=[idx3,{1}];
        excl(idx2{:}) = excl(idx2{:}) | excl(idx3{:});
        excl(idx3{:}) = true;
    end
end

if N3==0
    combos = {[]};
else
    combos = num2cell(idx_product(3*ones(1,N3)),2)';
end

chans = {};
for ci = 1:numel(combos)
    combo = combos{ci};
    if N3>0
        idx_full = num2cell(combo); idx_full = [idx_full,{1}];
        if excl(idx_full{:}), continue; end
    end

    out_counts = cellfun(@numel,cs.ZL);
    if M==0
        out_combos = zeros(1,0);
    else
        out_combos = idx_product(out_counts);
    end

    varying_in_pos = n2_idx;
    varying_in_counts = cellfun(@numel,cs.ZR(n2_idx));
    for d=1:N3
        if combo(d)~=1
            varying_in_pos(end+1) = shared_in_pos(d); %#ok<AGROW>
            varying_in_counts(end+1) = numel(cs.ZR{shared_in_pos(d)}); %#ok<AGROW>
        end
    end
    if isempty(varying_in_pos)
        in_combos = zeros(1,0);
    else
        in_combos = idx_product(varying_in_counts);
    end

    for oi = 1:size(out_combos,1)
        outmon = polynomial(1);
        for v=1:M
            outmon = outmon * pv_out(v)^cs.ZL{v}(out_combos(oi,v));
        end
        for ii = 1:size(in_combos,1)
            expr = outmon;
            for kk=1:numel(n2_idx)
                vpos = n2_idx(kk);
                col = find(varying_in_pos==vpos,1);
                deg = cs.ZR{vpos}(in_combos(ii,col));
                aB = cs.dom.in(vpos,1); bB = cs.dom.in(vpos,2);
                expr = expr * (bB^(deg+1)-aB^(deg+1))/(deg+1);
            end
            % Separable direction: sopvar duplicates the "lower" parameter into
            % the "upper" slot and apply_sopvar sums over all combo slots, so
            % the row picks up int_a^s + int_s^b = int_a^b, independent of s.
            for d=1:N3
                if combo(d)==1
                    continue   % multiplier: dummy pinned to degree 0, factor 1
                end
                vpos = shared_in_pos(d);
                col = find(varying_in_pos==vpos,1);
                deg = cs.ZR{vpos}(in_combos(ii,col));
                sfree = pv_out(shared_out_pos(d));
                aB = cs.dom.out(shared_out_pos(d),1); bB = cs.dom.out(shared_out_pos(d),2);
                if numel(cs.sep)>=d && cs.sep(d)
                    expr = expr * (bB^(deg+1)-aB^(deg+1))/(deg+1);
                elseif combo(d)==2
                    expr = expr * (sfree^(deg+1)-aB^(deg+1))/(deg+1);
                else
                    expr = expr * (bB^(deg+1)-sfree^(deg+1))/(deg+1);
                end
            end
            chans{end+1} = cleanpoly(expr,1e-12); %#ok<AGROW>
        end
    end
end
total = numel(chans);
end

function Pop2 = fill_empty_params(Pop)
% apply_sopvar (test utility in sopvar/Testfolder) does not handle the
% '[]'-as-zero-block shorthand: monopmat leaves an excluded combo's params{k}
% as a 0x0 double, and 'C*trm' then fails on the size mismatch. Expand empty
% params to correctly shaped zero sparse matrices before evaluation.
params = Pop.params;
nL = prod(cellfun(@numel,Pop.ZL));
nR = prod(cellfun(@numel,Pop.ZR));
nrow = Pop.dims(1)*nL; ncol = Pop.dims(2)*nR;
for k=1:numel(params)
    if isempty(params{k})
        params{k} = sparse(nrow,ncol);
    end
end
Pop2 = sopvar(params,Pop.vars,Pop.ZL,Pop.ZR,Pop.dom,Pop.dims);
end

function P = idx_product(counts)
% Returns prod(counts) x numel(counts) matrix of all 1-based index tuples.
n = numel(counts);
if n==0
    P = zeros(1,0);
    return
end
grids = cell(1,n);
[grids{:}] = ndgrid_cellcounts(counts);
P = zeros(numel(grids{1}),n);
for i=1:n
    P(:,i) = grids{i}(:);
end
end

function varargout = ndgrid_cellcounts(counts)
ranges = cellfun(@(c)1:c,num2cell(counts),'UniformOutput',false);
[varargout{1:numel(counts)}] = ndgrid(ranges{:});
end


%% ========================================================================
function case_C1()
% dj=2 (depmat.in=2), otype(1)=1 pointwise product of the 2 factors. The row
% count must equal the product of the per-factor channel counts (10 and 7;
% an index embedding, not a tensor product of values: Y(r)=Px1(r).*Px2(r)),
% and each of the 10x7 (factor1,factor2) channel pairs occurs in exactly one row.
vars = struct(); vars.in = {'s1'}; vars.out = {'s1'};
dom = struct(); dom.in = [0.3,1.7]; dom.out = [0.3,1.7];
depmat = struct(); depmat.in = 2; depmat.out = 1;
degs = struct();
degs.in = {{ [0;1]; [0;1;2] }};
degs.out = [1];
opts = struct(); opts.otype = 1; opts.exclude = {false(3,1)}; opts.sep = {false};

Zop = monopmat(1,vars,dom,depmat,degs,opts);
Zop_tv = Zop.ops{1,1};
Pop1 = Zop_tv.ops{1,1}; Pop2 = Zop_tv.ops{1,2};

s1 = polynomial({'s1'}); s1=s1(1);
a=dom.in(1); b=dom.in(2);
f1 = indep_leaf_channels(s1,a,b,[0;1],[0;1]);
f2 = indep_leaf_channels(s1,a,b,0,[0;1;2]);

X = {polynomial(1), polynomial(1)};
Y = cleanpoly(apply_op(Zop_tv, X),1e-12);
Px1 = cleanpoly(apply_sopvar(Pop1,polynomial(1)),1e-12);
Px2 = cleanpoly(apply_sopvar(Pop2,polynomial(1)),1e-12);
nr = size(Zop_tv,1);
assert(nr==numel(f1)*numel(f2), 'row count mismatch');

maxdiff=0;
for r=1:nr
    d = cleanpoly(Y(r)-Px1(r)*Px2(r),1e-10);
    if ~isempty(d.degmat), maxdiff = max(maxdiff, max(abs(full(d.coefficient)))); end
end
assert(maxdiff<1e-8, 'apply_op pointwise-product sanity check failed');

pair_hit = false(numel(f1),numel(f2));
for r=1:nr
    i1 = find_match_idx(Px1(r), f1);
    i2 = find_match_idx(Px2(r), f2);
    assert(i1>0 && i2>0, 'row %d did not match any independent factor channel', r);
    assert(~pair_hit(i1,i2), 'duplicate pair (%d,%d) at row %d', i1, i2, r);
    pair_hit(i1,i2) = true;
end
assert(all(pair_hit(:)), 'not all factor1xfactor2 pairs realized: %d/%d', sum(pair_hit(:)), numel(pair_hit));
end

function vals = indep_leaf_channels(s1,a,b,ZLdeg,ZRdeg)
vals = {};
for term=1:3
    for zl=1:numel(ZLdeg)
        outmon = s1^ZLdeg(zl);
        if term==1
            vals{end+1} = cleanpoly(outmon,1e-12); %#ok<AGROW>
        else
            for zr=1:numel(ZRdeg)
                m = ZRdeg(zr);
                if term==2, inner=(s1^(m+1)-a^(m+1))/(m+1); else, inner=(b^(m+1)-s1^(m+1))/(m+1); end
                vals{end+1} = cleanpoly(outmon*inner,1e-12); %#ok<AGROW>
            end
        end
    end
end
end

function idx = find_match_idx(p, vals)
idx = 0;
for i=1:numel(vals)
    d = cleanpoly(p-vals{i},1e-10);
    if nnz(full(d.coefficient))==0
        idx = i; return
    end
end
end


%% ========================================================================
function case_D1()
% q=2 states with different structure (nj, N, N3, degrees): blkdiag size, empty
% off-diagonal blocks, and a diagonal block identical (apply_sopvar output) to
% building that state alone.
dims = [2;1];
vars = struct(); vars.in={'s3','x1','x2'}; vars.out={'s3','y1'};
dom = struct(); dom.in=[-1,2; 0.3,1.7; -2,-0.5]; dom.out=[-1,2; 0.2,0.9];
depmat = struct();
depmat.in  = [1,0,0; 0,1,1];
depmat.out = [1,0; 0,1];
degs = struct();
degs.in = { {[0;1]}, {[0;1;2],[0;1]} };
degs.out = { {[0;2]}, {[0;1]} };
opts = struct();
opts.otype = 1;
opts.exclude = {false(3,1), false};
opts.sep = {false, false(1,0)};

Zop = monopmat(dims,vars,dom,depmat,degs,opts);
assert(isequal(size(Zop),[32,3]), 'blkdiag size mismatch: got %s', mat2str(size(Zop)));
assert(isempty(Zop.ops{1,2})||isequal(Zop.ops{1,2},0), 'off-diagonal block (1,2) should be empty');
assert(isempty(Zop.ops{2,1})||isequal(Zop.ops{2,1},0), 'off-diagonal block (2,1) should be empty');

vars1 = struct(); vars1.in = {'s3'}; vars1.out = {'s3'};
dom1 = struct(); dom1.in = [-1,2]; dom1.out = [-1,2];
depmat1 = struct(); depmat1.in = 1; depmat1.out = 1;
degs1 = struct(); degs1.in = {{[0;1]}}; degs1.out = {{[0;2]}};
opts1 = struct(); opts1.otype = 1; opts1.exclude = false(3,1); opts1.sep = false;
Zop1_alone = monopmat(dims(1),vars1,dom1,depmat1,degs1,opts1);
Pop1_alone = Zop1_alone.ops{1,1}.ops{1,1};
Pop1_inblk = Zop.ops{1,1}.ops{1,1};
assert(isequal(size(Pop1_alone),size(Pop1_inblk)), 'diagonal block size mismatch vs standalone build');
Px_alone = cleanpoly(apply_sopvar(Pop1_alone,polynomial(eye(2))),1e-12);
Px_inblk = cleanpoly(apply_sopvar(Pop1_inblk,polynomial(eye(2))),1e-12);
d = cleanpoly(Px_alone-Px_inblk,1e-10);
assert(nnz(full(d.coefficient))==0, 'diagonal block 1 does not match independent single-state build');
end


%% ========================================================================
%  Part E-F support: opts.otype(1)==2 (Kronecker/tensor product, dj>1)
%  ========================================================================
function p2 = rename_output_vars(p,varnames,suffix)
% Mirror the otype(1)==2 renaming of apply_op: every variable of this factor
% in 'varnames' (the state's output vars) is suffixed '_<suffix>'. All factors
% of a monopmat state share vars_tmp.out, so the suffix is the factor index.
p2 = p;
if isempty(p2.varname), return; end
is_v = ismember(p2.varname,varnames);
if any(is_v)
    p2.varname(is_v) = cellfun(@(nm)[nm,'_',num2str(suffix)],p2.varname(is_v),'UniformOutput',false);
end
end

function out = combine_factor_channels(chans_cell)
% Cartesian product of dj per-factor channel lists, multiplying across
% factors (matches the kron() over factors in apply_op).
out = {polynomial(1)};
for k=1:numel(chans_cell)
    new_out = cell(1,numel(out)*numel(chans_cell{k}));
    idx = 0;
    for i=1:numel(out)
        for j=1:numel(chans_cell{k})
            idx = idx+1;
            new_out{idx} = out{i}*chans_cell{k}{j};
        end
    end
    out = new_out;
end
end

function chans = tensor_ground_truth(cs)
% Independent ground truth for a dj>=1, otype(1)==2 state. 'cs' is the spec of
% enumerate_channels, except cs.ZR is dj x N (one row per factor) and cs.ZL is
% 1 x M (factor 1 only; factors 2..dj get the constant monomial, as in
% monopmat's "size(ZLj,1)==1 && dj>1" padding). Each factor is enumerated by
% enumerate_channels, then renamed and combined.
dj = size(cs.ZR,1);
M = numel(cs.vars.out);
chans_per_factor = cell(1,dj);
for k=1:dj
    cs_k = cs;
    cs_k.ZR = cs.ZR(k,:);
    if size(cs.ZL,1)>=k
        cs_k.ZL = cs.ZL(k,:);
    else
        cs_k.ZL = num2cell(zeros(1,M));
    end
    chans_k = enumerate_channels(cs_k);
    chans_k = cellfun(@(p)rename_output_vars(p,cs.vars.out,k),chans_k,'UniformOutput',false);
    chans_per_factor{k} = chans_k;
end
chans = combine_factor_channels(chans_per_factor);
end

function [matched,total] = check_tensor_block(Zop_tv,cs)
% Bijective match of a dj-factor otype(1)==2 tensopvar block against
% 'tensor_ground_truth', via apply_op (kron() of factors).
dj = size(cs.ZR,1);
chans = tensor_ground_truth(cs);
X = repmat({polynomial(1)},1,dj);
% apply_op calls apply_sopvar per factor: fill empty params (see fill_empty_params).
for k=1:size(Zop_tv.ops,2)
    Zop_tv.ops{1,k} = fill_empty_params(Zop_tv.ops{1,k});
end
Y = cleanpoly(apply_op(Zop_tv,X),1e-12);
nr = size(Zop_tv,1);
assert(nr==numel(chans), 'row count mismatch: monopmat=%d, independent=%d', nr, numel(chans));
row_used = false(nr,1);
matched = 0;
for c=1:numel(chans)
    found = 0;
    for r=1:nr
        if ~row_used(r)
            d = cleanpoly(Y(r)-chans{c},1e-10);
            if nnz(full(d.coefficient))==0
                found = r; break
            end
        end
    end
    if found>0
        row_used(found) = true;
        matched = matched+1;
    end
end
total = numel(chans);
assert(sum(row_used)==nr, 'only %d / %d rows accounted for by a channel', sum(row_used), nr);
end

function verify_block_matches_standalone(Zop_inblk,nj,cs,label)
% Check an in-block (q>1) dj-factor tensopvar equals (via apply_op at matching
% nj) a standalone (q=1) build with the same (nj,vars,dom,depmat,degs,exclude,
% sep). Channel-level truth is not re-derived at nj>1 (enumerate_channels
% assumes nj=1); it is established at nj=1 by verify_tensor_state(cs).
vars = struct(); vars.in = cs.vars.in; vars.out = cs.vars.out;
dom = struct(); dom.in = cs.dom.in; dom.out = cs.dom.out;
depmat = struct(); depmat.in = cs.depmat.in; depmat.out = cs.depmat.out;
degs = struct(); degs.in = {cs.ZR}; degs.out = {cs.ZL};
opts = struct(); opts.otype = [2,2]; opts.exclude = cs.exclude; opts.sep = cs.sep;
Zop_alone = monopmat(nj,vars,dom,depmat,degs,opts);
Zop_alone_tv = Zop_alone.ops{1,1};

assert(isequal(size(Zop_inblk),size(Zop_alone_tv)), '%s: size mismatch vs standalone build: %s vs %s', label, mat2str(size(Zop_inblk)), mat2str(size(Zop_alone_tv)));
dj = size(Zop_inblk.ops,2);
for k=1:dj
    Zop_inblk.ops{1,k} = fill_empty_params(Zop_inblk.ops{1,k});
    Zop_alone_tv.ops{1,k} = fill_empty_params(Zop_alone_tv.ops{1,k});
end
Xnj = repmat({polynomial(ones(nj,1))},1,dj);
Y_inblk = cleanpoly(apply_op(Zop_inblk,Xnj),1e-12);
Y_alone = cleanpoly(apply_op(Zop_alone_tv,Xnj),1e-12);
d = cleanpoly(Y_inblk-Y_alone,1e-10);
assert(nnz(full(d.coefficient))==0, '%s: in-block operator does not match standalone build', label);
end

function verify_tensor_state(cs)
% Check monopmat's otype=[2,2] basis for a single dj>1 state against the
% independent Kronecker ground truth.
vars = struct(); vars.in = cs.vars.in; vars.out = cs.vars.out;
dom = struct(); dom.in = cs.dom.in; dom.out = cs.dom.out;
depmat = struct(); depmat.in = cs.depmat.in; depmat.out = cs.depmat.out;
degs = struct(); degs.in = {cs.ZR}; degs.out = {cs.ZL};
opts = struct(); opts.otype = [2,2]; opts.exclude = cs.exclude; opts.sep = cs.sep;

Zop = monopmat(1,vars,dom,depmat,degs,opts);
Zop_tv = Zop.ops{1,1};
[matched,total] = check_tensor_block(Zop_tv,cs);
assert(matched==total, 'only %d / %d channels matched a unique row', matched, total);
end


%% ========================================================================
function case_E0()
% At dj=1, otype=1 and otype=[2,2] take the same code path in monopmat and must
% give identical operators. Reuses Part B's "crossed in/out order, partial
% exclude" configuration.
cs = struct();
cs.vars.in = {'s1','s2'}; cs.vars.out = {'s2','s1'};
cs.dom.in = [0.4,2.1; -1.3,-0.2]; cs.dom.out = [-1.3,-0.2; 0.4,2.1];
cs.depmat.in = [1,1]; cs.depmat.out = [1,1];
cs.ZR = {[0;1],[0;1;2]}; cs.ZL = {[0;1;2;3],[0;5]};
excl = false([3,3,1]); excl(1,1,1)=true; excl(3,3,1)=true;
cs.exclude = excl; cs.sep = [false,false];

vars = struct(); vars.in = cs.vars.in; vars.out = cs.vars.out;
dom = struct(); dom.in = cs.dom.in; dom.out = cs.dom.out;
depmat = struct(); depmat.in = cs.depmat.in; depmat.out = cs.depmat.out;
degs = struct(); degs.in = {cs.ZR}; degs.out = {cs.ZL};

opts1 = struct(); opts1.otype = 1; opts1.exclude = cs.exclude; opts1.sep = cs.sep;
opts2 = struct(); opts2.otype = [2,2]; opts2.exclude = cs.exclude; opts2.sep = cs.sep;
Zop1 = monopmat(1,vars,dom,depmat,degs,opts1);
Zop2 = monopmat(1,vars,dom,depmat,degs,opts2);
Pop1 = Zop1.ops{1,1}.ops{1,1};
Pop2 = Zop2.ops{1,1}.ops{1,1};
assert(isequal(size(Pop1),size(Pop2)), 'dj=1: otype(1)=1 vs 2 size mismatch: %s vs %s', mat2str(size(Pop1)), mat2str(size(Pop2)));
Pop1 = fill_empty_params(Pop1); Pop2 = fill_empty_params(Pop2);
Px1 = cleanpoly(apply_sopvar(Pop1,polynomial(1)),1e-12);
Px2 = cleanpoly(apply_sopvar(Pop2,polynomial(1)),1e-12);
d = cleanpoly(Px1-Px2,1e-10);
assert(nnz(full(d.coefficient))==0, 'dj=1: otype(1)=1 and otype(1)=2 built DIFFERENT operators');
end


%% ========================================================================
function cases = build_tensor_sweep_cases()
% Structural sweep for the single-state (q=1,dj>1) otype=[2,2] verifier.
cases = {};

c = struct(); c.label = 'dj=2, N=1,N3=1, all terms included';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.3,1.7]; c.dom.out = [0.3,1.7];
c.depmat.in = 2; c.depmat.out = 1;
c.ZR = {[0;1]; [0;1;2]}; c.ZL = {[0;1]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

c = struct(); c.label = 'dj=2, N=2 (N3=1,N2=1), heterogeneous per-var degree';
c.vars.in = {'s1','x1'}; c.vars.out = {'s1'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1];
c.depmat.in = [2,1]; c.depmat.out = 1;
c.ZR = {[0;1],[0;1]; [0;1;2],[0]}; c.ZL = {[0;3]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

c = struct(); c.label = 'dj=2, N3=1, multiplier excluded';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.3,1.7]; c.dom.out = [0.3,1.7];
c.depmat.in = 2; c.depmat.out = 1;
c.ZR = {[0;1]; [0;1;2]}; c.ZL = {[0;2]};
c.exclude = [true;false;false]; c.sep = false;
cases{end+1} = c;

c = struct(); c.label = 'dj=2, N3=1, separable (upper==lower)';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.3,1.7]; c.dom.out = [0.3,1.7];
c.depmat.in = 2; c.depmat.out = 1;
c.ZR = {[0;1]; [0;1;2]}; c.ZL = {[0;2]};
c.exclude = false(3,1); c.sep = true;
cases{end+1} = c;

c = struct(); c.label = 'dj=2, N3=2 (two shared vars), asymmetric degrees';
c.vars.in = {'s1','s2'}; c.vars.out = {'s1','s2'};
c.dom.in = [0.4,2.1; -1.3,-0.2]; c.dom.out = [0.4,2.1; -1.3,-0.2];
c.depmat.in = [2,2]; c.depmat.out = [1,1];
c.ZR = {[0;1],[0;1]; [0;1;2],[0;1]}; c.ZL = {[0;2],[0;3]};
c.exclude = false([3,3,1]); c.sep = [false,false];
cases{end+1} = c;

c = struct(); c.label = 'dj=3 (three factors), N=1,N3=1, asymmetric degrees';
c.vars.in = {'s1'}; c.vars.out = {'s1'};
c.dom.in = [0.3,1.7]; c.dom.out = [0.3,1.7];
c.depmat.in = 3; c.depmat.out = 1;
c.ZR = {[0;1]; [0;1;2]; [0;1;2;3]}; c.ZL = {[0;1]};
c.exclude = false(3,1); c.sep = false;
cases{end+1} = c;

end


%% ========================================================================
function case_F1()
% q=2: state1 dj=2 (shared var 's3' + N2 var 'x1'), state2 dj=1 (N2 var 'x2'
% to a disjoint output 'y1'). blkdiag under mixed dj with otype(1)=2, empty
% off-diagonal blocks, each diagonal block equal to a standalone build.
dims = [2;1];
vars = struct(); vars.in={'s3','x1','x2'}; vars.out={'s3','y1'};
dom = struct(); dom.in=[-1,2; 0.3,1.7; -2,-0.5]; dom.out=[-1,2; 0.2,0.9];
depmat = struct();
depmat.in  = [2,1,0; 0,0,1];
depmat.out = [1,0; 0,1];
degs = struct();
degs.in = { {[0;1],[0;1]; [0;1;2],[0]}, {[0;1]} };
degs.out = { {[0;2]}, {[0;1]} };
opts = struct();
opts.otype = [2,2];
opts.exclude = {false(3,1), false};
opts.sep = {false, false(1,0)};

Zop = monopmat(dims,vars,dom,depmat,degs,opts);

cs1 = struct(); cs1.vars.in={'s3','x1'}; cs1.vars.out={'s3'};
cs1.dom.in=[-1,2; 0.3,1.7]; cs1.dom.out=[-1,2];
cs1.depmat.in=[2,1]; cs1.depmat.out=1;
cs1.ZR = {[0;1],[0;1]; [0;1;2],[0]};
cs1.ZL = {[0;2]};
cs1.exclude=false(3,1); cs1.sep=false;

cs2 = struct(); cs2.vars.in={'x2'}; cs2.vars.out={'y1'};
cs2.dom.in=[-2,-0.5]; cs2.dom.out=[0.2,0.9];
cs2.depmat.in=1; cs2.depmat.out=1;
cs2.ZR = {[0;1]};
cs2.ZL = {[0;1]};
cs2.exclude=false; cs2.sep=false(1,0);

% A dj-factor state under otype=[2,2] scales its block as nj^dj, not nj
% (@tensopvar class doc: "(2,2): Rop: L2^{n^2} --> L2^{m^2}"): rows and
% columns both Kronecker-combine across factors.
dj1 = size(cs1.ZR,1); dj2 = size(cs2.ZR,1);
nr1 = numel(tensor_ground_truth(cs1))*dims(1)^dj1;
nr2 = numel(enumerate_channels(cs2))*dims(2)^dj2;
ncol = dims(1)^dj1 + dims(2)^dj2;
assert(isequal(size(Zop),[nr1+nr2,ncol]), 'blkdiag size mismatch: got %s, expected [%d,%d]', mat2str(size(Zop)), nr1+nr2, ncol);
assert(isempty(Zop.ops{1,2})||isequal(Zop.ops{1,2},0), 'off-diagonal block (1,2) should be empty');
assert(isempty(Zop.ops{2,1})||isequal(Zop.ops{2,1},0), 'off-diagonal block (2,1) should be empty');

% Channel-level correctness of each state's own structure (nj=1); cs1's exact
% domains differ from any Part E case.
verify_tensor_state(cs1);
verify_tensor_state(cs2);

% In-block operators at the real nj match standalone builds at the same nj.
verify_block_matches_standalone(Zop.ops{1,1},dims(1),cs1,'state1 (dj=2)');
verify_block_matches_standalone(Zop.ops{2,2},dims(2),cs2,'state2 (dj=1)');
end


%% ========================================================================
function case_F2()
% q=2, both states dj>1 with different dj (2 vs 3), shared variables, exclude
% patterns and per-factor degrees, under otype(1)=2.
dims = [2;1];
vars = struct(); vars.in={'s1','s2'}; vars.out={'s1','s2'};
dom = struct(); dom.in=[0.3,1.7; -1.6,-0.4]; dom.out=[0.3,1.7; -1.6,-0.4];
depmat = struct();
depmat.in  = [2,0; 0,3];
depmat.out = [1,0; 0,1];
degs = struct();
degs.in = { {[0;1]; [0;1;2]}, {[0;1]; [0;1]; [0;1;2;3]} };
degs.out = { {[0;1]}, {[0;1;2]} };
opts = struct();
opts.otype = [2,2];
opts.exclude = {false(3,1), [true;false;false]};
opts.sep = {false, false};

Zop = monopmat(dims,vars,dom,depmat,degs,opts);

cs1 = struct(); cs1.vars.in={'s1'}; cs1.vars.out={'s1'};
cs1.dom.in=[0.3,1.7]; cs1.dom.out=[0.3,1.7];
cs1.depmat.in=2; cs1.depmat.out=1;
cs1.ZR = {[0;1]; [0;1;2]};
cs1.ZL = {[0;1]};
cs1.exclude=false(3,1); cs1.sep=false;

cs2 = struct(); cs2.vars.in={'s2'}; cs2.vars.out={'s2'};
cs2.dom.in=[-1.6,-0.4]; cs2.dom.out=[-1.6,-0.4];
cs2.depmat.in=3; cs2.depmat.out=1;
cs2.ZR = {[0;1]; [0;1]; [0;1;2;3]};
cs2.ZL = {[0;1;2]};
cs2.exclude=[true;false;false]; cs2.sep=false;

% Row/column scaling nj^dj: see case_F1.
dj1 = size(cs1.ZR,1); dj2 = size(cs2.ZR,1);
nr1 = numel(tensor_ground_truth(cs1))*dims(1)^dj1;
nr2 = numel(tensor_ground_truth(cs2))*dims(2)^dj2;
ncol = dims(1)^dj1 + dims(2)^dj2;
assert(isequal(size(Zop),[nr1+nr2,ncol]), 'blkdiag size mismatch: got %s, expected [%d,%d]', mat2str(size(Zop)), nr1+nr2, ncol);
assert(isempty(Zop.ops{1,2})||isequal(Zop.ops{1,2},0), 'off-diagonal block (1,2) should be empty');
assert(isempty(Zop.ops{2,1})||isequal(Zop.ops{2,1},0), 'off-diagonal block (2,1) should be empty');

% Channel-level correctness (nj=1) and consistency at the real nj (see case_F1).
verify_tensor_state(cs1);
verify_tensor_state(cs2);
verify_block_matches_standalone(Zop.ops{1,1},dims(1),cs1,'state1 (dj=2)');
verify_block_matches_standalone(Zop.ops{2,2},dims(2),cs2,'state2 (dj=3)');
end


%% ========================================================================
function case_G1()
% Mixed spatial dependence, all states mapped to the same output space
% (depmat.out all true): a state's output may depend on variables its input
% does not (output-only, constant extension, ZL monomials only). The diagonal
% block of each state must equal the single-state operator built with that
% state's variables alone (checked by verify_single_state), and the row count
% a hand count: each output variable contributes its number of ZL monomials,
% each input variable a further factor 1+2*(dR+1) (multiplier, lower, upper).
vn = {'s1','s2'};
dm = [0.4,2.1; -1.3,-0.2];
dR = 1;         dL = [1,0];
% input variables of each state
cfgs = {{1,[1,2]}, {1,2}, {[],[1,2]}, {[],1}, {2,[1,2]}, {[1,2],[1,2]}};
for c=1:numel(cfgs)
    check_common_space(vn,dm,cfgs{c},dR,dL,c);
end
end

function check_common_space(vn,dm,deps,dR,dL,c)
q = numel(deps);    N = numel(vn);
depin = zeros(q,N);
for i=1:q
    depin(i,deps{i}) = 1;
end
vars = struct('in',{vn},'out',{vn});
dom = struct('in',dm,'out',dm);
depmat = struct('in',depin,'out',true(q,N));
degs = struct('in',dR,'out',dL);
Zop = monopmat(ones(q,1),vars,dom,depmat,degs,struct('otype',1));
assert(isequal(size(Zop.ops),[q,q]),'G1 case %d: wrong block structure.',c);
for i=1:q
    Ni = numel(deps{i});
    nr = prod(dL+1)*(1+2*(dR+1))^Ni;
    assert(Zop.dims.out(i)==nr,'G1 case %d, state %d: %d rows, expected %d.',c,i,Zop.dims.out(i),nr);
    % the same state alone, with its own input variables
    cs = struct();
    cs.vars.in = vn(deps{i});       cs.vars.out = vn;
    cs.dom.in = dm(deps{i},:);      cs.dom.out = dm;
    cs.depmat.in = ones(1,Ni);      cs.depmat.out = ones(1,N);
    cs.ZR = repmat({(0:dR)'},1,Ni);
    cs.ZL = arrayfun(@(d) (0:d)',dL,'UniformOutput',false);
    cs.exclude = false([3*ones(1,Ni),1]);   cs.sep = false(1,Ni);
    [matched,total] = verify_single_state(cs);
    assert(matched==total,'G1 case %d, state %d: standalone operator does not match the channels.',c,i);
    vars_i = struct('in',{cs.vars.in},'out',{cs.vars.out});
    dom_i = struct('in',cs.dom.in,'out',cs.dom.out);
    dep_i = struct('in',cs.depmat.in,'out',cs.depmat.out);
    degs_i = struct('in',{{cs.ZR}},'out',{{cs.ZL}});
    Zi = monopmat(1,vars_i,dom_i,dep_i,degs_i,struct('otype',1,'exclude',cs.exclude,'sep',cs.sep));
    Pi = Zi.ops{1,1}.ops{1,1};
    Pblk = Zop.ops{i,i}.ops{1,1};
    assert(isequal(size(Pi),size(Pblk)),'G1 case %d, state %d: block size differs from standalone build.',c,i);
    Px_i = cleanpoly(apply_sopvar(fill_empty_params(Pi),polynomial(1)),1e-12);
    Px_b = cleanpoly(apply_sopvar(fill_empty_params(Pblk),polynomial(1)),1e-12);
    d = cleanpoly(Px_i-Px_b,1e-10);
    assert(nnz(full(d.coefficient))==0,'G1 case %d, state %d: block differs from standalone build.',c,i);
    % the variables of the tensor-PI operator are stored sorted, so that int accepts it
    vo = Zop.ops{i,i}.vars.out;
    assert(isequal(vo,sort(vo)),'G1 case %d, state %d: output variables are not sorted.',c,i);
    Kop = int(Zop.ops{i,i});
    assert(isequal(size(Kop),size(Zop.ops{i,i})),'G1 case %d, state %d: int gave wrong size.',c,i);
end
end


%% ========================================================================
% Part J: 3 x 1 exclude over several variables (OR)
function case_J1()
% States x1(s1), x2(s1,s2) in the list (s1,s2), dummy degree 1: per variable
% the types give 1 (multiplier), 2 (lower), 2 (upper) rows. A 3 x 1 exclude
% excludes a term if its type is excluded along ANY variable of the state, and
% variables the state does not depend on are ignored: rows by hand
vars = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
dom = struct('in',[0.4,2.1; -1.3,-0.2],'out',[0.4,2.1; -1.3,-0.2]);
depmat = struct('in',[1,0;1,1],'out',[1,1;1,1]);
degs = struct('in',1,'out',0);
excls = {[0;0;1],[0;1;0],[1;0;0],[0;1;1],[1;1;0]};
rows =  {[3,9],   [3,9],   [4,16],  [1,1],   [2,4]};
for k=1:numel(excls)
    e = logical(excls{k});
    % x1: types allowed in s1; x2: product over s1,s2
    Zop = monopmat([1;1],vars,dom,depmat,degs,struct('exclude',e));
    assert(isequal(Zop.dims.out(:)',rows{k}),'exclude=%s: rows %s expected, got %s.',mat2str(e'),mat2str(rows{k}),mat2str(Zop.dims.out(:)'));
    % Same for the cell form with the 3-vector for each state
    Zop_c = monopmat([1;1],vars,dom,depmat,degs,struct('exclude',{{e,e}}));
    assert(isequal(Zop_c.dims.out(:)',rows{k}),'cell exclude=%s: rows %s expected, got %s.',mat2str(e'),mat2str(rows{k}),mat2str(Zop_c.dims.out(:)'));
    % Same as the explicit arrays, built from the allowed types per variable
    al = double(~e);
    Z1 = monopmat([1;1],vars,dom,depmat,degs,struct('exclude',{{~logical(al),~logical(al*al')}}));
    assert(isequal(Z1.dims.out(:)',rows{k}),'full arrays, exclude=%s: rows %s expected, got %s.',mat2str(e'),mat2str(rows{k}),mat2str(Z1.dims.out(:)'));
end
end

function case_J2()
% One state in (s1,s2): the kept terms (rows: type in s1, columns: type in s2;
% 1 multiplier, 2 lower, 3 upper) for a 3 x 1 exclude, by hand
vars = struct('in',{{'s1','s2'}},'out',{{'s1','s2'}});
dom = struct('in',[0.4,2.1; -1.3,-0.2],'out',[0.4,2.1; -1.3,-0.2]);
depmat = struct('in',[1,1],'out',[1,1]);
degs = struct('in',1,'out',0);
excls = {[0;0;1],                  [0;1;0],                  [1;0;0]};
kept  = {logical([1 1 0;1 1 0;0 0 0]), logical([1 0 1;0 0 0;1 0 1]), logical([0 0 0;0 1 1;0 1 1])};
for k=1:numel(excls)
    Zop = monopmat(1,vars,dom,depmat,degs,struct('exclude',excls{k}));
    P = Zop.ops{1,1}.ops{1,1};
    got = ~cellfun(@isempty,P.params);
    assert(isequal(got,kept{k}),'exclude=%s: kept terms differ from the hand-derived pattern.',mat2str(excls{k}'));
end
end

%% ========================================================================
% Part H: state depending on the 2nd variable only
function case_H1()
% sep=true for a state that depends on the 2nd variable of the list only
% (depmat.in=[0 1] for s1,s2), with own and common output space. Rows by hand:
% per output monomial, a multiplier row and one integral of dR+1 rows
% (upper == lower), 1+(dR+1) in all.
vn = {'s1','s2'};   dm = [0.5,2; -1,0.7];
% the call as reported (dR=1, dL=0: 1+2 = 3 rows)
Zop = monopmat(1,struct('in',{vn},'out',{vn}),struct('in',dm,'out',dm),struct('in',[0 1],'out',logical([0 1])),struct('in',1,'out',0),struct('sep',true));
assert(isequal(size(Zop),[3,1]),'H1: reported call gives size %s, expected [3 1].',mat2str(size(Zop)));
for dR=0:2
    for dL=0:1
        for common=[false,true]
            Zop = check_embedded_state(vn,dm,2,dR,dL,false(3,3),[false,true],common,'array');
            Mi = 1+common;          % number of output variables, each with dL+1 monomials
            nr_exp = (dL+1)^Mi*(1+(dR+1));
            assert(nrows(Zop)==nr_exp,'H1: dR=%d, dL=%d, common=%d: %d rows, expected %d.',dR,dL,common,nrows(Zop),nr_exp);
        end
    end
end
end

function case_H2()
% Exclude for a state depending on the 2nd variable only, passed as the 3x3
% array (the state takes slice (1,:), multiplier in s1; other rows are set true
% to catch a wrong slice) and as the cell with its own 3-vector, with and
% without sep. Rows by hand, per output monomial: multiplier 1, lower and upper
% integral dR+1 each; sep excludes the upper integral, and the lower one too
% if the upper was.
vn = {'s1','s2'};   dm = [0.5,2; -1,0.7];
pats = {[0;1;0], false, @(dR) 1+(dR+1);    % lower excluded: multiplier + upper
        [1;0;0], false, @(dR) 2*(dR+1);    % multiplier excluded
        [0;0;1], false, @(dR) 1+(dR+1);    % upper excluded
        [1;0;1], false, @(dR) dR+1;
        [0;1;1], false, @(dR) 1;
        [0;1;0], true,  @(dR) 1;           % sep: lower excluded, so multiplier only
        [0;0;1], true,  @(dR) 1;           % sep: upper excluded forces lower excluded
        [1;0;0], true,  @(dR) dR+1};       % sep: one integral, used for both
for form={'array','cell'}
    for ip=1:size(pats,1)
        e = pats{ip,1};     is_sep = pats{ip,2};
        excl_full = true(3,3);
        excl_full(1,:) = e';
        for dRL=[1,0; 2,1]'
            dR = dRL(1);    dL = dRL(2);
            Zop = check_embedded_state(vn,dm,2,dR,dL,excl_full,[false,is_sep],false,form{1});
            nr_exp = (dL+1)*pats{ip,3}(dR);
            assert(nrows(Zop)==nr_exp,'H2 (%s), exclude %s, sep=%d, dR=%d, dL=%d: %d rows, expected %d.',form{1},mat2str(e'),is_sep,dR,dL,nrows(Zop),nr_exp);
        end
    end
end
end

function case_H3()
% Three variables, every nonempty dependence set, several separability
% patterns, and an asymmetric 3x3x3 exclude array; blocks compared to the same
% state built alone (check_embedded_state).
vn = {'s1','s2','s3'};   dm = [0.5,2; -1,0.7; 0,1];
[A,B,C] = ndgrid(1:3,1:3,1:3);
excl_asym = mod(A+2*B+3*C,4)==0;    % entry (1,1,1) is not excluded: each slice keeps its multiplier
excls = {false(3,3,3), excl_asym};
seps = {[true,true,true], [false,true,false], [false,false,false]};
subsets = {1,2,3,[1 2],[1 3],[2 3],[1 2 3]};
for ie=1:numel(excls)
    for is=1:numel(seps)
        for id=1:numel(subsets)
            deps = subsets{id};
            dR = double(numel(deps)<3);     % keep the 3-variable states small
            check_embedded_state(vn,dm,deps,dR,0,excls{ie},seps{is},false,'array');
        end
    end
end
end

function case_H4()
% Configuration of piesos_poslpivar: one state per variable, depmat.in=[1 0;0 1],
% sep=true, own and common output space. Blocks must equal the single-state
% operators and have the hand-counted sizes.
vn = {'s1','s2'};   dm = [0.5,2; -1,0.7];   dR = 1;     dL = 1;
for common=[false,true]
    depin = [1 0; 0 1];
    depout = depin>0;
    if common, depout = true(2,2); end
    Zop = monopmat([1;1],struct('in',{vn},'out',{vn}),struct('in',dm,'out',dm),struct('in',depin,'out',depout),struct('in',dR,'out',dL),struct('otype',1,'sep',true));
    nr = (dL+1)^(1+common)*(1+(dR+1));
    assert(isequal(size(Zop),[2*nr,2]),'H4: common=%d: size %s, expected [%d 2].',common,mat2str(size(Zop)),2*nr);
    assert(isempty(Zop.ops{1,2})||isequal(Zop.ops{1,2},0),'H4: off-diagonal block (1,2) should be empty');
    assert(isempty(Zop.ops{2,1})||isequal(Zop.ops{2,1},0),'H4: off-diagonal block (2,1) should be empty');
    for i=1:2
        Zi = check_embedded_state(vn,dm,i,dR,dL,false(3,3),[true,true],common,'array');
        assert_same_sopvar(Zi.ops{1,1}.ops{1,1},Zop.ops{i,i}.ops{1,1},sprintf('H4: common=%d, state %d',common,i));
    end
end
end

function Zop = check_embedded_state(vn,dm,deps,dR,dL,excl_full,sep_full,common,form)
% One state (dims 1) depending on vn(deps), built from the full list: excl_full
% is the 3^N exclude array, sep_full the 1 x N separability; form 'array'
% passes the exclude as that array, 'cell' as the cell with the state's own
% array; common=true maps the state to all variables of vn, else to its own.
% The block must equal the same state built alone, which in turn is matched to
% the channel enumeration of verify_single_state.
N = numel(vn);      Ni = numel(deps);
depin = zeros(1,N);     depin(deps) = 1;
if common, depout = true(1,N); else, depout = logical(depin); end
ex_i = slice_exclude(excl_full,deps,N);
sep_i = sep_full(deps);
vars = struct('in',{vn},'out',{vn});
dom = struct('in',dm,'out',dm);
depmat = struct('in',depin,'out',depout);
degs = struct('in',dR,'out',dL);
opts = struct('otype',1,'sep',sep_full);
if strcmp(form,'cell')
    opts.exclude = {ex_i};
else
    opts.exclude = excl_full;
end
Zop = monopmat(1,vars,dom,depmat,degs,opts);
if common, iout = 1:N; else, iout = deps; end
cs = struct();
cs.vars.in = vn(deps);      cs.vars.out = vn(iout);
cs.dom.in = dm(deps,:);     cs.dom.out = dm(iout,:);
cs.depmat.in = ones(1,Ni);  cs.depmat.out = ones(1,numel(iout));
cs.ZR = repmat({(0:dR)'},1,Ni);
cs.ZL = repmat({(0:dL)'},1,numel(iout));
cs.exclude = ex_i;          cs.sep = sep_i;
[matched,total] = verify_single_state(cs);
assert(matched==total,'state in %s: standalone operator does not match the channels.',mat2str(deps));
Zi = monopmat(1,struct('in',{cs.vars.in},'out',{cs.vars.out}),struct('in',cs.dom.in,'out',cs.dom.out), ...
    struct('in',cs.depmat.in,'out',cs.depmat.out),struct('in',{{cs.ZR}},'out',{{cs.ZL}}),struct('otype',1,'exclude',cs.exclude,'sep',cs.sep));
assert_same_sopvar(Zi.ops{1,1}.ops{1,1},Zop.ops{1,1}.ops{1,1},sprintf('state in %s (common=%d, %s exclude)',mat2str(deps),common,form));
end

function ex = slice_exclude(excl_full,deps,N)
% Exclude array of a state that depends on the variables deps (of N): entries
% of excl_full at all combinations of the indices in those variables, index 1
% (multiplier) in all other variables.
Ni = numel(deps);
ex = false([3*ones(1,Ni),1,1]);
combos = idx_product(3*ones(1,Ni));
for r=1:size(combos,1)
    full_idx = num2cell(ones(1,N));
    full_idx(deps) = num2cell(combos(r,:));
    sub = num2cell(combos(r,:));
    ex(sub{:}) = excl_full(full_idx{:});
end
end

function assert_same_sopvar(P1,P2,label)
% Same size, and the same polynomial in each row for input 1.
assert(isequal(size(P1),size(P2)),'%s: sizes differ (%s vs %s).',label,mat2str(size(P1)),mat2str(size(P2)));
Px1 = cleanpoly(apply_sopvar(fill_empty_params(P1),polynomial(1)),1e-12);
Px2 = cleanpoly(apply_sopvar(fill_empty_params(P2),polynomial(1)),1e-12);
d = cleanpoly(Px1-Px2,1e-10);
assert(nnz(full(d.coefficient))==0,'%s: operators differ.',label);
end

function nr = nrows(Zop)
sz = size(Zop);
nr = sz(1);
end


%% ========================================================================
% Part I: omitted inputs
function case_I1()
% opts (6th input) and degs (5th) may be omitted. Defaults: no exclude, no sep,
% degree 0 if degs is omitted too. Rows by hand, one state in both variables,
% each with a multiplier and a lower and upper integral of dR+1 rows, and dL+1
% output monomials: (dL+1)^2*(1+2*(dR+1))^2.
vn = {'s1','s2'};   dm = [0.5,2; -1,0.7];
vars = struct('in',{vn},'out',{vn});
dom = struct('in',dm,'out',dm);
depmat = struct('in',[1 1],'out',true(1,2));
% 5 inputs: same as opts=struct()
degs = struct('in',1,'out',0);
Z0 = monopmat(1,vars,dom,depmat,degs,struct());
Z5 = monopmat(1,vars,dom,depmat,degs);
assert(nrows(Z5)==(0+1)^2*(1+2*(1+1))^2,'I1: 5 inputs: %d rows, expected %d.',nrows(Z5),25);
assert_same_sopvar(Z0.ops{1,1}.ops{1,1},Z5.ops{1,1}.ops{1,1},'I1: opts omitted');
% 4 inputs: degree 0
Z0 = monopmat(1,vars,dom,depmat,struct('in',0,'out',0),struct());
Z4 = monopmat(1,vars,dom,depmat);
assert(nrows(Z4)==(0+1)^2*(1+2*(0+1))^2,'I1: 4 inputs: %d rows, expected %d.',nrows(Z4),9);
assert_same_sopvar(Z0.ops{1,1}.ops{1,1},Z4.ops{1,1}.ops{1,1},'I1: degs and opts omitted');
end
