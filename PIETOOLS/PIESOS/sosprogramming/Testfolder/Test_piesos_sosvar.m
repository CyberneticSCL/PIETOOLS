function Test_piesos_sosvar()
% TEST_PIESOS_SOSVAR tests, using assertions, piesos_sosvar.
%  Part A: input validation (Z, opts, psatz, states of different
%          dependence, different operators on the left and right).
%  Part B: nZ/Zop/Pmat/monomial-dedup structure for a sweep of bases (one
%          and two states, Z=[x;x.*x], cross term, opts.exclude, opts.sep,
%          degree 3), via a prefix replica ('piesos_sosvar_prefix') and
%          independent monopmat/unique enumeration.
%  Part C: the real piesos_sosvar for the Part B sweep, with concrete
%          decision variable values and state, checked two independent
%          ways:
%          (1) pairwise: the sum over all (ridx,cidx) of the block
%              contributions, rebuilt from Zop/Pmat, against
%              <Zop*Z,Pmat*Zop*Z> by apply_op(Zop,..) and plain
%              polynomial integration (no int()/intop);
%          (2) full V: the returned V (V.Z/V.C) evaluated at the same
%              state.
%  Part D: an unvalidated (chkval=false), inconsistent Z fails fast in
%          initialize.
%  Part E: opts.psatz accepted/rejected, number of SOS variables declared.
%  Part F: V against the sum of squares of Cholesky factors of random PSD
%          Gram matrices (so V>0), and against the quadratic form with
%          Pmat.
%  Part G: the same for states of different dependence (2 and 3
%          variables), psatz, exclude, sep.
%  Part H: SDPs: a psatz multiplier certifies a weight that is positive
%          on the interval only.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_sosvar
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
% DJ (CC), 10/22/2026: Z_by_rows: Zop and Pmat are indexed by the entries
%                        of Z (rows of Z.C), not by the rows of Z.Z;
% DJ (CC), 10/23/2026: Bases from otimes (tensor products), not x1.*x1, which
%                        piesos_quadvar rejects;

% NOTES
% - piesos_sosvar declares V by piesos_quadvar: the loop and kernels are
%   tested in Test_piesos_quadvar. Here only Part D targets the own code
%   of piesos_sosvar (the initialize guard); the Part B prefix replica is
%   an independent check of the structure of Zop and Pmat;
% - The Part C pairwise check loops over the full 1:nZ^2 pairs and does
%   not share int() calls between pairs, so it also checks any symmetry
%   shortcut in the function;
% - Part C/F/G ground truths do not use piesos_quadvar; seeds
%   rng(20261002), rng(11);
% - Zop from opts.exclude has params{k}=[], which apply_sopvar cannot
%   handle: the ground truths expand it to zeros (fix_empty_params_zop);
% - Part H calls sossolve (needs an SDP solver).

nfail = 0;
ntotal = 0;

%% ---- Part A: input validation ----
[ok,msg] = run_case('A1: non-polyopvar Z errors', @case_A1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('A2: non-struct opts errors', @case_A2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('A3: invalid opts.psatz errors (piesos_poslpivar)', @case_A3);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('A4: monomial of states of different dependence errors', @case_A4);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

[ok,msg] = run_case('A5: different operators on the left and right error', @case_A5);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part B: structural sweep, prefix replica (nZ/Zop/Pmat/dedup) ----
cases = build_sweep_cases();
for i=1:numel(cases)
    cs = cases{i};
    [ok,msg] = run_case(sprintf('B%d: %s',i,cs.label), @() verify_prefix(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part C: real piesos_sosvar gives a NUMERICALLY CORRECT V ----
for i=1:numel(cases)
    cs = cases{i};
    [ok,msg] = run_case(sprintf('C%d: %s -- end-to-end numeric check',i,cs.label), @() verify_end_to_end(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part D: robustness against an inconsistent Z ----
[ok,msg] = run_case('D1: inconsistent (chkval=false) Z fails fast via initialize', @case_D1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part E: opts.psatz end-to-end ----
[ok,msg] = run_case('E1: invalid psatz errors, psatz=0 accepted', @case_E1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

%% ---- Part F: positivity, V against the sum of squares of the Cholesky factors ----
fcases =build_positivity_cases();
for i=1:numel(fcases)
    cs = fcases{i};
    [ok,msg] = run_case(sprintf('F%d: %s',i,cs.label), @() verify_positivity(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part G: independent checks for states of different dependence ----
% (verify_positivity, with states in their own variables, 2 spatial variables)
gcases = build_mixed_cases();
for i=1:numel(gcases)
    cs = gcases{i};
    [ok,msg] = run_case(sprintf('G%d: %s',i,cs.label), @() verify_positivity(cs));
    [nfail,ntotal] = tally(nfail,ntotal,ok,msg);
end

%% ---- Part H: SDPs with psatz multipliers ----
[ok,msg] = run_case('H1: psatz certifies int g*x^2 (g>=0 on [a,b] only), 1 variable', @case_H1);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);
[ok,msg] = run_case('H2: psatz in s2 certifies int g(s2)*y2^2, states of different dependence', @case_H2);
[nfail,ntotal] = tally(nfail,ntotal,ok,msg);

fprintf('\n================ Test_piesos_sosvar: %d / %d cases passed ================\n', ntotal-nfail, ntotal);
if nfail>0
    error('Test_piesos_sosvar: %d / %d cases FAILED.', nfail, ntotal);
else
    fprintf('Test_piesos_sosvar: ALL CASES PASSED.\n');
end

end


%% ========================================================================
function [ok,msg] = run_case(name,fn)
try
    fn();
    ok = true; msg = 'OK';
catch ME
    ok = false; msg = ME.message;
end
fprintf('%-68s %s\n', name, tern(ok,'PASS',['FAIL: ' msg]));
end

function [nfail,ntotal] = tally(nfail,ntotal,ok,~)
ntotal = ntotal+1;
if ~ok, nfail = nfail+1; end
end

function s = tern(c,a,b)
if c, s=a; else, s=b; end
end

function assert_err_contains(fn,needle,ctxmsg)
% Assert that calling fn() throws an error whose message contains 'needle'.
threw = false; msg = '';
try
    fn();
catch ME
    threw = true; msg = ME.message;
end
if ~threw
    error('%s: expected an error containing "%s", but no error was thrown.', ctxmsg, needle);
elseif ~contains(msg,needle)
    error('%s: expected error containing "%s", got "%s".', ctxmsg, needle, msg);
end
end


%% ========================================================================
%  Part A: input validation
%% ========================================================================
function case_A1()
polyopvar x1(s1);
prog = piesos_program(x1);
assert_err_contains(@() piesos_sosvar(prog,5), ...
    "'polyopvar'", 'A1');
end

function case_A2()
polyopvar x1(s1);
prog = piesos_program(x1);
assert_err_contains(@() piesos_sosvar(prog,x1,5), ...
    "'struct'", 'A2');
end

function case_A3()
% opts.psatz is passed on to piesos_poslpivar, which rejects invalid specifications
polyopvar x1(s1);
prog = piesos_program(x1);
assert_err_contains(@() piesos_sosvar(prog,x1,struct('psatz',[0 1])), ...
    '1 x N array', 'A3 (2 orders for 1 variable)');
assert_err_contains(@() piesos_sosvar(prog,x1,struct('psatz',-1)), ...
    'nonnegative integers', 'A3 (negative order)');
assert_err_contains(@() piesos_sosvar(prog,x1,struct('psatz',1.5)), ...
    'nonnegative integers', 'A3 (non-integer order)');
end

function case_A4()
% A monomial of states depending on different variables: monopmat gives all
% factors both variables, so V would have too many dummy variables (the guard
% errors before the loop). Same variables, or linear monomials, are accepted (Part G)
y1 = polyopvar('y1',polynomial({'s1'}),[0 1]);
y2 = polyopvar('y2',polynomial({'s1','s2'}),[0 1;0 1]);
Z = [y1;y2;otimes(y1,y2)];
prog = piesos_program([y1;y2]);
% sep: few rows, so without the guard the error would come only from combine_terms
assert_err_contains(@() piesos_sosvar(prog,Z,struct('sep',true)), 'different variables', 'A4');
end

function case_A5()
% piesos_quadvar allows different operators on the left and right (opts.L, opts.R,
% degs.L, degs.R), but V is then not SOS: error. Differing operators of equal row
% count: exclude lower vs upper integral, or dummy monomials 1,s vs 1,s^2.
% Equal operators are accepted, and match the common options.
polyopvar x1(s1);
prog = piesos_program(x1);
exL = struct('exclude',[false;true;false]);
exR = struct('exclude',[false;false;true]);
dL = struct('in',{{{[0;1]}}},'out',{{{[0;1]}}});
dR = struct('in',{{{[0;2]}}},'out',{{{[0;1]}}});
assert_err_contains(@() piesos_sosvar(prog,x1,struct('L',exL,'R',exR)), ...
    'same operators', 'A5 (opts.L~=opts.R)');
assert_err_contains(@() piesos_sosvar(prog,x1,struct('degs',struct('L',dL,'R',dR))), ...
    'same operators', 'A5 (degs.L~=degs.R)');
[progA,VA,PmatA,ZopA] = piesos_sosvar(prog,x1,struct('L',exL,'R',exL));
[progB,VB,PmatB,ZopB] = piesos_sosvar(prog,x1,exL);
assert(isequal(progA,progB) && isequal(ZopA,ZopB) && isequal(PmatA{1}.C,PmatB{1}.C) && ...
    isequal(VA.Z,VB.Z) && isequal(VA.C,VB.C), 'A5: opts.L=opts.R must match common options.');
[progA,VA,PmatA,ZopA] = piesos_sosvar(prog,x1,struct('degs',struct('L',dL,'R',dL)));
[progB,VB,PmatB,ZopB] = piesos_sosvar(prog,x1,struct('degs',dL));
assert(isequal(progA,progB) && isequal(ZopA,ZopB) && isequal(PmatA{1}.C,PmatB{1}.C) && ...
    isequal(VA.Z,VB.Z) && isequal(VA.C,VB.C), 'A5: degs.L=degs.R must match common degrees.');
end


%% ========================================================================
%  Part B-C: structural sweep
%% ========================================================================
function cases = build_sweep_cases()
cases = {};

c = struct(); c.label = 'single scalar state, degree 2, sep (orig. example)';
c.nstates = 1; c.Zexpr = 'x1'; c.degs = 2; c.opts_extra = struct('sep',true);
cases{end+1} = c;

c = struct(); c.label = 'two independent states, Z=[x1;x2]';
c.nstates = 2; c.Zexpr = '[x1;x2]'; c.degs = 1; c.opts_extra = struct();
cases{end+1} = c;

c = struct(); c.label = 'multi-monomial single state, Z=[x1;x1.*x1]';
c.nstates = 1; c.Zexpr = '[x1;otimes(x1,x1)]'; c.degs = 1; c.opts_extra = struct();
cases{end+1} = c;

c = struct(); c.label = 'cross term, Z=[x1;x2;x1.*x2]';
c.nstates = 2; c.Zexpr = '[x1;x2;otimes(x1,x2)]'; c.degs = 1; c.opts_extra = struct();
cases{end+1} = c;

c = struct(); c.label = 'opts.exclude sET (drop multiplier term)';
c.nstates = 1; c.Zexpr = 'x1'; c.degs = 2; c.opts_extra = struct('exclude',[true;false;false]);
cases{end+1} = c;

c = struct(); c.label = 'opts.sep true, single state';
c.nstates = 1; c.Zexpr = 'x1'; c.degs = 2; c.opts_extra = struct('sep',true);
cases{end+1} = c;

c = struct(); c.label = 'opts.sep false, single state';
c.nstates = 1; c.Zexpr = 'x1'; c.degs = 2; c.opts_extra = struct('sep',false);
cases{end+1} = c;

c = struct(); c.label = 'degree 3, single scalar state';
c.nstates = 1; c.Zexpr = 'x1'; c.degs = 3; c.opts_extra = struct('sep',true);
cases{end+1} = c;

end

function [prog,Z,opts] = build_case(cs)
% Build (prog,Z,opts) for a sweep case; x1 (x2) declared fresh in the function
% workspace, so no state leaks between cases.
polyopvar x1(s1); %#ok<NASGU>
if cs.nstates>=2
    polyopvar x2(s1); %#ok<NASGU>
end
Z = eval(cs.Zexpr); %#ok<NASGU>
if cs.nstates==1
    prog = piesos_program(x1);
else
    prog = piesos_program([x1;x2]);
end
opts = cs.opts_extra;
opts.degs = cs.degs;
end

function [prefix,Zop_expected] = independent_prefix(Z,opts)
% Ground truth: expected Zop from monopmat called directly (not via piesos_sosvar)
if ~Z.chkval
    Z = initialize(Z);
end
dims = Z.C.dims; vars = Z.C.vars; dom = Z.C.dom; depmat = Z.C.depmat;
% Every state mapped to the variables on which some state depends; primary
% degrees (degs.out) split off from the kernels of Zop (set to 0)
depmat.out = repmat(any(depmat.in,1),numel(dims.in),1);
if ~isfield(opts,'degs'), degs = 0; else, degs = opts.degs; end
Zop_expected = monopmat(dims.in,vars,dom,depmat,struct('in',degs,'out',0),opts);
prefix = struct('Z',Z,'dims',dims,'nZ_expected',numel(dims.in));
end

function verify_prefix(cs)
[prog,Z,opts] = build_case(cs); %#ok<ASGLU>

[gt,Zop_expected] = independent_prefix(Z,opts);
nZ_expected = gt.nZ_expected;
assert(nZ_expected==size(gt.Z.Z,1), ...
    'ground truth inconsistency: numel(dims.in)=%d ~= size(Z.Z,1)=%d (should always agree for a validated Z)', ...
    nZ_expected, size(gt.Z.Z,1));
assert(numel(Zop_expected.dims.out)==nZ_expected, ...
    'monopmat block count (%d) does not match expected nZ (%d) -- the (ridx,cidx) loop in piesos_sosvar indexes Zop.ops{ridx,ridx} up to nZ, so these must agree', ...
    numel(Zop_expected.dims.out), nZ_expected);

pfx = piesos_sosvar_prefix(prog,Z,opts);
assert(pfx.nZ==nZ_expected, 'piesos_sosvar prefix nZ=%d, expected %d', pfx.nZ, nZ_expected);
assert(isequal(size(pfx.Zop),size(Zop_expected)), ...
    'piesos_sosvar prefix Zop size %s ~= independently-built monopmat size %s', ...
    mat2str(size(pfx.Zop)), mat2str(size(Zop_expected)));
% Pmat: q x q cell, Pmat{i,j} a dpvar of size (rows of block i) x (rows of block j)
blk = Zop_expected.dims.out(:);
assert(iscell(pfx.Pmat) && isequal(size(pfx.Pmat),[nZ_expected,nZ_expected]), ...
    'Pmat should be a %d x %d cell', nZ_expected, nZ_expected);
for i=1:nZ_expected
    for j=1:nZ_expected
        assert(isa(pfx.Pmat{i,j},'dpvar') && isequal(size(pfx.Pmat{i,j}),[blk(i),blk(j)]), ...
            'Pmat{%d,%d} should be a dpvar of size %dx%d, got %s %s',i,j,blk(i),blk(j), ...
            class(pfx.Pmat{i,j}),mat2str(size(pfx.Pmat{i,j})));
    end
end

% Z_lin against an independently built expected set: all (i,j) in 1..nZ x 1..nZ,
% concatenating each pair's nonzero factors with zeros pushed to the end. Two
% pairs collapse to one Z_lin row only if their factor lists match exactly, so
% "x1 then x2" and "x2 then x1" are different rows (combine_terms merges them
% later). Plain double loop and unique(...,'rows'), independent of
% uniquerows_integerTable.
Zraw = full(gt.Z.Z);
d = size(Zraw,2);
Zfull_gt = zeros(nZ_expected^2,2*d);
row = 0;
for i=1:nZ_expected
    for j=1:nZ_expected
        row = row+1;
        combo = [Zraw(i,:), Zraw(j,:)];
        combo = [combo(combo>0), combo(combo==0)]; %#ok<AGROW>
        Zfull_gt(row,:) = combo;
    end
end
Zlin_gt = unique(Zfull_gt,'rows');
nZ_lin_expected = size(Zlin_gt,1);
assert(size(pfx.Z_lin,1)==nZ_lin_expected, ...
    'Z_lin has %d distinct monomials, expected %d from independent (i,j) enumeration over all %d^2 pairs', ...
    size(pfx.Z_lin,1), nZ_lin_expected, nZ_expected);
% Every independently built row appears exactly once among the rows of Z_lin
Zlin_actual = full(pfx.Z_lin);
if size(Zlin_actual,2)<size(Zlin_gt,2)
    Zlin_actual = [Zlin_actual, zeros(size(Zlin_actual,1),size(Zlin_gt,2)-size(Zlin_actual,2))];
elseif size(Zlin_gt,2)<size(Zlin_actual,2)
    Zlin_gt = [Zlin_gt, zeros(size(Zlin_gt,1),size(Zlin_actual,2)-size(Zlin_gt,2))]; %#ok<AGROW>
end
assert(isequal(sortrows(Zlin_actual),sortrows(Zlin_gt)), ...
    'Z_lin content does not bijectively match the independently-enumerated monomial set');
% old2new_idcs: one entry per (i,j) pair (nZ^2), each a real row of Z_lin
assert(numel(pfx.old2new_idcs)==nZ_expected^2, ...
    'old2new_idcs should have one entry per (i,j) pair considered (nZ^2=%d), got %d', ...
    nZ_expected^2, numel(pfx.old2new_idcs));
assert(all(pfx.old2new_idcs(pfx.blk_nums)>=1 & pfx.old2new_idcs(pfx.blk_nums)<=size(pfx.Z_lin,1)), ...
    'every (ridx,cidx) pair should map to a valid row of Z_lin');
assert(numel(unique(pfx.old2new_idcs(pfx.blk_nums)))==size(pfx.Z_lin,1), ...
    'every row of Z_lin should be reached by at least one (ridx,cidx) pair (the tril+doubling defect this replaced left some rows unreached)');

end


%% ========================================================================
%  Part C: real piesos_sosvar gives a NUMERICALLY CORRECT V (end-to-end)
%% ========================================================================
function verify_end_to_end(cs)
% Build concrete prog/Z/opts, call the real piesos_sosvar, and check the result
% two independent ways.
%
% (1) PAIRWISE (pairwise_formula_total): for every (ridx,cidx) in 1:nZ^2 recompute
%     the block's contribution from Zop and Pmat only (the outputs of this call,
%     not piesos_sosvar's internals), evaluated by apply_op at a concrete state
%     and decision variable assignment, and summed. This equals
%     <Zop*Z,Pmat*Zop*Z>, checked against 'ground_truth_via_Zop' (apply_op(Zop,..)
%     and plain polynomial integration, no int()/intop).
%
% (2) FULL V (eval_V_at_state): the returned V, evaluated by walking V.Z row by
%     row and applying the matching V.C.ops{1,l} 'intop' via apply_op. Relies on
%     V.C.ops{1,l} having its dummy order equal to the factor order of V.Z(l,:),
%     for every row also after merging pairs with different factor counts (x1.*x2
%     against x1 or x2).
[prog,Z,opts] = build_case(cs); %#ok<ASGLU>
[prog,V,Pmat,Zop] = piesos_sosvar(prog,Z,opts); %#ok<ASGLU>

[s1var,xvals,d0map] = concrete_values(Z,Pmat);
Zraw = Z_by_rows(Z);
dom_ab = Z.C.dom.in(1,:);

Vtrue = ground_truth_via_Zop(Zop,Pmat,xvals,Zraw,s1var,dom_ab,d0map);
tol = 1e-6*max(1,abs(Vtrue));

% ---- Check (1): pairwise formula ----
Vpairwise = pairwise_formula_total(Zop,Pmat,xvals,Zraw,s1var,d0map);
assert(abs(Vpairwise-Vtrue) < tol, ...
    'pairwise formula total %.10f does not match ground truth %.10f (kron/otimes fix or nZ^2 accumulation is broken)', ...
    Vpairwise, Vtrue);

% ---- Check (2): full V ----
[Vcand,ok2,errmsg2] = try_eval_V_at_state(V,Z,xvals,d0map,s1var);
if ~ok2
    error('full-V evaluation errored unexpectedly: %s', errmsg2);
end
assert(abs(Vcand-Vtrue) < tol, ...
    'full-V evaluation %.10f does not match ground truth %.10f', Vcand, Vtrue);
end

function Zr = Z_by_rows(Z)
% Z.Z with its rows in the order of the entries of Z (rows of Z.C), the order of
% the blocks of Zop and Pmat: the monomial of entry i is the row of Z.Z with a
% block in row i of Z.C (Z.Z can order the monomials differently, e.g. for
% Z=[x1;x2;x1.*x2])
[~,kmon] = max(~cellfun(@isempty,Z.C.ops),[],2);
Zr = full(Z.Z);
Zr = Zr(kmon,:);
end

function [s1var,xvals,d0map] = concrete_values(Z,Pmat)
% Reproducible concrete values: one distinct low-degree polynomial per state, and
% one decision variable value per Pmat.dvarname (containers.Map by name, as
% eval_dpvar_at expects).
if ~Z.chkval, Z = initialize(Z); end
vn = Z.C.vars.in;
s1var = polynomial(vn(1)); s1var = s1var(1);
nstates = numel(Z.varname);
xvals = cell(1,nstates);
for i=1:nstates
    xvals{i} = (2+i)*s1var.^2 - i*s1var + (1+0.5*i);
end
rng(20261002);
% Pmat is a cell of dpvar blocks: all decision variables appearing in any block
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)];                                    %#ok<AGROW>
end
names = unique(names,'stable');
d0map = containers.Map(names, num2cell(0.5+rand(numel(names),1)));
end

function Vtot = ground_truth_via_Zop(Zop,Pmat,xvals,Zraw,s1var,dom_ab,d0map)
% <Zop*Z,Pmat*Zop*Z>, via apply_op(Zop,...) directly, not through int()/intop
Zop = fix_empty_params_zop(Zop);
Pmat_num = assemble_numeric(Pmat,d0map);
nZst = size(Zraw,1);
Xcell_zop = cell(1,nZst);
for k=1:nZst
    idcs0 = Zraw(k,Zraw(k,:)~=0);
    if isempty(idcs0), Xcell_zop{k} = {1}; else, Xcell_zop{k} = xvals(idcs0); end
end
Yzop = apply_op(Zop,Xcell_zop);
Vtot = double(int(Yzop.'*Pmat_num*Yzop, s1var, dom_ab(1), dom_ab(2)));
end

function Pnum = assemble_numeric(Pmat,d0map)
% The (sum(mi) x sum(mj)) 'polynomial' matrix of the cell of dpvar blocks Pmat,
% at the decision variable values d0map (by name)
[q1,q2] = size(Pmat);
rows = cell(q1,1);
for i=1:q1
    blks = cell(1,q2);
    for j=1:q2
        blks{j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
    rows{i} = horzcat(blks{:});
end
Pnum = vertcat(rows{:});
end

function Zop2 = fix_empty_params_zop(Zop)
% monopmat leaves params{k}=[] for terms excluded via opts.exclude, which
% apply_sopvar cannot handle ("Matrix dimensions must agree" in its mtimes, for
% opts.exclude=[true;false;false]): expand each [] to a zero sparse matrix of the
% right size. Needed for ground_truth_via_Zop only (apply_op on intops avoids it).
Zop2 = Zop;
[p,d] = size(Zop.ops);
for i=1:p
    for j=1:d
        if ~isempty(Zop.ops{i,j})
            Zop2.ops{i,j} = fix_empty_params_tensopvar(Zop.ops{i,j});
        end
    end
end
end

function T2 = fix_empty_params_tensopvar(T)
T2 = T;
[nt,nf] = size(T.ops);
for i=1:nt
    for j=1:nf
        entry = T.ops{i,j};
        if isa(entry,'sopvar')
            T2.ops{i,j} = fill_empty_params_sopvar(entry);
        elseif isa(entry,'tensopvar')
            T2.ops{i,j} = fix_empty_params_tensopvar(entry);
        end
    end
end
end

function Pop2 = fill_empty_params_sopvar(Pop)
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

function Vtot = pairwise_formula_total(Zop,Pmat,xvals,Zraw,s1var,d0map)
% Sum over every (ridx,cidx) in 1:nZ^2 of the block contribution, built from
% Zop/Pmat with the construction of the piesos_quadvar loop: Kij =
% int(otimes(repmat(CL,n,1),repelem(CR,m,1),[1,2]),Pij(1:end)), whose native
% dummy order [CL factors, CR factors] is reordered (reorder_vars) to Z_lin's
% [cidx factors, ridx factors] before substituting the states. The domain is
% already in int(). Deliberately a full, independent 1:nZ^2 loop (no sharing of
% int() calls between pairs), as a check of any symmetry shortcut in the function.
nZst = size(Zraw,1);
Vtot = 0;
for ridx=1:nZst
    for cidx=1:nZst
        CL = Zop.ops{ridx,ridx};   CR = Zop.ops{cidx,cidx};
        m = size(CL,1);            n = size(CR,1);
        Pij = Pmat{ridx,cidx};
        Kij = int(otimes(repmat(CL,n,1),repelem(CR,m,1),[1,2]),Pij(1:end));
        A = CL.deps.in;   B = CR.deps.in;
        var_order = cell(1,numel(A));
        for jj=1:numel(A)
            var_order{jj} = [A(jj)+1:A(jj)+B(jj), 1:A(jj)];
        end
        Kij = reorder_vars(Kij,var_order);
        idcs_r = Zraw(ridx,Zraw(ridx,:)~=0);
        idcs_c = Zraw(cidx,Zraw(cidx,:)~=0);
        dummy_src_idcs = [idcs_c, idcs_r];
        vterm = build_vterm(Kij,dummy_src_idcs,xvals,s1var);
        Knum = eval_dpvar_at(apply_op(Kij,vterm),d0map);
        Vtot = Vtot + double(Knum);
    end
end
end

function vterm = build_vterm(Kij,dummy_src_idcs,xvals,s1var)
% Dummy-substituted product that apply_op expects for a single Kij: one factor
% xvals{idx} at the corresponding dummy name, per entry of dummy_src_idcs, in order.
dumnames = Kij.varnames.in_dum{1};
if isempty(dummy_src_idcs)
    vterm = polynomial(1);
    return
end
if numel(dumnames)~=numel(dummy_src_idcs)
    error('build_vterm: dummy count %d ~= factor count %d', numel(dumnames), numel(dummy_src_idcs));
end
vterm = polynomial(1);
for f=1:numel(dummy_src_idcs)
    xf = xvals{dummy_src_idcs(f)};
    dpf = polynomial(dumnames(f)); dpf = dpf(1);
    vterm = vterm*subs(xf,s1var,dpf);
end
end

function [Vnum,ok,errmsg] = try_eval_V_at_state(V,Z,xvals,d0map,s1var)
ok = true; errmsg = '';
try
    Vnum = eval_V_at_state(V,Z,xvals,d0map,s1var);
catch ME
    ok = false; errmsg = ME.message; Vnum = NaN;
end
end

function Ynum = eval_V_at_state(V,Z,xvals,d0map,s1var) %#ok<INUSD>
% Evaluator for the returned V: walks V.Z row by row and applies the matching
% V.C.ops{1,l} 'intop' via apply_op. V.C.ops{1,l}'s dummy order is assumed to
% match the factor order of V.Z(l,:) position for position (Z is unused).
if ~V.chkval, V = initialize(V); end
Vmat = full(V.Z);
nZ_lin_V = size(Vmat,1);
total = polynomial(0);
for l=1:nZ_lin_V
    idcs_full_V = Vmat(l, Vmat(l,:)~=0);
    Kent = V.C.ops{1,l};
    if isempty(Kent)
        continue
    end
    if isnumeric(Kent)
        if all(Kent(:)==0), continue; end
        error('eval_V_at_state: unexpected nonzero numeric V.C entry at row %d',l);
    end

    vars_in = Kent.vars.in;
    if numel(vars_in)~=1
        error('eval_V_at_state: expected 1 shared spatial var, got %d',numel(vars_in));
    end
    vterm = build_vterm(Kent,idcs_full_V,xvals,s1var);
    Knum = eval_dpvar_at(apply_op(Kent,vterm),d0map);
    total = total + Knum;
end
Ynum = double(total);
end


%% ========================================================================
%  Part D: robustness against an inconsistent Z
%% ========================================================================
function case_D1()
polyopvar x1(s1);
Z = combine_terms([x1; otimes(x1,x1)]);    % 2 distinct monomials: x1, x1^2
assert(Z.chkval, 'test setup: Z should be chkval=true after combine_terms');
assert(size(Z.Z,1)==2 && size(Z.C.ops,2)==2, 'test setup: expected a 2-monomial, 2-column Z');

% Corrupt: truncate Z.Z to 1 row, leave Z.C at 2 columns, mark unvalidated. Without
% the 'initialize' call piesos_sosvar would silently under-count the monomial
% basis; with it, initialize's consistency check fires before Zop/Pmat are built.
Zbad = Z;
Zbad.Z = Z.Z(1,:);
Zbad.chkval = false;

prog = piesos_program(x1);
opts = struct(); opts.degs = 1;
assert_err_contains(@() piesos_sosvar(prog,Zbad,opts), ...
    'number of monomials must match', 'D1');
end


%% ========================================================================
%  Part E: opts.psatz, end-to-end
%% ========================================================================
function case_E1()
% psatz=0 is the plain SOS case with one declared SOS variable, psatz=1 declares
% one for each multiplier (2); values of V: Part F/G
polyopvar x1(s1);
prog = piesos_program(x1);
assert_err_contains(@() piesos_sosvar(prog,x1,struct('degs',1,'psatz',[0,1])), 'psatz', 'E1');
n0 = prog.var.num;
prog0 = piesos_sosvar(prog,x1,struct('degs',1,'psatz',0));
assert(prog0.var.num==n0+1,'psatz=0 must declare a single SOS variable.');
prog1 = piesos_sosvar(prog,x1,struct('degs',1,'psatz',1));
assert(prog1.var.num==n0+2,'psatz=1 must declare an SOS variable for each of the 2 multipliers.');
end


%% ========================================================================
%  Parts F-H: positivity of V, states of different dependence, SDPs
%% ========================================================================
function c = mkcase(label,vn,dom,states,Zexpr,din,dout,opts_extra,psatz)
% Case for verify_positivity: spatial variables vn (domains dom); states{i} =
% indices (in vn) of the variables of the state x<i>; Zexpr the basis in terms
% of x1,x2,..; degrees din (dummy variables, kernels of Zop) and dout (primary
% variables, Gram basis of Pmat), the same in all variables; opts_extra
% (exclude, sep); psatz the 1 x N orders of the psatz multipliers, or []
c = struct('label',label,'vn',{vn},'dom',dom,'states',{states},'Zexpr',Zexpr,...
    'din',din,'dout',dout,'opts_extra',opts_extra,'psatz',psatz);
end

function cases = build_positivity_cases()
% One spatial variable, s1 in [0.5,2]
cases = {};
vn = {'s1'};    dom = [0.5 2];
cases{end+1} = mkcase('single state, in 2, out 2, sep',vn,dom,{1},'x1',2,2,struct('sep',true),[]);
cases{end+1} = mkcase('two states, Z=[x1;x2]',vn,dom,{1,1},'[x1;x2]',1,1,struct(),[]);
cases{end+1} = mkcase('Z=[x1;x1.*x1]',vn,dom,{1},'[x1;otimes(x1,x1)]',1,1,struct(),[]);
cases{end+1} = mkcase('cross term Z=[x1;x2;x1.*x2]',vn,dom,{1,1},'[x1;x2;otimes(x1,x2)]',1,1,struct(),[]);
cases{end+1} = mkcase('constant Pmat (out 0), in 2',vn,dom,{1},'x1',2,0,struct(),[]);
cases{end+1} = mkcase('kernels of degree 0, Pmat of degree 4 (in 0, out 2)',vn,dom,{1},'x1',0,2,struct(),[]);
cases{end+1} = mkcase('exclude the multiplier term',vn,dom,{1},'x1',2,1,struct('exclude',[true;false;false]),[]);
cases{end+1} = mkcase('psatz 1',vn,dom,{1},'x1',1,2,struct(),1);
cases{end+1} = mkcase('psatz 2',vn,dom,{1},'x1',0,2,struct(),2);
cases{end+1} = mkcase('two states, psatz 1',vn,dom,{1,1},'[x1;x2]',1,2,struct(),1);
cases{end+1} = mkcase('Z=[x1;x1.*x1], psatz 1',vn,dom,{1},'[x1;otimes(x1,x1)]',0,1,struct(),1);
end

function cases = build_mixed_cases()
% Two spatial variables; states depending on different variables
cases = {};
vn = {'s1','s2'};    dom = [0.5 2;-1 0.7];
cases{end+1} = mkcase('y1(s1), y2(s1,s2)',vn,dom,{1,[1 2]},'[x1;x2]',1,1,struct(),[]);
cases{end+1} = mkcase('y1(s1), y2(s1,s2), psatz [0 1] in s2',vn,dom,{1,[1 2]},'[x1;x2]',1,2,struct(),[0 1]);
cases{end+1} = mkcase('y1(s1), y2(s1,s2), psatz [1 1] (4 terms)',vn,dom,{1,[1 2]},'[x1;x2]',0,1,struct(),[1 1]);
cases{end+1} = mkcase('y1(s1), y2(s1,s2), y3(s2)',vn,dom,{1,[1 2],2},'[x1;x2;x3]',1,1,struct(),[]);
cases{end+1} = mkcase('y1(s1) and y2(s1,s2), exclude multiplier',vn,dom,{1,[1 2]},'[x1;x2]',1,1,struct('exclude',[true;false;false]),[]);
cases{end+1} = mkcase('y1(s1) and y2(s1,s2), sep',vn,dom,{1,[1 2]},'[x1;x2]',1,1,struct('sep',true),[]);
cases{end+1} = mkcase('y1(s1), y2(s1,s2), Z=[x1;x2;x1.*x1]',vn,dom,{1,[1 2]},'[x1;x2;otimes(x1,x1)]',0,0,struct('sep',true),[]);
cases{end+1} = mkcase('disjoint y1(s1), y2(s2), psatz [0 1]',vn,dom,{1,2},'[x1;x2]',1,2,struct(),[0 1]);
% The first state not depending on s1 (the left factor CL of the pair has other
% variables than the right one, listed after them in Kij)
cases{end+1} = mkcase('y1(s2), y2(s1)',vn,dom,{2,1},'[x1;x2]',1,1,struct(),[]);
cases{end+1} = mkcase('y1(s2), y2(s1,s2)',vn,dom,{2,[1 2]},'[x1;x2]',1,1,struct(),[]);
% Three spatial variables, kernels of degree 0 (mixed states in 3 variables
% rely on misc/otimes_params and @intop/clean_op as fixed 10/09/2026)
vn3 = {'s1','s2','s3'};    dom3 = [0.5 2;-1 0.7;1 3];
cases{end+1} = mkcase('3 variables: y1(s1,s2), y2(s3)',vn3,dom3,{[1 2],3},'[x1;x2]',0,0,struct(),[]);
cases{end+1} = mkcase('3 variables: y1(s1,s2), y2(s1,s2,s3)',vn3,dom3,{[1 2],[1 2 3]},'[x1;x2]',0,0,struct(),[]);
cases{end+1} = mkcase('3 variables: y1(s1,s2), y2(s3), psatz [0 0 1]',vn3,dom3,{[1 2],3},'[x1;x2]',0,1,struct(),[0 0 1]);
end

function [Z,prog0,opts,dep_state] = build_generic(cs)
% Basis Z in the states x<i>, program, options, and the dependence of each
% state on the variables of Z.C.vars.in (rows in the order of Z.varname)
N = numel(cs.vn);
for i=1:numel(cs.states)
    idx = cs.states{i};
    S = polyopvar(sprintf('x%d',i),polynomial(cs.vn(idx)'),cs.dom(idx,:)); %#ok<NASGU>
    eval(sprintf('x%d = S;',i));
end
Z = eval(cs.Zexpr);
if ~Z.chkval, Z = initialize(Z); end
allst = eval(['[',strjoin(arrayfun(@(i) sprintf('x%d',i),1:numel(cs.states),'UniformOutput',false),';'),']']);
prog0 = piesos_program(allst);
opts = cs.opts_extra;
opts.degs = struct('in',cs.din,'out',cs.dout);
if ~isempty(cs.psatz)
    opts.psatz = cs.psatz;
end
assert(isequal(reshape(Z.C.vars.in,1,[]),reshape(cs.vn,1,[])),'Variables of Z are not ordered as in the case.');
dep_state = false(numel(Z.varname),N);
for i=1:numel(cs.states)
    dep_state(strcmp(Z.varname,sprintf('x%d',i)),cs.states{i}) = true;
end
end

function Z = zvec_poly(E,varnames)
% Monomial vector with exponents E{k} (column vectors) in the variable
% varnames{k}, variable 1 slowest, from the mixed-radix expansion of the
% index (independent of the code under test)
if isempty(E)
    Z = polynomial(1);
    return
end
nE = cellfun(@numel,E);
nZ = prod(nE);
pv = polynomial(varnames(:));
Zc = cell(nZ,1);
for idx=0:nZ-1
    r = idx;
    mon = polynomial(1);
    for k=numel(E):-1:1
        ik = mod(r,nE(k))+1;
        r = floor(r/nE(k));
        mon = mon*pv(k)^E{k}(ik);
    end
    Zc{idx+1} = mon;
end
Z = vertcat(Zc{:});
end

function g = psatz_terms(dom,ps,lin)
% Multipliers prod_k (s_k-a_k)^jk*(b_k-s_k)^(qk-jk) of the terms with linear
% index lin (first variable fastest) in the array of size ps+1, as functions
% of the vector s of the values of the variables
g = cell(1,numel(lin));
for t=1:numel(lin)
    jj = cell(1,numel(ps));
    [jj{:}] = ind2sub([ps+1,ones(1,max(0,2-numel(ps)))],lin(t));   % size vector of >=2 elements
    j = cell2mat(jj)-1;
    g{t} = @(s) psatz_prod(s,dom,ps,j);
end
end

function val = psatz_prod(s,dom,ps,j)
% prod_k (s_k-a_k)^jk*(b_k-s_k)^(qk-jk), for s numeric or a polynomial vector
val = 1;
for k=1:numel(ps)
    val = val*(s(k)-dom(k,1))^j(k)*(dom(k,2)-s(k))^(ps(k)-j(k));
end
end

function verify_positivity(cs)
% V=<Zop*Z,Pmat*Zop*Z> of the real piesos_sosvar, for the Gram matrices set to
% random positive definite matrices Qk (names from direct sosquadvar calls,
% for monomials Zs written out by hand), against
%   sum_k int g_k |Rk*(I x Zs)*Y|^2,   Qk = Rk'*Rk,   Y = Zop*(states),
% with Y by polynomial arithmetic (apply_op, independent of int/intop)
% and the cholesky factors Rk (without Pmat), and against the integral of the
% quadratic form with Pmat at the same values. In particular V>0. States are
% polynomials in their own variables, the integral over the variables of Z.
[Z,prog0,opts,dep_state] = build_generic(cs);
[prog,V,Pmat,Zop] = piesos_sosvar(prog0,Z,opts);
N = numel(cs.vn);
vnZ = Z.C.vars.in;      domZ = Z.C.dom.in;      dep = Z.C.depmat.in;
used = find(any(dep>0,1));          nu = numel(used);
vu = vnZ(used);         domu = domZ(used,:);
m = Zop.dims.out(:);    q = numel(m);
assert(iscell(Pmat) && isequal(size(Pmat),[q,q]),'Pmat must be a q x q cell.');
% psatz terms and monomials of the basis
ps = zeros(1,N);
if ~isempty(cs.psatz), ps = cs.psatz; end
ps_u = ps(used);
if isempty(used), error('No spatial variables.'); end
lin = 1:prod(ps_u+1);
g = psatz_terms(domu,ps_u,lin);
T = numel(g);
E = cell(1,q);
for i=1:q
    E{i} = cell(1,nu);
    for k=1:nu
        e = (0:cs.dout)';
        E{i}{k} = e(e<=max(e)-ceil(ps_u(k)/2));
    end
end
Zs = cell(q,1);
for i=1:q
    Zs{i} = zvec_poly(E{i},vu);
end
% Gram variables of direct sosquadvar calls, in the order of declaration
pr = prog0;     Qn = cell(1,T);
for t=1:T
    [pr,~,Qn{t}] = sosquadvar(pr,Zs,Zs,m,m,'pos');
end
assert(prog.var.num==prog0.var.num+T && isequal(prog.var.idx,pr.var.idx), ...
    'piesos_sosvar must declare %d SOS variable(s) as sosquadvar does.',T);
nz = cellfun(@(e) prod(cellfun(@numel,e)),E);
off = [0;cumsum(m(:).*nz(:))];
rng(11);
names = cell(0,1);      vals = zeros(0,1);
R = cell(1,T);
for t=1:T
    Rt = randn(off(end));
    Qt = Rt'*Rt/off(end) + 0.1*eye(off(end));       Qt = (Qt+Qt')/2;
    R{t} = chol(Qt);
    for i=1:q
        for j=1:q
            nm = Qn{t}{i,j};
            vl = Qt(off(i)+1:off(i+1),off(j)+1:off(j+1));
            names = [names;nm(:)];          %#ok<AGROW>
            vals = [vals;vl(:)];            %#ok<AGROW>
        end
    end
end
[names,iu] = unique(names,'stable');
d0map = containers.Map(names,num2cell(vals(iu)));
% Random states in their own variables
sv = polynomial(vnZ(:));
X = cell(1,numel(Z.varname));
for i=1:numel(Z.varname)
    own = find(dep_state(i,:));
    p = polynomial(1+0.5*i);
    for k=own
        p = p + (0.3*i)*sv(k) - (0.2*i)*sv(k)^2;
    end
    if numel(own)>1, p = p + 0.15*sv(own(1))*sv(own(2)); end
    X{i} = p;
end
% Zop*(monomials) by polynomial arithmetic
Zf = Z_by_rows(Z);
Xc = cell(1,size(Zf,1));
for k=1:size(Zf,1)
    idc = Zf(k,Zf(k,:)~=0);
    if isempty(idc), Xc{k} = {1}; else, Xc{k} = X(idc); end
end
Y = apply_op(fix_empty_params_zop(Zop),Xc);
% (a) sum of squares of the Cholesky factors
offm = [0;cumsum(m(:))];
blocks = cell(q,1);
for i=1:q
    yi = Y(offm(i)+1:offm(i+1));
    parts = cell(m(i),1);
    for r=1:m(i)
        parts{r} = yi(r)*Zs{i};
    end
    blocks{i} = vertcat(parts{:});
end
Wv = vertcat(blocks{:});
Vsos = 0;
for t=1:T
    Wt = R{t}*Wv;
    Ft = g{t}(sv(used))*(Wt.'*Wt);
    for k=used
        Ft = int(Ft,sv(k),domZ(k,1),domZ(k,2));
    end
    Vsos = Vsos + double(Ft);
end
% (b) the quadratic form with Pmat at the same values
Pnum = assemble_numeric(Pmat,d0map);
F = Y.'*Pnum*Y;
for k=used
    F = int(F,sv(k),domZ(k,1),domZ(k,2));
end
Vtrue = double(F);
% (c) the returned functional V
vals_x = struct();
for i=1:numel(Z.varname)
    vals_x.(Z.varname{i}) = X{i};
end
Vd = eval_V_generic(V,vals_x,dep_state,vnZ);
Vlib = double(eval_dpvar_at(Vd,d0map));
tol = 1e-8*(1+abs(Vsos));
assert(Vsos>0,'The sum of squares must be positive, got %g.',Vsos);
assert(abs(Vtrue-Vsos)<=tol,'<Zop*x,Pmat*Zop*x> = %.10g differs from the sum of squares %.10g.',Vtrue,Vsos);
assert(abs(Vlib-Vsos)<=tol,'V(x) = %.10g differs from the sum of squares %.10g.',Vlib,Vsos);
end

function Y = eval_V_generic(V,vals,dep_state,vnames)
% Evaluate the functional V (polyopvar of intop kernels) at the states
% vals.(name), polynomials in their own variables vnames. The dummy variables
% of s_k are numbered over the factors of the monomial that depend on s_k, and
% named s_k_j (as in apply_op) if the monomial has more than 1 factor.
% dep_state(i,k): state i (in the order of V.varname) depends on s_k.
Y = 0;
nZ = size(V.Z,1);
Zf = full(V.Z);
for k=1:nZ
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        Y = Y + apply_op(op,{1});
    elseif isa(op,'intop')
        slot = zeros(1,numel(vnames));
        v = 1;
        for j=1:d
            xj = vals.(V.varname{idcs(j)});
            for mm=1:numel(vnames)
                if dep_state(idcs(j),mm)
                    slot(mm) = slot(mm)+1;
                    if d>1
                        isv = strcmp(xj.varname,vnames{mm});
                        xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                    end
                end
            end
            v = kron(v,xj);
        end
        Y = Y + apply_op(op,v);
    else
        xv = cell(1,d);
        for j=1:d, xv{j} = vals.(V.varname{idcs(j)}); end
        Y = Y + apply_op(op,xv);
    end
end
end

function feas = h_feasible(Z,opts,Vt)
% Feasibility of V==Vt for V the distributed SOS functional of piesos_sosvar
prog = piesos_program(Z);
[prog,V] = piesos_sosvar(prog,Z,opts);
prog = piesos_eq(prog,V-Vt);
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
feas = (info.pinf==0 && info.dinf==0);
end

function case_H1()
% V = int g(s)*x(s)^2 ds, g=(s-a)(b-s), is nonnegative on [a,b] only: it is
% an SOS functional with a psatz multiplier (g_1 of order 2), not without,
% and -V is not for any multiplier. Control: int x^2 needs no multiplier
a = 0.5;    b = 2;
x1 = polyopvar('x1',polynomial({'s1'}),[a b]);
s = polynomial({'s1'});
g = (s-a)*(b-s);
Vg = innerprod(x1,g,[]);
Vc = innerprod(x1,1,[]);
o0 = struct('degs',struct('in',0,'out',2));
o2 = o0;    o2.psatz = 2;
assert(h_feasible(x1,o0,Vc),'int x^2 must be feasible without psatz.');
assert(h_feasible(x1,o2,Vg),'int g*x^2 must be feasible with psatz 2.');
assert(~h_feasible(x1,o0,Vg),'int g*x^2 is not an SOS functional without psatz.');
assert(~h_feasible(x1,o2,-Vg),'-int g*x^2 is not nonnegative: infeasible with psatz 2.');
end

function case_H2()
% States of different dependence, y1(s1) and y2(s1,s2): int g(s2)*y2^2, with
% g=(s2-a2)(b2-s2), is certified with a psatz multiplier of order 2 in s2
% (in the common space (s1,s2) the multiplier acts on all blocks, also that of
% y1 which does not depend on s2), not without; -int g*y2^2 never
dom = [0.5 2;-1 0.7];
y1 = polyopvar('y1',polynomial({'s1'}),dom(1,:));
y2 = polyopvar('y2',polynomial({'s1','s2'}),dom);
Z = [y1;y2];
s = polynomial({'s1','s2'});
g = (s(2)-dom(2,1))*(dom(2,2)-s(2));
Vg = innerprod(y2,g,[]);
o0 = struct('degs',struct('in',0,'out',[0 2]));
o2 = o0;    o2.psatz = [0 2];
assert(h_feasible(Z,o2,Vg),'int g(s2)*y2^2 must be feasible with psatz [0 2].');
assert(~h_feasible(Z,o0,Vg),'int g(s2)*y2^2 is not an SOS functional without psatz.');
assert(~h_feasible(Z,o2,-Vg),'-int g(s2)*y2^2 is not nonnegative: infeasible with psatz [0 2].');
end


%% ========================================================================
%  Prefix replica of piesos_sosvar.m (input validation through the nZ/Z_full/
%  uniquerows monomial deduplication), up to but not including the (ridx,cidx)
%  loop. Kept so Zop/Pmat/nZ/Z_lin can be inspected independent of V.
%% ========================================================================
function pfx = piesos_sosvar_prefix(prog,Z,opts) %#ok<INUSD>
if ~isa(Z,'polyopvar')
    error("Basis of distributed monomials must be specified as 'polyopvar' object.")
end
if ~Z.chkval
    Z = initialize(Z);
end
if nargin<=2
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
end
if ~isfield(opts,'degs')
    degs = 0;
else
    degs = opts.degs;
end

dims = Z.C.dims; vars = Z.C.vars; dom = Z.C.dom; depmat = Z.C.depmat;

[~,Pmat,Zop] = piesos_poslpivar(prog,dims.in,vars.in,dom.in,depmat.in,degs,opts);

nZ = numel(dims.in);
Z_full = [repelem(Z.Z,nZ,1), repmat(Z.Z,nZ,1)];
for i=1:size(Z_full,1)
    Z_i = Z_full(i,:);
    Z_full(i,:) = [Z_i(Z_i>0), Z_i(Z_i==0)];
end
[M,Z_lin] = uniquerows_integerTable([sum(Z_full>0,2),Z_full]);
Z_lin = Z_lin(:,2:end);
old2new_idcs = M*(1:size(Z_lin,1))';
blk_nums = 1:nZ^2;

pfx = struct('nZ',nZ,'Zop',Zop,'Pmat',{Pmat},'Z_lin',Z_lin, ...
    'old2new_idcs',old2new_idcs,'blk_nums',blk_nums);
end
