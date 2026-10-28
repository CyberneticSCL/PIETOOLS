function T = bench_getsol_sop(opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T = BENCH_GETSOL_SOP(OPTS) measures the wall time and memory of
% getsol_lpivar_sop(prog,P) over a sweep of the number q of decision
% variables of P and of the number nv of spatial variables, and compares
% it with the legacy routines in one variable. The values of the decision
% variables are random (prog.solinfo.RRx, drawn per name); no LPI is
% solved. With s1, ..., s_nv on [0,1] and m the component count that gives
% about q decision variables, P is
%   'lpivar'  for each nv in OPTS.nv and q in OPTS.q, a square decision
%             operator on L2^m[s1,...,s_nv] from lpivar_cdopvar at degree
%             OPTS.deg;
%   'multi'   when OPTS.multi, for each nv and q = OPTS.q(end), a square
%             decision operator on R^m x L2^m[s1,...,s_nv] from
%             lpivar_cdopvar;
%   'gram'    for each nv in OPTS.gram_nv and q in OPTS.gram_q, a positive
%             decision operator on L2^m[s1,...,s_nv] from poscopvar at
%             degree OPTS.deg, whose names appear twice in decvartable;
% each call timed as the least of OPTS.reps wall times, with the
% decvartable of the program ('own') and with that table extended by 2q
% names P does not use, shuffled ('pad'); the two results must be equal
% (isequal), and the MATLAB profiler records the PeakMem and
% TotalMemAllocated of one more 'pad' call. Then
%   'legacy'  for each q in OPTS.legacy_q, on L2^m[s1] at degrees
%             OPTS.deg*[1 1 1], each timed once: getsol_lpivar_sop on the
%             'cdopvar' from lpivar_cdopvar, getsol_lpivar on the 'dopvar'
%             of the same family from lpivar, and lpigetsol on the 'dpvar'
%             of the names of the 'cdopvar', which must return the last q
%             entries of RRx. Once the time projected for the next q
%             exceeds OPTS.legacy_cap (linear in q for lpivar and
%             getsol_lpivar, as q^2.25 for lpigetsol), that route is
%             skipped.
%
% INPUT
% - opts:   (optional) 'struct' with optional fields (defaults): q ([1e4
%           1e5 1e6 3e6]), nv ([1 2 3]), deg (1), reps (3), gram_q ([1e4
%           1e5 1e6]), gram_nv (1), legacy_q ([1e3 1e4 1e5 1e6 3e6]),
%           legacy_cap (120, seconds), multi (true); the entries of nv and
%           gram_nv in 1..3. An empty q, gram_q or legacy_q skips its
%           part; 'multi' needs a nonempty q.
%
% OUTPUT
% - T:      'struct' with 'table' fields lpivar, multi and gram, columns
%           family, nv, q (numel(P.Zd)), N (the length of the extended
%           decvartable), nnzB and ncell (the nonzeros and the number of
%           the B parameters of P), build_s (seconds, the declaration of
%           P), t_own and t_pad (seconds), peak_MB and alloc_MB (PeakMem
%           and TotalMemAllocated, MB), in_MB (the memory of the B
%           parameters and of P.Zd, MB); and 'table' field legacy, columns
%           q_dopvar and q_cdopvar (the numbers of decision variables),
%           t_sop, t_getsol_lpivar, t_dpvar_names and t_lpivar_build
%           (seconds, NaN when skipped); one printed line per row.
%
% NOTES
% Two results of one size that differ, or lpigetsol values other than the
% last q entries of RRx, are an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - bench_getsol_sop
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
% Initial coding MMP, 09/29/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was a measured run time
%                and memory, the reasons no solve is needed and the
%                profiled call is separate, and references to CLAUDE.md and
%                to another test.

% Defaults for the fields OPTS lacks
dflt = struct('q',[1e4 1e5 1e6 3e6],'nv',[1 2 3],'deg',1,'reps',3, ...
              'gram_q',[1e4 1e5 1e6],'gram_nv',1,'legacy_q',[1e3 1e4 1e5 1e6 3e6],'legacy_cap',120,'multi',true);
if nargin<1,    opts = struct();    end
fn = fieldnames(dflt);
for k = 1:numel(fn),    if ~isfield(opts,fn{k}),    opts.(fn{k}) = dflt.(fn{k});    end,    end
% Variable names (nv <= 3); the output
V = {'s1','s2','s3'};
T = struct();

% % % 1. lpivar_cdopvar sweep over (nv, q), one square block.
rows = {};
for nv = opts.nv
    % Decision variables at one component
    vars = V(1:nv);     dom = repmat([0 1],nv,1);
    q1 = count_q(vars,dom,1,opts.deg);
    for qt = opts.q
        % Component count for about qt variables; P declared (timed);
        % getsol measured
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        t0 = tic;   [prog,P] = lpivar_cdopvar(prog,m,{vars},dom,opts.deg);  tb = toc(t0);
        rows(end+1,:) = measure(sprintf('lpivar nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
end
if ~isempty(rows),  T.lpivar = cell2table(rows,'VariableNames',hdr());  end

% % % 2. Multi-block: R^m x L2^m[vars], 2x2 blocks on one list.
if opts.multi
    rows = {};
    for nv = opts.nv
        % Decision variables at one component per space; the component
        % count for the largest q; P declared (timed); getsol measured
        vars = V(1:nv);     dom = repmat([0 1],nv,1);
        q1 = count_q_multi(vars,dom,1,opts.deg);
        qt = opts.q(end);
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        sp = {{},vars};
        t0 = tic;
        [prog,P] = lpivar_cdopvar(prog,struct('out',[m;m],'in',[m;m]), ...
                                  struct('out',{sp},'in',{sp}),dom,opts.deg);
        tb = toc(t0);
        rows(end+1,:) = measure(sprintf('2x2 multi nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
    T.multi = cell2table(rows,'VariableNames',hdr());
end

% % % 3. Gram form (poscopvar): several nonzeros per variable, names
% % % listed twice in decvartable. q ~ m^2 for a Gram of m*(basis) rows.
rows = {};
for nv = opts.gram_nv
    % Decision variables of the Gram form at one component
    vars = V(1:nv);     dom = repmat([0 1],nv,1);
    if isempty(opts.gram_q),    break,  end
    prog = mkprog(vars,dom);
    [~,P1] = poscopvar(prog,1,{vars},dom,opts.deg);    q1 = numel(P1.Zd);
    for qt = opts.gram_q
        % Component count for about qt variables; P declared (timed);
        % getsol measured
        m = max(1,round(sqrt(qt/q1)));
        prog = mkprog(vars,dom);
        t0 = tic;   [prog,P] = poscopvar(prog,m,{vars},dom,opts.deg);   tb = toc(t0);
        rows(end+1,:) = measure(sprintf('poscopvar nv=%d',nv),nv,prog,P,tb,opts.reps); %#ok<AGROW>
        print_row(rows(end,:));
        clear prog P
    end
end
if ~isempty(rows),  T.gram = cell2table(rows,'VariableNames',hdr());    end

% % % 4. Legacy paths, 1-D. A route stops where its projected cost at the
% % % next size exceeds legacy_cap: lpivar build and getsol_lpivar ~q,
% % % the dpvar-of-names route ~q^2.25 (measured 1e3..3e4, 09/29).
rows = {};  skipA = false;  skipB = false;
for iq = 1:numel(opts.legacy_q)
    % Ratio of the next size to this one (inf at the last); the component
    % count and the 1-D program
    qt = opts.legacy_q(iq);
    rq = inf;   if iq<numel(opts.legacy_q),  rq = opts.legacy_q(iq+1)/qt;   end
    q1 = count_q({'s1'},[0 1],1,opts.deg);
    m = max(1,round(sqrt(qt/q1)));
    pvar s1 s1_dum
    progL = lpiprogram(s1,s1_dum,[0 1]);
    d3 = opts.deg*[1 1 1];
    % The legacy 'dopvar' (route A, unless skipped), then the 'cdopvar' of
    % the same family; random values
    tA = NaN;   tB = NaN;   tbL = NaN;  qd = NaN;
    if ~skipA
        t0 = tic;   [progL,Pd] = lpivar(progL,[0 0;m m],d3);    tbL = toc(t0);  % legacy dopvar
        qd = numel(progL.decvartable);
    end
    [progL,Pc] = lpivar_cdopvar(progL,m,{{'s1'}},[0 1],d3);   % same family, own variables
    NL = numel(progL.decvartable);
    progL.solinfo.RRx = randn(NL,1);    progL.solinfo.info = struct('pinf',0);
    qc = numel(Pc.Zd);
    % getsol_lpivar_sop; getsol_lpivar (route A) and lpigetsol on the names
    % (route B), each skipped from the size its projected time exceeds
    % legacy_cap
    t0 = tic;   Sc = getsol_lpivar_sop(progL,Pc);   tn = toc(t0); %#ok<NASGU>
    if ~skipA
        t0 = tic;   Sd = getsol_lpivar(progL,Pd);   tA = toc(t0); %#ok<NASGU>
        skipA = rq*max(tA,tbL) > opts.legacy_cap;
    end
    if ~skipB
        t0 = tic;   dv = double(lpigetsol(progL,dpvar(cellstr(Pc.Zd(:)))));   tB = toc(t0);
        skipB = tB*rq^2.25 > opts.legacy_cap;
        assert(isequal(dv(:),progL.solinfo.RRx(end-qc+1:end)),'recover_dvals route read other values')
    end
    % Record and print
    rows(end+1,:) = {qd,qc,tn,tA,tB,tbL}; %#ok<AGROW>
    fprintf('legacy q=%8d (cdopvar q=%8d): getsol_lpivar_sop %8.4f s | getsol_lpivar(dopvar) %8.3f s | lpigetsol(dpvar of q names) %8.3f s | lpivar build %7.2f s\n',...
        qd,qc,tn,tA,tB,tbL);
    clear progL Pd Pc Sd Sc
end
if ~isempty(rows),  T.legacy = cell2table(rows,'VariableNames',{'q_dopvar','q_cdopvar','t_sop','t_getsol_lpivar','t_dpvar_names','t_lpivar_build'});  end
end


% ------------------------------------------------------------------------
function h = hdr()
% H = HDR() returns the column names of the tables lpivar, multi and gram.
% INPUT  none.
% OUTPUT h: 1x12 'cellstr'.
h = {'family','nv','q','N','nnzB','ncell','build_s','t_own','t_pad','peak_MB','alloc_MB','in_MB'};
end

function print_row(r)
% PRINT_ROW(R) prints one row of the tables lpivar, multi and gram.
% INPUT  r: 1x12 'cell', the entries of the columns of HDR.
% OUTPUT none; one printed line.
fprintf('%-16s nv=%d q=%9d N=%9d nnzB=%10d cells=%4d | build %7.1f s | getsol own %8.4f s, pad %8.4f s | peak %8.1f MB, alloc %8.1f MB (B+Zd %8.1f MB)\n',r{:});
end

function row = measure(lbl,nv,prog,P,tb,reps)
% ROW = MEASURE(LBL,NV,PROG,P,TB,REPS) times getsol_lpivar_sop(PROG,P) at
% random values drawn per name (a name listed twice has one value), on
% the decvartable of PROG ('own') and on that table extended by 2q names
% P lacks and shuffled ('pad'), each the least of REPS wall times; the two
% results must be equal; one more profiled 'pad' call gives PeakMem and
% TotalMemAllocated.
% INPUT  lbl: 'char', the row label; nv: 'double', the number of spatial
%        variables; prog: 'struct', the LPI program that declares P; P:
%        'cdopvar' object; tb: 'double', the declaration time (s); reps:
%        'double'.
% OUTPUT row: 1x12 'cell', the entries of the columns of HDR.
% Random values per name; the extended, shuffled table
dt = prog.decvartable(:);   q = numel(P.Zd);
[u,~,j] = unique(dt);   val = randn(numel(u),1);
prog.solinfo.RRx = val(j);      prog.solinfo.info = struct('pinf',0);
pad = strcat('zz_pad_',cellstr(string((1:2*q)')));
dt2 = [dt; pad];    rr2 = [val(j); randn(2*q,1)];   pm = randperm(numel(dt2));
prog2 = prog;   prog2.decvartable = dt2(pm);    prog2.solinfo.RRx = rr2(pm);
N = numel(dt2);
[nnzB,ncell,inMB] = op_size(P);
% Least of REPS wall times on each table; equal results
tc = inf;   ts = inf;
for r = 1:reps
    t0 = tic;   S1 = getsol_lpivar_sop(prog,P);    tc = min(tc,toc(t0));
    t0 = tic;   S2 = getsol_lpivar_sop(prog2,P);   ts = min(ts,toc(t0));
end
assert(isequal(S1,S2),'%s: own and padded/shuffled tables disagree',lbl)
clear S1 S2
% One profiled call on the extended table: PeakMem and TotalMemAllocated
profile('-memory','on');
S2 = getsol_lpivar_sop(prog2,P); %#ok<NASGU>
profile('off');
info = profile('info');     ft = info.FunctionTable;
k = find(strcmp({ft.FunctionName},'getsol_lpivar_sop'),1);
pk = ft(k).PeakMem/2^20;    al = ft(k).TotalMemAllocated/2^20;
profile('clear');
row = {lbl,nv,q,N,nnzB,ncell,tb,tc,ts,pk,al,inMB};
end

function [nnzB,ncell,MB] = op_size(P)
% [NNZB,NCELL,MB] = OP_SIZE(P) returns the nonzeros and the number of the
% B parameters of the 'sdopvar' blocks of P, and the memory of those
% parameters and of P.Zd.
% INPUT  P: 'cdopvar' object.
% OUTPUT nnzB, ncell: 'double'; MB: 'double', in MB.
nnzB = 0;   ncell = 0;  by = 0;
for k = 1:numel(P.C)
    b = P.C{k};
    if ~isa(b,'sdopvar'),   continue,   end
    prm = b.params;
    for g = 1:numel(prm.B)
        nnzB = nnzB + nnz(prm.B{g});    ncell = ncell + 1;
        Bg = prm.B{g};  w = whos('Bg');     by = by + w.bytes;
    end
end
Zd = P.Zd;  w = whos('Zd');     MB = (by + w.bytes)/2^20;
end

function q = count_q(vars,dom,m,deg)
% Q = COUNT_Q(VARS,DOM,M,DEG) returns the number of decision variables of
% the square operator on L2^M[VARS] that lpivar_cdopvar declares at degree
% DEG.
% INPUT  vars: 'cellstr', the variables; dom: nv x 2 'double'; m, deg:
%        'double'.
% OUTPUT q: 'double'.
prog = mkprog(vars,dom);
[~,P] = lpivar_cdopvar(prog,m,{vars},dom,deg);
q = numel(P.Zd);
end

function q = count_q_multi(vars,dom,m,deg)
% Q = COUNT_Q_MULTI(VARS,DOM,M,DEG) returns the number of decision
% variables of the square operator on R^M x L2^M[VARS] that lpivar_cdopvar
% declares at degree DEG.
% INPUT  vars: 'cellstr', the variables; dom: nv x 2 'double'; m, deg:
%        'double'.
% OUTPUT q: 'double'.
prog = mkprog(vars,dom);    sp = {{},vars};
[~,P] = lpivar_cdopvar(prog,struct('out',[m;m],'in',[m;m]),struct('out',{sp},'in',{sp}),dom,deg);
q = numel(P.Zd);
end

function prog = mkprog(vars,dom)
% PROG = MKPROG(VARS,DOM) returns an LPI program on the spatial variables
% VARS with the dummy variables VARS_dum and domain DOM: from lpiprogram
% for at most 2 variables, otherwise from sosprogram with that vartable
% and dom.
% INPUT  vars: 'cellstr', the variables; dom: nv x 2 'double'.
% OUTPUT prog: 'struct'.
nv = numel(vars);
if nv<=2
    prog = lpiprogram(polynomial(vars(:)),[],dom);
else
    prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));
    prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))];
    prog.dom = dom;
end
end
