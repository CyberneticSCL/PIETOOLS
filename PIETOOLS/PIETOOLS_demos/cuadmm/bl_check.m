function R = bl_check(suite,mode)
% bl_check(suite)          run one suite against the banked expectations
% bl_check(suite,'bank')   run it and RECORD expectations for any case not yet
%                          banked (never overwrites an existing expectation)
%
% Lean by design: no dump to disk, no second Mosek solve -- those exist for the
% GPU comparison, not for regression -- so a 1-D suite costs startup plus
% conversion and nothing else.
%
% WHAT IS COMPARED, AND HOW STRICTLY
%   m, K.f, block list, nnz(At)   EXACT. They are properties of the assembled
%                                 program, independent of solver luck, so ANY
%                                 difference means the program changed. This
%                                 is the check that cleared the lpi_eq nnz fix.
%   rel_b                         within 10x of the banked value AND on the
%                                 same side of 1e-4. Mosek varies ~1e-6
%                                 relative run to run (measured), while a real
%                                 change moves rel_b by orders of magnitude, so
%                                 a tighter band would report noise and a
%                                 looser one would miss a degraded certificate.
%
% VERDICTS
%   PASS    structure exact, rel_b in band
%   FAIL    structure changed, or a passing case stopped passing
%   FIXED   a KNOWN FAILURE now passes -- news, and the reason should be found
%   MOVED   a known failure changed but still fails -- also news
%   ERR     did not build or solve
%   NEW     no banked expectation (run with 'bank' to record one)
% A suite of mode 'runs' checks only that each case builds and solves.
%
% CC, 09/26/2026: 'sentinel_trivial' now means rel_b == 1 (the X=0 residual),
%   not rel_b > 0.5, which had mislabelled stab2_rd (rel_b 2.08, Mosek failing,
%   ||x|| = 0.028) as a trivial-point return. The verdict logic treats both
%   sentinel labels alike, so no verdict changes; the label is what a reader
%   uses to tell "no certificate" from "a bad one".

if nargin < 2, mode = 'check'; end
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
EXP  = fullfile(HERE,'bl_expect.tsv');     % banked expectations: a package INPUT, kept in the repo
OUT  = fullfile(cuadmm_outdir(),'baseline',sprintf('check_%s.tsv',suite));
S = bl_suites(suite);
C = bl_cases();  cid = {C.id};
E = readexp(EXP);

fprintf('\n=== %s: %s  (%d cases) ===\n',S.name,S.question,numel(S.ids));
fid = fopen(OUT,'w');
fprintf(fid,'id\tverdict\tm\tKf\tKs\tnnzAt\trel_b\tref_rel_b\tt_wall\tnote\n');
R = struct('id',{},'verdict',{});
t_all = tic;
for k = 1:numel(S.ids)
    id = S.ids{k};
    c = C(strcmp(cid,id));
    if isempty(c)
        v = 'ERR'; note = 'not in bl_cases';  m=NaN;Kf=NaN;Ks='';nz=NaN;rb=NaN;tw=NaN;
    else
        fprintf('  start %-18s ',id);
        try
            t0 = tic;
            evalc('[sol,~] = feval(c.builder,c.args{:});');
            tw = toc(t0);
            Sh = sdpshape(sol);
            Atf=[]; bf=[];
            for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end %#ok<AGROW>
            x = sol.solinfo.RRx(:);
            rb = norm(full(Atf'*x-bf))/max(norm(full(bf)),eps);
            m = Sh.m; Kf = Sh.Kf; Ks = mat2str(Sh.Ks); nz = Sh.nnzAt;
            [v,note] = judge(S.check,id,E,m,Kf,Ks,nz,rb);
        catch ME
            v = 'ERR'; note = regexprep(ME.message,'[\t\r\n]+',' ');
            m=NaN;Kf=NaN;Ks='';nz=NaN;rb=NaN;tw=NaN;
        end
    end
    if strcmp(v,'NEW') && strcmp(mode,'bank') && ~isnan(m)
        bank(EXP,id,m,Kf,Ks,nz,rb,c.kind);
        note = [note ' -- BANKED'];
    end
    ref = NaN; if isKey(E,id), ref = E(id).rel_b; end
    fprintf('%-6s rel_b=%-10.3e t=%5.1fs %s\n',v,rb,tw,note);
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',id,v,num2str(m),num2str(Kf), ...
            Ks,num2str(nz),num2str(rb,'%.6g'),num2str(ref,'%.6g'),num2str(tw,'%.2f'),note);
    R(end+1) = struct('id',id,'verdict',v); %#ok<AGROW>
end
fclose(fid);
vv = {R.verdict};
fprintf('--- %s: %d PASS, %d FAIL, %d FIXED, %d MOVED, %d ERR, %d NEW  (%.0f s) ---\n', ...
    S.name,nnz(strcmp(vv,'PASS')),nnz(strcmp(vv,'FAIL')),nnz(strcmp(vv,'FIXED')), ...
    nnz(strcmp(vv,'MOVED')),nnz(strcmp(vv,'ERR')),nnz(strcmp(vv,'NEW')),toc(t_all));
fprintf('CHECKDONE %s\n',S.name);
end


function [v,note] = judge(mode,id,E,m,Kf,Ks,nz,rb)
note = '';
if strcmp(mode,'runs')
    if rb < 1e-4, v = 'PASS'; else, v = 'FAIL'; note = 'rel_b >= 1e-4'; end
    return
end
if ~isKey(E,id), v = 'NEW'; note = 'no banked expectation'; return; end
e = E(id);
% --- structure: exact
d = {};
if m  ~= e.m,     d{end+1} = sprintf('m %d->%d',e.m,m); end
if Kf ~= e.Kf,    d{end+1} = sprintf('Kf %d->%d',e.Kf,Kf); end
if ~strcmp(Ks,e.Ks), d{end+1} = sprintf('Ks %s->%s',e.Ks,Ks); end
if ~isnan(e.nnzAt) && nz ~= e.nnzAt, d{end+1} = sprintf('nnz %d->%d',e.nnzAt,nz); end
structural = ~isempty(d);
% --- residual: 10x band, same side of 1e-4
okside = (rb < 1e-4) == (e.rel_b < 1e-4);
inband = abs(log10(max(rb,1e-300)) - log10(max(e.rel_b,1e-300))) <= 1;
wasfail = strncmp(e.state,'sentinel',8);
if structural
    v = 'FAIL'; note = ['STRUCTURE CHANGED: ' strjoin(d,', ')];
elseif wasfail && rb < 1e-4
    v = 'FIXED'; note = sprintf('known failure now passes (rel_b %.2e -> %.2e)',e.rel_b,rb);
elseif wasfail && ~inband
    v = 'MOVED'; note = sprintf('known failure moved (rel_b %.2e -> %.2e)',e.rel_b,rb);
elseif ~okside || ~inband
    v = 'FAIL'; note = sprintf('rel_b %.2e -> %.2e',e.rel_b,rb);
else
    v = 'PASS';
    if wasfail, note = 'still fails as documented'; end
end
end


function E = readexp(fn)
E = containers.Map();
if ~exist(fn,'file'), return; end
fid = fopen(fn,'r');  fgetl(fid);
while true
    l = fgetl(fid);  if ~ischar(l), break; end
    p = strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    if numel(p) < 9, continue; end
    E(p{1}) = struct('m',str2double(p{2}),'Kf',str2double(p{3}),'Ks',p{4}, ...
                     'nnzAt',str2double(p{5}),'rel_b',str2double(p{6}), ...
                     'numerr',str2double(p{7}),'kind',p{8},'state',p{9});
end
fclose(fid);
end


function bank(fn,id,m,Kf,Ks,nz,rb,kind)
st = 'pass';
% Trivial = the X=0 residual, rel_b EXACTLY 1 -- the same test bl_run records.
% rel_b > 1 is a solver returning a point worse than zero (stab2_rd, 2.08), not X=0.
%if rb > 0.5, st = 'sentinel_trivial';                                      % CC, 09/26/2026 (was)
%elseif rb > 1e-4, st = 'sentinel_poor'; end                                % CC, 09/26/2026 (was)
if abs(rb-1) <= 1e-6                                                        % CC, 09/26/2026
    st = 'sentinel_trivial';                                                % CC, 09/26/2026
elseif rb > 1e-4                                                            % CC, 09/26/2026
    st = 'sentinel_poor';                                                   % CC, 09/26/2026
end                                                                         % CC, 09/26/2026
fid = fopen(fn,'a');
fprintf(fid,'%s\t%d\t%d\t%s\t%d\t%.6g\t%s\t%s\t%s\n',id,m,Kf,Ks,nz,rb,'NaN',kind,st);
fclose(fid);
end
