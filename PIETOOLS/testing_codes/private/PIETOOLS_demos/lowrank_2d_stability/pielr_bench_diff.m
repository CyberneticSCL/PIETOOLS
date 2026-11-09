function D = pielr_bench_diff(Ra,Rb)                                       % CC, 09/23/2026
% PIELR_BENCH_DIFF  Compare two pielr_bench runs, case by case.
%
%   Ra = pielr_bench(struct('tag','baseline'));
%   ... change one thing ...
%   Rb = pielr_bench(struct('tag','with-the-change'));
%   pielr_bench_diff(Ra,Rb)
%
% THE POINT OF THE WHOLE SUITE IS THIS FUNCTION.  A solver change is worth
% shipping only if it moves a column here, and the direction it moves each
% case is what tells you whether it is a real improvement or a trade.  Reading
% two tables side by side does not scale past a handful of cases and quietly
% hides the cases that got worse.
%
% WHAT COUNTS AS BETTER
%   verdict   no -> yes is better;  yes -> no is a REGRESSION, printed as such
%   gamma     smaller is better, since every certified gamma is an upper bound
%   rel/ipm   THE SCORE (CC, 09/23/2026): the low-rank residual over the
%             reference residual on the same program.  Smaller is better, and
%             unlike the bare residual it is comparable ACROSS programs --
%             which the bare residual is not, since the reference itself
%             ranges over three and a half orders of magnitude here.
%   rel       printed for information, never scored on its own.
%   time      reported, never scored: these runs are not timed under control.
%
% Cases present in only one run are listed separately rather than dropped.

na = {Ra.name};  nb = {Rb.name};
all = unique([na nb],'stable');
ta = tagof(Ra);  tb = tagof(Rb);
fprintf('\npielr_bench_diff:  A = ''%s''   B = ''%s''\n',ta,tb);
fprintf('%-22s %-12s %-12s %-20s %-22s %s\n', ...
    'case','A','B','rel/ipm A -> B','gamma A -> B','');
fprintf('%s\n',repmat('-',1,112));

D = struct('name',{},'change',{},'gam_a',{},'gam_b',{},'ratio_a',{},'ratio_b',{});
nreg = 0;  nfix = 0;  ngam = 0;
for i = 1:numel(all)
    ia = find(strcmp(na,all{i}),1);
    ib = find(strcmp(nb,all{i}),1);
    if isempty(ia) || isempty(ib)
        fprintf('%-22s %-12s %-12s   (present in only one run)\n',all{i}, ...
            tern(isempty(ia),'-','present'),tern(isempty(ib),'-','present'));
        continue
    end
    a = Ra(ia);  b = Rb(ib);
    ch = '';
    if ~a.ok && b.ok, ch = 'FIXED';   nfix = nfix+1;
    elseif a.ok && ~b.ok, ch = 'REGRESSION'; nreg = nreg+1;
    end
    gs = '';
    if ~isnan(a.ratio) || ~isnan(b.ratio)
        gs = sprintf('%.4f -> %.4f',a.ratio,b.ratio);
        if ~isnan(a.ratio) && ~isnan(b.ratio)
            if b.ratio < a.ratio - 1e-6, gs = [gs '  tighter']; ngam = ngam+1;
            elseif b.ratio > a.ratio + 1e-6, gs = [gs '  LOOSER'];
            end
        end
    end
    ss = '';
    if isfinite(getf2(a,'score')) || isfinite(getf2(b,'score'))
        ss = sprintf('%.4g -> %.4g',getf2(a,'score'),getf2(b,'score'));
        if isfinite(getf2(a,'score')) && isfinite(getf2(b,'score'))
            if getf2(b,'score') < getf2(a,'score'), ss = [ss ' +'];
            elseif getf2(b,'score') > getf2(a,'score'), ss = [ss ' -'];
            end
        end
    end
    fprintf('%-22s %-12s %-12s %-20s %-22s %s\n',all{i},vd(a),vd(b),ss,gs,ch);
    D(end+1) = struct('name',all{i},'change',ch,'gam_a',a.gam,'gam_b',b.gam, ...
                      'ratio_a',a.ratio,'ratio_b',b.ratio); %#ok<AGROW>
end
fprintf('%s\n',repmat('-',1,112));
fprintf('A certifies %d/%d   B certifies %d/%d   fixed %d   REGRESSIONS %d   gamma tighter in %d\n\n', ...
    sum([Ra.ok]),numel(Ra),sum([Rb.ok]),numel(Rb),nfix,nreg,ngam);
if nreg>0
    fprintf(['A regression means a case that used to certify no longer does.  ' ...
             'That is a hard stop: the change is not an improvement whatever ' ...
             'else it moved.\n\n']);
end
end

function s = vd(r)
if r.ok, s = sprintf('yes %s',mat2str(r.r)); else, s = 'no'; end
if numel(s)>12, s = s(1:12); end
end
function t = tagof(R)
t = '';
if isfield(R,'meta') && ~isempty(R(1).meta), t = R(1).meta.tag; end
if isempty(t), t = '(untagged)'; end
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function v = getf2(s,f)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
