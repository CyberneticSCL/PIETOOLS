function E = lb_estimate(file, skip)
% E = LB_ESTIMATE(FILE, SKIP) projects the wall time of the full baseline, per
% library and settings tier, from the rows already measured in FILE.
%
% INPUT
%   file:   results file (default results.csv in this folder)
%   skip:   cellstr of 'KEY@tier' left out, as opts.skip of lb_run (default {})
%
% OUTPUT
%   E:  struct with fields 'tiers', 'libs', 'hours' (libs x tiers, 1-D
%       rows), 'solves' (same shape, expected solves), 'hours2d' (2 x
%       tiers, low and high bound for the 2-D rows). Also printed.
%
% NOTES
%   - 1-D rows, two classes: S (stability, gain and norm executives) and
%     K (Hinf/H2 control and estimator synthesis).
%       S: per tier T, log t = a_T + p_T*log(q), least squares over the
%          measured S rows of that tier; q = n1, the number of PIE
%          states in L2 (DDE16 and DDE17 have n0+n1 = 12 and 11 but
%          differ 7x in time, n1 = 9 and 5).
%       K: log t = a_K + p_K*log(q) over the measured K rows on light,
%          scaled to tier T by the S ratio t_S(q,T)/t_S(q,light).
%     t is the time per solve, seconds/n_solves of a row.
%   - Solves per row: measured if the row was run on that tier; else 1
%     for a single solve, and the median measured bisection count for a
%     bisection row.
%   - 2-D rows (lb_run projects each from its lighter tiers and skips it
%     above 600 s per solve, growth 10 per 2-D tier step):
%       low  - every 2-D row costs what PDE 35 costs: measured on light
%              and heavy, light on extreme/stripped, light*10 on custom,
%              heavy and veryheavy skipped (projected > 600 s);
%       high - every 2-D row has a light solve of 60 s: 60 s each on
%              extreme/stripped/light, 600 s on custom, heavy and
%              veryheavy skipped.
%     Neither bound covers a 2-D row whose extreme solve alone exceeds
%     60 s; the first tier of a row is never projected.

here = fileparts(mfilename('fullpath'));
if nargin<1 || isempty(file), file = fullfile(here,'results.csv'); end
if nargin<2, skip = {}; end
S = load(fullfile(here,'catalog.mat'));  C = S.C;
C = C([C.run] & ~strcmp({C.exec},'ERROR'));
tiers = {'extreme','stripped','light','heavy','veryheavy','custom'};
libs = {'PDE','DDE','NDSDDF'};
isK = @(c) ~isempty(regexp(c.type,'control|estimator','once'));
R = lb_read(file);

% Measured 1-D points: tier index, class, q, time per solve
P = zeros(0,4);
for k = 1:numel(R)
    c = C(strcmp({C.key},[R(k).lib,':',R(k).id]));
    if isempty(c) || c.dim>=2 || ~(R(k).n>=1) || strcmp(R(k).status,'error'), continue, end
    if strcmp(R(k).tier,'extreme') && R(k).sec>2 && R(k).n==1 && strcmp(R(k).id,'DDE05/stability_dual')
        continue          % first solve of a MATLAB session: includes warm-up
    end
    P(end+1,:) = [find(strcmp(tiers,R(k).tier)), isK(c), lb_q(c.size), R(k).sec/R(k).n]; %#ok<AGROW>
end
a = nan(1,numel(tiers));  p = nan(1,numel(tiers));
for t = 1:numel(tiers)
    Q = P(P(:,1)==t & P(:,2)==0,:);
    if numel(unique(Q(:,3)))>=2
        b = [ones(size(Q,1),1), log(Q(:,3))] \ log(Q(:,4));
        a(t) = b(1);  p(t) = b(2);
    elseif ~isempty(Q)
        p(t) = 0;  a(t) = mean(log(Q(:,4)));
    end
end
QK = P(P(:,1)==3 & P(:,2)==1,:);
if size(QK,1)>=2 && numel(unique(QK(:,3)))>=2
    bK = [ones(size(QK,1),1), log(QK(:,3))] \ log(QK(:,4));
else
    bK = [NaN; NaN];
end
fprintf('\n1-D time per solve [s] = exp(a + p*log(n1)), measured %d rows:\n', size(P,1));
for t = 1:numel(tiers)
    fprintf('  S %-9s a = %6.2f  p = %5.2f   (n1=1: %.2f s, n1=9: %.2f s)\n', tiers{t}, a(t), p(t), ...
        exp(a(t)), exp(a(t)+p(t)*log(9)));
end
fprintf('  K light     a = %6.2f  p = %5.2f   (n1=2: %.1f s, n1=9: %.1f s)\n', bK(1), bK(2), ...
    exp(bK(1)+bK(2)*log(2)), exp(bK(1)+bK(2)*log(9)));

nb = [R(strcmp({R.kind},'bisect')).n];
nbis = median(nb);
H = zeros(numel(libs),numel(tiers));  N = H;
for t = 1:numel(tiers)
    for k = 1:numel(C)
        c = C(k);
        if c.dim>=2 || ismember([c.key,'@',tiers{t}],skip), continue, end
        L = strcmp(libs,c.lib);
        r = R(strcmp({R.lib},c.lib) & strcmp({R.id},c.key(numel(c.lib)+2:end)) & strcmp({R.tier},tiers{t}));
        if ~isempty(r) && r(1).n>=1
            H(L,t) = H(L,t) + r(1).sec;   N(L,t) = N(L,t) + r(1).n;   continue
        end
        if c.bisect, ns = nbis; else, ns = 1; end
        q = lb_q(c.size);
        ts = exp(a(t) + p(t)*log(q));
        if isK(c)
            ts = max(1, exp(bK(1) + bK(2)*log(q)) / exp(a(3) + p(3)*log(q))) * ts;   % K >= S
        end
        H(L,t) = H(L,t) + ns*ts;   N(L,t) = N(L,t) + ns;
    end
end
H = H/3600;
% 2-D bounds
n2 = nnz([C.dim]>=2);
t35 = struct();
for t = 1:numel(tiers)
    r = R(strcmp({R.lib},'PDE') & strcmp({R.id},'35') & strcmp({R.tier},tiers{t}));
    if ~isempty(r) && r(1).n>=1, t35.(tiers{t}) = r(1).sec; end
end
lo = [t35.light, t35.light, t35.light, 0, 0, 10*t35.light];
hi = [60, 60, 60, 0, 0, 600];
H2 = [lo; hi]*n2/3600;

fprintf('\nProjected hours, 1-D rows (expected solves):\n%-8s', '');
fprintf('%15s', tiers{:});  fprintf('%9s\n','total');
for L = 1:numel(libs)
    fprintf('%-8s', libs{L});
    for t = 1:numel(tiers), fprintf('%8.2f (%4.0f)', H(L,t), N(L,t)); end
    fprintf('%9.2f\n', sum(H(L,:)));
end
fprintf('%-8s', '2-D lo');  fprintf('%15.2f', H2(1,:));  fprintf('%9.2f\n', sum(H2(1,:)));
fprintf('%-8s', '2-D hi');  fprintf('%15.2f', H2(2,:));  fprintf('%9.2f\n', sum(H2(2,:)));
fprintf('1-D total %.2f h; 2-D (%d rows) %.2f to %.2f h; all %.2f to %.2f h\n', sum(H(:)), n2, ...
    sum(H2(1,:)), sum(H2(2,:)), sum(H(:))+sum(H2(1,:)), sum(H(:))+sum(H2(2,:)));
E = struct('tiers',{tiers},'libs',{libs},'hours',H,'solves',N,'hours2d',H2);
end

function q = lb_q(sz)
% n1 of the PIE size string '[n0 n1]', at least 1.
v = str2num(sz); %#ok<ST2NM>
q = max(1,v(min(2,end)));
end

function R = lb_read(file)
% Results rows as a struct array (lib, id, tier, status, n, sec, kind).
R = struct('lib',{},'id',{},'tier',{},'status',{},'n',{},'sec',{},'kind',{});
if ~exist(file,'file'), return, end
L = splitlines(fileread(file));
for k = 2:numel(L)
    f = strsplit(L{k},',','CollapseDelimiters',false);
    if numel(f)<10, continue, end
    kind = 'single';
    if any(strcmp(f{4},{'converged','max_solves','lo_fail','hi_pass'})), kind = 'bisect'; end
    R(end+1) = struct('lib',f{1},'id',f{2},'tier',f{3},'status',f{4}, ...
        'n',str2double(f{7}),'sec',str2double(f{8}),'kind',kind); %#ok<AGROW>
end
end
