function p0_compare(a, b)
% P0_COMPARE(A,B) compares two passes of the Phase-0 fixtures and writes
% ../compare_pass<A>_pass<B>.txt:
% - F1, per example: the largest relative difference of the captured
%   operators (Dop) and of the solution;
% - F2, per PIE: that of the operator outputs (container and opvar);
% - F3, per row: the status, the result and the SDP shape in the note; then
%   the number of rows whose m, or whose nx, ndv or Ks, differ, and every
%   row whose status, result or shape differs;
% - F4, per DEMO: the status, the largest relative difference of its
%   numeric variables, the number of its object variables (copvar,
%   pie_struct, ...) that differ (isequaln), and its numeric scalars that
%   differ.
% Relative difference: max|x-y| / max(1,max|x|), NaN = NaN counted equal;
% Inf when sizes or classes differ. Objects nested in structs or cells are
% not counted. A fixture recorded in one of the two passes only is reported
% as such, not compared.
%
% NOTES
% - Sets the path with P0_PATH first: an object whose class is not on the
%   path loads as a 0x0 double without a warning. Without an export the
%   F4 objects are not compared (reported as NaN).
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
try
    p0_path();  okpath = true;
catch
    okpath = false;
end
rep = fullfile(root,sprintf('compare_pass%d_pass%d.txt',a,b));
fid = fopen(rep,'w');
fprintf(fid,'Phase 0: pass %d vs pass %d, %s\n', a, b, char(datetime('now')));
one = @(F) ~exist(fullfile(root,F,sprintf('pass%d',a)),'dir') || ~exist(fullfile(root,F,sprintf('pass%d',b)),'dir');
% F1
fprintf(fid,'\n== F1 (label, status a/b, Dop, solution)\n');
if one('F1')
    fprintf(fid,'F1: not recorded in both passes\n');
else
    [n, worst] = deal(0, 0);
    for f = dir(fullfile(root,'F1',sprintf('pass%d',a),'*.mat'))'
        g = fullfile(root,'F1',sprintf('pass%d',b),f.name);
        if ~exist(g,'file'),    fprintf(fid,'%s missing in pass %d\n',f.name,b);  continue,   end
        A = load(fullfile(f.folder,f.name));    B = load(g);
        dD = NaN;   dS = NaN;
        if isstruct(A.cap) && isstruct(B.cap),  dD = p0_maxrel(A.cap.Dop, B.cap.Dop);   end
        if isstruct(A.solution) && isstruct(B.solution),    dS = p0_maxrel(A.solution, B.solution);    end
        fprintf(fid,'%s %s/%s %.3g %.3g\n', A.label, A.rec.status, B.rec.status, dD, dS);
        n = n + 1;  worst = max([worst, dD, dS],[],'omitnan');
    end
    fprintf(fid,'F1: %d examples, largest relative difference %.3g\n', n, worst);
end
% F2
fprintf(fid,'\n== F2 (label, status a/b, outputs, opvar outputs)\n');
if one('F2')
    fprintf(fid,'F2: not recorded in both passes\n');
else
    [n, worst] = deal(0, 0);
    for f = dir(fullfile(root,'F2',sprintf('pass%d',a),'*.mat'))'
        g = fullfile(root,'F2',sprintf('pass%d',b),f.name);
        if ~exist(g,'file'),    fprintf(fid,'%s missing in pass %d\n',f.name,b);  continue,   end
        A = load(fullfile(f.folder,f.name));    B = load(g);
        dY = 0;     dO = 0;
        if numel(A.fix)~=numel(B.fix)
            dY = Inf;
        else
            for k = 1:numel(A.fix)
                dY = max(dY, p0_maxrel(A.fix(k).Y, B.fix(k).Y));
                dO = max(dO, p0_maxrel(A.fix(k).Yop, B.fix(k).Yop));
            end
        end
        fprintf(fid,'%s %s/%s %.3g %.3g\n', A.label, A.rec.status, B.rec.status, dY, dO);
        n = n + 1;  worst = max([worst, dY, dO],[],'omitnan');
    end
    fprintf(fid,'F2: %d PIEs, largest relative difference %.3g\n', n, worst);
end
% F3
fprintf(fid,'\n== F3 (row, status a/b, result a/b, shape equal)\n');
fa = fullfile(root,'F3',sprintf('pass%d',a),'sop_results.csv');
fb = fullfile(root,'F3',sprintf('pass%d',b),'sop_results.csv');
if exist(fa,'file') && exist(fb,'file')
    RA = p0_readcsv(fa);    RB = p0_readcsv(fb);
    [ns, nr, nsh, nrow, nm, nx, ma, mb] = deal(0, 0, 0, 0, 0, 0, 0, 0);
    dl = {};
    pat = 'm=\d+ nx=\d+ ndv=\d+ Ks=\[[^\]]*\]|m=\d+ nx=\d+ ndv=\d+ Ks=\d+';
    for k = 1:size(RA,1)
        j = find(strcmp(RB(:,1),RA{k,1}) & strcmp(RB(:,2),RA{k,2}) & strcmp(RB(:,3),RA{k,3}),1);
        if isempty(j),  continue,   end
        nrow = nrow + 1;
        sa = RA{k,4};   sb = RB{j,4};
        ra = str2double(RA{k,5});   rb = str2double(RB{j,5});
        sha = regexp(RA{k,10},pat,'match','once');
        shb = regexp(RB{j,10},pat,'match','once');
        dst = ~strcmp(sa,sb);
        drs = ~(isequaln(ra,rb) || abs(ra-rb) <= 1e-6*max(1,abs(ra)));
        dsh = ~strcmp(sha,shb);
        ns = ns + dst;  nr = nr + drs;  nsh = nsh + dsh;
        % m (equality rows) apart from the rest of the shape (nx, ndv, Ks)
        mka = str2double(regexp(sha,'(?<=^m=)\d+','match','once'));
        mkb = str2double(regexp(shb,'(?<=^m=)\d+','match','once'));
        nm = nm + ~isequaln(mka,mkb);
        nx = nx + ~strcmp(regexprep(sha,'^m=\d+ ',''), regexprep(shb,'^m=\d+ ',''));
        if ~isnan(mka) && ~isnan(mkb),  ma = ma + mka;  mb = mb + mkb;  end
        fprintf(fid,'%s|%s|%s %s/%s %s/%s %d\n', RA{k,1}, RA{k,2}, RA{k,3}, sa, sb, RA{k,5}, RB{j,5}, ~dsh);
        if dst || drs || dsh
            dl{end+1} = sprintf('  %s|%s|%s  %s/%s  %s/%s  [%s] / [%s]', RA{k,1}, RA{k,2}, RA{k,3}, ...
                sa, sb, RA{k,5}, RB{j,5}, sha, shb);                                      %#ok<AGROW>
        end
    end
    fprintf(fid,'F3: %d rows in both; status differs %d, result differs (>1e-6 rel) %d, shape differs %d\n', nrow, ns, nr, nsh);
    fprintf(fid,'F3: m differs on %d rows (sum of m %d / %d); nx, ndv or Ks differs on %d rows\n', nm, ma, mb, nx);
    % Rows recorded in one pass only (a skipped or killed row) are not compared: list them
    ka = strcat(RA(:,1),'|',RA(:,2),'|',RA(:,3));   kb = strcat(RB(:,1),'|',RB(:,2),'|',RB(:,3));
    oa = setdiff(ka,kb);    ob = setdiff(kb,ka);
    fprintf(fid,'F3: pass %d only (%d): %s\nF3: pass %d only (%d): %s\n', ...
        a, numel(oa), strjoin(oa(:)',' '), b, numel(ob), strjoin(ob(:)',' '));
    fprintf(fid,'F3: rows with a different status, result or shape (%d):\n', numel(dl));
    if ~isempty(dl),    fprintf(fid,'%s\n', strjoin(dl, newline));     end
else
    fprintf(fid,'F3: results of one pass missing\n');
end
% F4
fprintf(fid,'\n== F4 (DEMO, status a/b, numeric variables, objects differing/total, scalars a/b differing)\n');
if one('F4')
    fprintf(fid,'F4: not recorded in both passes\n');
else
    [n, worst, nod, ntot] = deal(0, 0, 0, 0);
    for f = dir(fullfile(root,'F4',sprintf('pass%d',a),'DEMO*.mat'))'
        g = fullfile(root,'F4',sprintf('pass%d',b),f.name);
        if ~exist(g,'file'),    fprintf(fid,'%s missing in pass %d\n',f.name,b);  continue,   end
        A = load(fullfile(f.folder,f.name));    B = load(g);
        fa4 = fieldnames(A.vars);   fb4 = fieldnames(B.vars);
        fn = intersect(fa4, fb4);
        [dV, no, nd] = deal(0, 0, 0);
        if numel(fn) < max(numel(fa4),numel(fb4)),  dV = Inf;   end   % variable names differ
        sc = {};
        for j = 1:numel(fn)
            x = A.vars.(fn{j});     y = B.vars.(fn{j});
            if isobject(x) || isobject(y)
                no = no + 1;
                if okpath,  nd = nd + ~isequaln(x,y);   else,   nd = NaN;   end
            else
                dV = max(dV, p0_maxrel(x,y));
                if isnumeric(x) && isscalar(x) && isnumeric(y) && isscalar(y) && ...
                        ~(isequaln(x,y) || abs(x-y) <= 1e-6*max(1,abs(x)))
                    sc{end+1} = sprintf('%s %.10g/%.10g', fn{j}, x, y);                %#ok<AGROW>
                end
            end
        end
        fprintf(fid,'%s %s/%s %.3g %g/%d %s\n', f.name(1:end-4), A.rec.status, B.rec.status, dV, nd, no, strjoin(sc,'; '));
        n = n + 1;  worst = max(worst, dV);     nod = nod + nd;     ntot = ntot + no;
    end
    fprintf(fid,'F4: %d DEMOs, largest relative difference of numeric variables %.3g, object variables differing %g of %d\n', ...
        n, worst, nod, ntot);
end
fclose(fid);
fprintf('\nDONE compare %d %d: %s\n', a, b, rep);
end

function d = p0_maxrel(x, y)
% Largest relative difference over two nested values (see P0_COMPARE).
if (isnumeric(x) || islogical(x)) && (isnumeric(y) || islogical(y))
    if ~isequal(size(x),size(y)),   d = Inf;    return,     end
    if isempty(x),  d = 0;  return,     end
    x = double(full(x(:)));     y = double(full(y(:)));
    nn = isnan(x) & isnan(y);
    if any(isnan(x) ~= isnan(y)),   d = Inf;    return,     end
    d = max(abs(x(~nn)-y(~nn)),[],'omitnan') / max(1, max(abs(x(~nn)),[],'omitnan'));
    if isempty(d),  d = 0;  end
elseif iscell(x) && iscell(y)
    if ~isequal(size(x),size(y)),   d = Inf;    return,     end
    d = 0;
    for k = 1:numel(x),     d = max(d, p0_maxrel(x{k}, y{k}));      end
elseif isstruct(x) && isstruct(y)
    if ~isequal(sort(fieldnames(x)),sort(fieldnames(y))) || numel(x)~=numel(y),     d = Inf;    return,     end
    d = 0;
    fn = fieldnames(x);
    for i = 1:numel(x)
        for k = 1:numel(fn),    d = max(d, p0_maxrel(x(i).(fn{k}), y(i).(fn{k})));  end
    end
else
    d = 0;
    if ~isequaln(x,y),  d = NaN;    end
end
end

function R = p0_readcsv(f)
% Rows of a results CSV (header dropped) as a cell, 10 columns, the note
% keeping any further commas.
L = splitlines(strtrim(fileread(f)));
R = cell(0,10);
for k = 2:numel(L)
    c = strsplit(L{k},',','CollapseDelimiters',false);
    if numel(c) < 10,   continue,   end
    R(end+1,:) = [c(1:9), {strjoin(c(10:end),',')}];                      %#ok<AGROW>
end
end
