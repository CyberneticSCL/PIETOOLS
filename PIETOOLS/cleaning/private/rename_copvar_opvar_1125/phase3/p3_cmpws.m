function p3_cmpws(dold, dnew, outfile)
% P3_CMPWS(DOLD,DNEW,OUTFILE), 11/27/2026: compares the class-free DEMO
% workspaces (d3_ws.m) of the two trees, file by file and leaf by leaf, and
% writes the differing leaves to OUTFILE. A class name stored by d3_ws in a
% field class_ compares under the rename map (copvar = opvar, cdopvar =
% dopvar); fields and variables whose names mark a time (time, cpusec,
% elapsed, toc, tic, date, clock) are skipped. Numeric leaves compare with
% isequaln, so equal means equal bit for bit (NaN = NaN).
fid = fopen(outfile,'w');
F = dir(fullfile(dold,'*.mat'));    ntot = 0;
for k = 1:numel(F)
    g = fullfile(dnew,F(k).name);
    if ~exist(g,'file'),    fprintf(fid,'%s: only in the old run\n',F(k).name);    continue,   end
    a = load(fullfile(dold,F(k).name));     b = load(g);
    d = walk(a,b,'',{});    ntot = ntot + numel(d);
    fprintf(fid,'%s: %d differing leaves\n',F(k).name,numel(d));
    for j = 1:min(numel(d),40),     fprintf(fid,'   %s\n',d{j});   end
end
fprintf(fid,'total differing leaves: %d in %d files\n',ntot,numel(F));
fclose(fid);
type(outfile);
end

function d = walk(x,y,p,d)
% D = WALK(X,Y,P,D): appends to D a line for each leaf of X and Y that differs.
if ischar(x) && ischar(y)
    if ~strcmp(rmap(x),rmap(y)),    d{end+1} = sprintf('%s: ''%s'' vs ''%s''',p,x,y);  end
    return
end
if ~strcmp(class(x),class(y)),  d{end+1} = sprintf('%s: class %s vs %s',p,class(x),class(y));  return,   end
if isstruct(x)
    if ~isequal(size(x),size(y)),   d{end+1} = sprintf('%s: struct size',p);   return,     end
    fx = fieldnames(x);     fy = fieldnames(y);
    if ~isequal(sort(fx),sort(fy)),     d{end+1} = sprintf('%s: fields differ',p);     end
    f = intersect(fx,fy);
    for i = 1:numel(x)
        for j = 1:numel(f)
            if ~isempty(regexpi(f{j},'time|cpusec|elapsed|toc|tic|date|clock','once')),  continue,   end
            d = walk(x(i).(f{j}),y(i).(f{j}),sprintf('%s(%d).%s',p,i,f{j}),d);
        end
    end
elseif iscell(x)
    if ~isequal(size(x),size(y)),   d{end+1} = sprintf('%s: cell size',p);     return,     end
    for i = 1:numel(x),     d = walk(x{i},y{i},sprintf('%s{%d}',p,i),d);    end
elseif isnumeric(x) || islogical(x)
    if ~isequaln(x,y)
        m = Inf;
        if isequal(size(x),size(y)),    m = max(abs(double(full(x(:)))-double(full(y(:)))));  end
        d{end+1} = sprintf('%s: numeric, max|diff| %.3g',p,m);
    end
elseif ~isequaln(x,y)
    d{end+1} = sprintf('%s: %s differs',p,class(x));
end
end

function s = rmap(s)
% S = RMAP(S): a class name under the rename map.
s = regexprep(s,'^cdopvar$','dopvar');      s = regexprep(s,'^copvar$','opvar');
end
