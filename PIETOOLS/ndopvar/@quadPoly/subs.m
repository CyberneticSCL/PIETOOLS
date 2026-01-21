function H = subs(F, vars, vals)

% vars: char/string or cellstr of variable names like 's1','theta2','t3'
% vals: numeric 0/1 or char/string or cellstr like 0,1,'s2','theta2'

nvar = size(F.Zs,2);

[modeS, toS, modeT, toT] = buildMaps(nvar, vars, vals);

m  = F.dim(1);
n  = F.dim(2);
ds = size(F.Zs,1);
dt = size(F.Zt,1);

[I,J,V] = find(F.C);
if isempty(V)
    H = quadPoly(sparse(m*1, n*1), zeros(1,nvar), zeros(1,nvar), F.dim, F.ns, F.nt);
    return;
end

br = floor((I-1)/ds) + 1;   ls = I - (br-1)*ds;
bc = floor((J-1)/dt) + 1;   lt = J - (bc-1)*dt;

ZsNew = zeros(0,nvar);
ZtNew = zeros(0,nvar);
mapS  = containers.Map('KeyType','char','ValueType','double');
mapT  = containers.Map('KeyType','char','ValueType','double');

brH = zeros(numel(V),1);
bcH = zeros(numel(V),1);
isH = zeros(numel(V),1);
itH = zeros(numel(V),1);
vH  = zeros(numel(V),1);

kout = 0;

for k = 1:numel(V)
    es = F.Zs(ls(k),:);
    et = F.Zt(lt(k),:);

    [es2, et2, ok] = applyMaps(es, et, modeS, toS, modeT, toT);
    if ~ok
        continue;
    end

    ks = expKey(es2);
    kt = expKey(et2);

    if isKey(mapS, ks)
        is = mapS(ks);
    else
        ZsNew(end+1,:) = es2; %#ok<AGROW>
        is = size(ZsNew,1);
        mapS(ks) = is;
    end

    if isKey(mapT, kt)
        it = mapT(kt);
    else
        ZtNew(end+1,:) = et2; %#ok<AGROW>
        it = size(ZtNew,1);
        mapT(kt) = it;
    end

    kout = kout + 1;
    brH(kout) = br(k);
    bcH(kout) = bc(k);
    isH(kout) = is;
    itH(kout) = it;
    vH(kout)  = V(k);
end

if kout == 0
    ZsNew = zeros(1,nvar);
    ZtNew = zeros(1,nvar);
    CH = sparse(m*size(ZsNew,1), n*size(ZtNew,1));
    H  = quadPoly(CH, ZsNew, ZtNew, F.dim, F.ns, F.nt);
    return;
end

brH = brH(1:kout);
bcH = bcH(1:kout);
isH = isH(1:kout);
itH = itH(1:kout);
vH  = vH(1:kout);

dsN = size(ZsNew,1);
dtN = size(ZtNew,1);

Inew = (brH-1)*dsN + isH;
Jnew = (bcH-1)*dtN + itH;

% Order triplets for faster sparse build
nRows = m*dsN;
lin = Inew + (Jnew-1)*nRows;
[~,p] = sort(lin);
Inew = Inew(p); Jnew = Jnew(p); vH = vH(p);

CH = sparse(Inew, Jnew, vH, m*dsN, n*dtN);
H  = quadPoly(CH, ZsNew, ZtNew, F.dim, F.ns, F.nt);

end

% ---------------- local helpers ----------------

function [modeS,toS,modeT,toT] = buildMaps(nvar, vars, vals)
% mode: 0 -> set to 0, 1 -> set to 1, 2 -> map to s(to), 3 -> map to theta(to)

% defaults: identity
modeS = 2*ones(1,nvar); toS = 1:nvar;      % s_i -> s_i
modeT = 3*ones(1,nvar); toT = 1:nvar;      % theta_i -> theta_i

if nargin < 2 || isempty(vars)
    return;
end

[vars, vals] = normalizePairs(vars, vals);

for k = 1:numel(vars)
    [side, idx] = parseVar(vars{k});
    [rMode, rSide, rIdx] = parseRhs(vals{k});

    if side == 's'
        if rMode <= 1
            modeS(idx) = rMode; toS(idx) = idx;
        elseif rSide == 's'
            modeS(idx) = 2; toS(idx) = rIdx;
        else
            modeS(idx) = 3; toS(idx) = rIdx;
        end
    else % theta side
        if rMode <= 1
            modeT(idx) = rMode; toT(idx) = idx;
        elseif rSide == 's'
            modeT(idx) = 2; toT(idx) = rIdx;   % theta_i -> s_rIdx
        else
            modeT(idx) = 3; toT(idx) = rIdx;   % theta_i -> theta_rIdx
        end
    end
end
end

function [vars, vals] = normalizePairs(vars, vals)
if ischar(vars) || isstring(vars), vars = {char(vars)}; end
if ischar(vals) || isstring(vals) || isnumeric(vals) || islogical(vals)
    vals = {vals};
end
if isstring(vars), vars = cellstr(vars); end
if isstring(vals), vals = cellstr(vals); end
if ~iscell(vars) || ~iscell(vals) || numel(vars) ~= numel(vals)
    error('subs: vars and vals must be same-length scalars/cells.');
end
for i = 1:numel(vars)
    vars{i} = char(string(vars{i}));
end
end

function [side, idx] = parseVar(v)
v = lower(strtrim(char(string(v))));
v = strrep(v, ' ', '');

tok = regexp(v, '^(s|t|theta)(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('subs: bad variable name "%s" (use s1, theta2, t3).', v);
end
tag = tok{1};
idx = str2double(tok{2});
if strcmp(tag,'s')
    side = 's';
else
    side = 't'; % theta side
end
end

function [rMode, rSide, rIdx] = parseRhs(x)
% rMode: 0 or 1 for constants; otherwise 2 (maps to s) or 3 (maps to theta)
% rSide: 's' or 't' if mapping to variable; rIdx target index
rSide = 's'; rIdx = 1;

if isnumeric(x) || islogical(x)
    rMode = (x~=0);  % 0 or 1
    return;
end

s = lower(strtrim(char(string(x))));
s = strrep(s, ' ', '');

if strcmp(s,'0')
    rMode = 0; return;
elseif strcmp(s,'1')
    rMode = 1; return;
end

tok = regexp(s, '^(s|t|theta)(\d+)$', 'tokens', 'once');
if isempty(tok)
    error('subs: RHS must be 0,1,s#,theta#,t#; got "%s".', s);
end
tag = tok{1};
rIdx = str2double(tok{2});
if strcmp(tag,'s')
    rMode = 2; rSide = 's';
else
    rMode = 3; rSide = 't';
end
end

function [es2, et2, ok] = applyMaps(es, et, modeS, toS, modeT, toT)
nvar = numel(es);
es2 = zeros(1,nvar);
et2 = zeros(1,nvar);
ok  = true;

for i = 1:nvar
    a = es(i);
    if a ~= 0
        ms = modeS(i);
        if ms == 0
            ok = false; return;          % s_i -> 0 kills term
        elseif ms == 1
            % s_i -> 1 : drop exponent
        elseif ms == 2
            es2(toS(i)) = es2(toS(i)) + a;  % s_i -> s_j
        else
            et2(toS(i)) = et2(toS(i)) + a;  % s_i -> theta_j
        end
    end

    b = et(i);
    if b ~= 0
        mt = modeT(i);
        if mt == 0
            ok = false; return;          % theta_i -> 0 kills term
        elseif mt == 1
            % theta_i -> 1 : drop exponent
        elseif mt == 2
            es2(toT(i)) = es2(toT(i)) + b;  % theta_i -> s_j
        else
            et2(toT(i)) = et2(toT(i)) + b;  % theta_i -> theta_j
        end
    end
end
end

function k = expKey(e)
k = sprintf('%d,', e);
end
