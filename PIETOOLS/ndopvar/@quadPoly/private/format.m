function s = format(obj)
% obj is a quadPoly scalar

m  = obj.dim(1); n  = obj.dim(2);
ds = size(obj.Zs,1);
dt = size(obj.Zt,1);

E = strings(m,n);
for i = 1:m
    for j = 1:n
        r = (i-1)*ds + (1:ds);
        c = (j-1)*dt + (1:dt);
        B = obj.C(r,c);

        if nnz(B) == 0
            E(i,j) = "0";
        else
            [ii,jj,vv] = find(B);
            terms = strings(numel(vv),1);
            for t = 1:numel(vv)
                ms = monomStr(obj.Zs(ii(t),:), obj.ns, false);
                mt = monomStr(obj.Zt(jj(t),:), obj.nt, true); % theta-side
                terms(t) = termStr(vv(t), ms, mt);
            end
            E(i,j) = joinTerms(terms);
        end
    end
end

% column widths
w = zeros(1,n);
for j = 1:n
    w(j) = max(strlength(E(:,j)));
end

lines = strings(0,1);
lines(end+1) = sprintf('quadPoly  (%dx%d)', m, n);
lines(end+1) = "[";

for i = 1:m
    rowParts = strings(1,n);
    for j = 1:n
        rowParts(j) = pad(E(i,j), w(j), 'right');   % <-- edit here
    end
    rowStr = strjoin(rowParts, " , ");
    if i < m
        lines(end+1) = "  " + rowStr + " ;";
    else
        lines(end+1) = "  " + rowStr;
    end
end

lines(end+1) = "]";
s = strjoin(lines, newline);

end

function ms = monomStr(e, names, greekTheta)
if isempty(names)
    names = arrayfun(@(k) sprintf('x%d',k), 1:numel(e), 'UniformOutput', false);
end
if all(e==0), ms = "1"; return; end

parts = strings(0,1);
for k = 1:numel(e)
    ek = e(k);
    if ek==0, continue; end
    v = string(names{k});
    if greekTheta
        if startsWith(v,"t")
            v = "θ" + extractAfter(v,1);
        elseif startsWith(lower(v),"theta")
            v = "θ" + extractAfter(v,5);
        end
    end
    v = toSubscriptTrailingDigits(v);

    if ek == 1
        parts(end+1) = v;
    else
        parts(end+1) = v + toSuperscriptInt(ek);
    end
end
ms = strjoin(parts, "·"); % centered dot
end

function v = toSubscriptTrailingDigits(v)
txt = char(v);
d = regexp(txt, '\d+$', 'match', 'once');
if isempty(d), v = string(txt); return; end
base = txt(1:end-numel(d));
v = string(base) + digitsToSub(d);
end

function out = digitsToSub(digits)
table = '₀₁₂₃₄₅₆₇₈₉';
digits = char(digits);
out = "";
for i = 1:numel(digits)
    out = out + string(table(digits(i)-'0'+1));
end
end

function out = toSuperscriptInt(k)
table = '⁰¹²³⁴⁵⁶⁷⁸⁹';
digits = char(num2str(k));
out = "";
for i = 1:numel(digits)
    out = out + string(table(digits(i)-'0'+1));
end
end

function ts = termStr(a, ms, mt)
if ms=="1" && mt=="1"
    varPart = "";
elseif ms=="1"
    varPart = mt;
elseif mt=="1"
    varPart = ms;
else
    varPart = ms + "·" + mt;
end

if isreal(a)
    coef = sprintf('%.2g', a);
else
    coef = sprintf('(%.2g%+.2gi)', real(a), imag(a));
end

if varPart==""
    ts = string(coef);
elseif a==1
    ts = varPart;
elseif a==-1
    ts = "-" + varPart;
else
    ts = string(coef) + "·" + varPart;
end
end

function expr = joinTerms(terms)
expr = "";
for t = 1:numel(terms)
    tt = terms(t);
    if tt=="0", continue; end
    if expr==""
        expr = tt;
    else
        if startsWith(tt,"-")
            expr = expr + " - " + extractAfter(tt,1);
        else
            expr = expr + " + " + tt;
        end
    end
end
if expr=="", expr="0"; end
end
