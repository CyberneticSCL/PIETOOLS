function out = latex(obj, mode)
if nargin < 2, mode = "string"; end

L = buildLatex(obj);

if mode == "show"
    figure('Color','w'); axis off;
    text(0.01, 0.5, char(L), 'Interpreter','latex', 'FontSize', 14, ...
        'VerticalAlignment','middle', 'HorizontalAlignment','left');
end

out = L;
end

function L = buildLatex(obj)
m  = obj.dim(1); n  = obj.dim(2);
ds = size(obj.Zs,1);
dt = size(obj.Zt,1);

E = strings(m,n);
for i = 1:m
    for j = 1:n
        r = (i-1)*ds + (1:ds);
        c = (j-1)*dt + (1:dt);
        B = obj.C(r,c);

        if nnz(B)==0
            E(i,j) = "0";
        else
            [ii,jj,vv] = find(B);
            terms = strings(numel(vv),1);
            for t = 1:numel(vv)
                ms = monomLatex(obj.Zs(ii(t),:), obj.ns, false);
                mt = monomLatex(obj.Zt(jj(t),:), obj.nt, true);
                terms(t) = termLatex(vv(t), ms, mt);
            end
            E(i,j) = joinTermsLatex(terms);
        end

        % Make sure no actual newline characters end up in entries
        E(i,j) = replace(E(i,j), newline, " ");
        E(i,j) = replace(E(i,j), "\n", " ");
    end
end

rows = strings(m,1);
for i = 1:m
    rows(i) = strjoin(E(i,:), " & ");
end

cols = repmat('c', 1, n);  % n columns
L = "$\left[\begin{array}{" + string(cols) + "}" + ...
    strjoin(rows, " \\ ") + ...
    "\end{array}\right]$";
end

function ms = monomLatex(e, names, isThetaSide)
if isempty(names)
    names = arrayfun(@(k) sprintf('x%d',k), 1:numel(e), 'UniformOutput', false);
end
if all(e==0), ms = "1"; return; end

parts = strings(0,1);
for k = 1:numel(e)
    ek = e(k);
    if ek==0, continue; end

    v = string(names{k});
    [base, sub] = splitNameDigits(v);

    if isThetaSide
        if lower(base)=="t" || lower(base)=="theta"
            base = "\theta";     % <-- single backslash
        end
    end

    if sub ~= ""
        vLatex = base + "_{" + sub + "}";
    else
        vLatex = base;
    end

    if ek==1
        parts(end+1) = vLatex;
    else
        parts(end+1) = vLatex + "^{" + string(ek) + "}";
    end
end

ms = strjoin(parts, " \cdot ");
end

function ts = termLatex(a, ms, mt)
if ms=="1" && mt=="1"
    prod = "";
elseif ms=="1"
    prod = mt;
elseif mt=="1"
    prod = ms;
else
    prod = ms + " \cdot " + mt;
end

if isreal(a)
    coef = string(sprintf('%.6g', a));
else
    coef = "(" + string(sprintf('%.6g', real(a))) + string(sprintf('%+.6g', imag(a))) + "i)";
end

if prod==""
    ts = coef;
elseif a == 1
    ts = prod;
elseif a == -1
    ts = "-" + prod;
else
    ts = coef + " \cdot " + prod;
end
end

function expr = joinTermsLatex(terms)
expr = "";
for t = 1:numel(terms)
    tt = terms(t);
    if tt=="" || tt=="0", continue; end
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

function [base, sub] = splitNameDigits(v)
txt = char(v);
d = regexp(txt, '\d+$', 'match', 'once');
if isempty(d)
    base = string(txt);
    sub  = "";
else
    base = string(txt(1:end-numel(d)));
    sub  = string(d);
end
end
