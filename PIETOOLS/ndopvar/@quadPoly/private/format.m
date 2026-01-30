function s = format(obj)
%FORMAT Pretty-print a quadPoly matrix with tensor-decomposed exponent bases.
%
% Representation:
%   F(s,t) = (I_m ⊗ Zs(s)^T) * C * (I_n ⊗ Zt(t))
% where
%   Zs(s) = Zs{1}(s1) ⊗ ... ⊗ Zs{ks}(s_ks),  Zs{i} stores exponents for variable ns{i}
%   Zt(t) = Zt{1}(t1) ⊗ ... ⊗ Zt{kt}(t_kt),  Zt{j} stores exponents for variable nt{j}
%
% This function formats each matrix entry F(i,j) as a sum of monomials
% in the left variables (ns) and right variables (nt).

m  = obj.dim(1);
n  = obj.dim(2);

% Tensor basis sizes
ds = prod(cellfun(@numel, obj.Zs));  % # monomials on s-side
dt = prod(cellfun(@numel, obj.Zt));  % # monomials on t-side

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
                % ii(t) is a linear index into the s tensor basis (1..ds)
                % jj(t) is a linear index into the t tensor basis (1..dt)
                es = tensorExpAt(obj.Zs, ii(t));  % exponent vector aligned to obj.ns
                et = tensorExpAt(obj.Zt, jj(t));  % exponent vector aligned to obj.nt

                ms = monomStr(es, obj.ns, false);
                mt = monomStr(et, obj.nt, true);  % theta-side styling (optional)
                terms(t) = termStr(vv(t), ms, mt);
            end

            E(i,j) = joinTerms(terms);
        end
    end
end

% Column widths for aligned display
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
        rowParts(j) = pad(E(i,j), w(j), 'right');
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

% -------------------------------------------------------------------------
% Helper: decode tensor-product basis linear index -> exponent vector
% -------------------------------------------------------------------------
function e = tensorExpAt(Zcell, idx)
%TENSOREXPAT Return exponent vector for tensor basis index idx (1-based).
%
% Zcell : 1×k cell, each Zcell{q} is an exponent list for variable q.
% idx   : linear index in the Kronecker basis ordering consistent with MATLAB kron,
%         i.e., last variable varies fastest.
%
% Output:
%   e : k×1 exponent vector aligned with Zcell order.

k = numel(Zcell);
e = zeros(k,1);

% Sizes of each 1D basis
sz = cellfun(@numel, Zcell);

% MATLAB kron ordering: last factor varies fastest.
% Convert (idx-1) to mixed radix with radices sz(end),...,sz(1).
x = idx - 1;
for q = k:-1:1
    rq = sz(q);
    iq = mod(x, rq) + 1;    % 1..rq index into Zcell{q}
    e(q) = Zcell{q}(iq);
    x = floor(x / rq);
end
end

% -------------------------------------------------------------------------
% Monomial string builder (unchanged, now fed exponent vectors from tensorExpAt)
% -------------------------------------------------------------------------
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
table = '₀₁₂₃₄₅₆₇₈⁹'; %#ok<NASGU>
% NOTE: corrected table below (previous line kept for readability in editors)
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
