function disp(obj)
%DISP Pretty-print a quadPoly matrix with tensor-decomposed exponent bases.
%
% Representation:
%   F(s,t) = (I_m ⊗ Zs(s)^T) * C * (I_n ⊗ Zt(t))
% where
%   Zs(s) = Zs{1}(s1) ⊗ ... ⊗ Zs{ks}(s_ks),
%   Zt(t) = Zt{1}(t1) ⊗ ... ⊗ Zt{kt}(t_kt),
% with Zs{i} the exponents of variable ns{i}, Zt{j} those of nt{j}.
%
% This function formats each matrix entry F(i,j) as a sum of monomials
% in the left variables (ns) and right variables (nt), and prints the
% matrix in the command window.
%
% INPUTS
% - obj:    'quadPoly' object;
%
% NOTES
% - Coefficients are shown with 2 significant digits. In the name of a
%   right variable, a leading t is shown as the symbol theta; trailing
%   digits of variable names are shown as subscripts.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - disp
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% SS, 01/16/2026: Initial coding;

% Print the formatted matrix
fprintf('%s\n', format(obj));
end

function s = format(obj)
% S = FORMAT(OBJ) returns the display of the 'quadPoly' OBJ, each entry a
% sum of monomials in the left and right variables.
%
% INPUTS
% - obj:    'quadPoly' object of dimension [m,n];
%
% OUTPUTS
% - s:      string, a title line and one line per row of the matrix.
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

% Title line, then one line per row: entries padded to the column widths,
% separated by commas, rows ended by semicolons
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
% INPUTS
% - Zcell:  1 x k cell, each Zcell{q} is an exponent list for variable q;
% - idx:    positive integer, linear index in the Kronecker basis
%           ordering of MATLAB kron (the last variable varies fastest);
%
% OUTPUTS
% - e:      k x 1 exponent vector aligned with Zcell order.

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
% Monomial string builder
% -------------------------------------------------------------------------
function ms = monomStr(e, names, greekTheta)
% MS = MONOMSTR(E,NAMES,GREEKTHETA) returns the string of the monomial
% with exponents E in the variables NAMES (x1, x2, ... if empty), with
% subscripts and superscripts.
%
% INPUTS
% - e:          k x 1 array of exponents;
% - names:      1 x k cellstr of variable names;
% - greekTheta: logical, true to display a leading t of a name as the
%               symbol theta;
%
% OUTPUTS
% - ms:         string, "1" for the constant monomial.
if isempty(names)
    names = arrayfun(@(k) sprintf('x%d',k), 1:numel(e), 'UniformOutput', false);
end
if all(e==0), ms = "1"; return; end

% One factor per variable with a nonzero exponent, a name starting with t
% displayed with theta if greekTheta
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
    % Trailing digits of the variable name as subscripts
    txt = char(v);
    dg = regexp(txt, '\d+$', 'match', 'once');
    if isempty(dg)
        v = string(txt);
    else
        sub = "";   tbl = '₀₁₂₃₄₅₆₇₈₉';
        for ii = 1:numel(dg),  sub = sub + string(tbl(dg(ii)-'0'+1));  end
        v = string(txt(1:end-numel(dg))) + sub;
    end

    if ek == 1
        parts(end+1) = v;
    else
        % The exponent as superscript digits
        tbl2 = '⁰¹²³⁴⁵⁶⁷⁸⁹';   dgs = char(num2str(ek));   sup = "";
        for ii = 1:numel(dgs),  sup = sup + string(tbl2(dgs(ii)-'0'+1));  end
        parts(end+1) = v + sup;
    end
end
ms = strjoin(parts, "·"); % centered dot
end

function ts = termStr(a, ms, mt)
% TS = TERMSTR(A,MS,MT) returns the string of the term a*ms*mt, for the
% coefficient A and the monomial strings MS (left) and MT (right).
%
% INPUTS
% - a:      scalar coefficient;
% - ms, mt: strings of the left and right monomials;
%
% OUTPUTS
% - ts:     string of the term.
% Product of the nonconstant monomials
if ms=="1" && mt=="1"
    varPart = "";
elseif ms=="1"
    varPart = mt;
elseif mt=="1"
    varPart = ms;
else
    varPart = ms + "·" + mt;
end

% Coefficient with 2 significant digits (real and imaginary part if complex)
if isreal(a)
    coef = sprintf('%.2g', a);
else
    coef = sprintf('(%.2g%+.2gi)', real(a), imag(a));
end

% The coefficient +-1 is shown only as a sign, unless the term is constant
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
% EXPR = JOINTERMS(TERMS) joins the term strings TERMS into a sum, with
% the signs of negative terms, skipping the terms "0".
%
% INPUTS
% - terms:  string array of terms;
%
% OUTPUTS
% - expr:   string of the sum.
% Skip zero terms; append a negative term as a subtraction
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
