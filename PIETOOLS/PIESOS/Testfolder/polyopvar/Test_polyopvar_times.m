function Test_polyopvar_times()
% Test_polyopvar_times.m -- one local function for each case:
%   case_constant_term      product of distributed polynomials, one of
%                           which has a constant (degree-0) term, e.g
%   case_single_nonzero     @polyopvar/times.m for B with a single nonzero
%                           in B.Z (e.g
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_polyopvar_times_constant_term,
%   Test_polyopvar_times_single_nonzero into this test, as cases;

cases = {@case_constant_term, @case_single_nonzero};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_polyopvar_times: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_polyopvar_times: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_constant_term()
% case_constant_term -- product of distributed
% polynomials, one of which has a constant (degree-0) term, e.g.
% lam*(R^2-<u,u>). The constant is stored as a tensopvar without
% variables, so the product needs otimes(intop,tensopvar) and
% otimes(tensopvar,intop).
% The value of the product at polynomial states is compared with that
% computed by direct integration (ground truth independent of the product
% code).
%

    tol = 1e-9;
    pvar s
    dom = [0,1];
    ncases = 0;

    % States a, b (functions of s), functionals built from innerprod
    x = polyopvar({'a';'b'},[1;1],{'s'},dom,true(2,1));
    a_x = polyopvar({'a'},1,{'s'},dom,true);
    b_x = polyopvar({'b'},1,{'s'},dom,true);
    rng(7);
    for icase = 1:3
        a = rand_poly([1,1],s,2,5);         b = rand_poly([1,1],s,2,5);
        vals = struct('a',a,'b',b);
        R2 = 3+icase;
        c1 = 0.5*icase;
        switch icase
            case 1
                lam = innerprod(a_x,a_x);                    lam_val = int(a^2,s,0,1);
                ipg = innerprod(a_x,a_x);                    ipg_val = int(a^2,s,0,1);
            case 2
                lam = innerprod(a_x,b_x)+2*innerprod(b_x,b_x);  lam_val = int(a*b,s,0,1)+2*int(b^2,s,0,1);
                ipg = innerprod(a_x,a_x);                    ipg_val = int(a^2,s,0,1);
            case 3
                lam = c1*innerprod(b_x,b_x);                 lam_val = c1*int(b^2,s,0,1);
                ipg = innerprod(a_x,b_x);                    ipg_val = int(a*b,s,0,1);
        end
        g = R2-ipg;                                          g_val = R2-ipg_val;
        expected = double(lam_val*g_val);
        if abs(expected)<1e-6
            error("case_constant_term (case %d): value is zero, the test is vacuous.",icase)
        end
        % lam*g, g*lam, lam.*g, g.*lam
        prods = {lam*g, g*lam, lam.*g, g.*lam};
        for k=1:numel(prods)
            val = double(eval_V(prods{k},vals));
            if abs(val-expected)>tol*max(1,abs(expected))
                error("case_constant_term failed (case %d, product %d): value %g, expected %g.",icase,k,val,expected)
            end
            ncases = ncases+1;
        end
    end

    % Product of two functionals which both have constant terms: (R2-<a,a>)*(R3-<b,b>)
    a = rand_poly([1,1],s,2,5);             b = rand_poly([1,1],s,2,5);
    vals = struct('a',a,'b',b);
    g1 = 2-innerprod(a_x,a_x);              g2 = 5-innerprod(b_x,b_x);
    expected = double((2-int(a^2,s,0,1))*(5-int(b^2,s,0,1)));
    val = double(eval_V(g1*g2,vals));
    if abs(val-expected)>tol*max(1,abs(expected))
        error("case_constant_term failed (two constant terms): value %g, expected %g.",val,expected)
    end
    ncases = ncases+1;

    fprintf('case_constant_term passed (%d cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_single_nonzero()
% case_single_nonzero -- @polyopvar/times.m for B with a
% single nonzero in B.Z (e.g. a single state x) and A with several
% monomials, e.g. (x+2*x^2).*x: the value is compared with the explicit
% product of the polynomials in the values of the states.
%

    tol = 1e-10;
    pvar s
    polyopvar x(s) y(s)
    xv = 2+3*s;     yv = 1-s^2;
    ncases = 0;

    cases = {   (x+2*(x.*x)).*x,                    (xv+2*xv^2)*xv
                x.*(x+2*(x.*x)),                    xv*(xv+2*xv^2)
                (x+(x.*x)+3*(x.*x.*x)).*x,          (xv+xv^2+3*xv^3)*xv
                (x+y).*x,                           (xv+yv)*xv
                (x+y).*y,                           (xv+yv)*yv
                (x+2*y+x.*y).*y,                    (xv+2*yv+xv*yv)*yv
                (x+y).*(x+y),                       (xv+yv)^2     };
    for k=1:size(cases,1)
        Y = apply_poly(cases{k,1},{xv,yv});
        D = Y-cases{k,2};
        if ~isempty(nonzeros(D.coefficient)) && max(abs(nonzeros(D.coefficient)))>tol
            error("case_single_nonzero failed (%d): residual %.3e.",k,max(abs(nonzeros(D.coefficient))))
        elseif isempty(nonzeros(polynomial(cases{k,2}).coefficient))
            error("case_single_nonzero (%d): reference value is zero; test would be vacuous.",k)
        end
        ncases = ncases+1;
    end

    fprintf('case_single_nonzero passed (%d cases).\n',ncases)
end

%% ---- helpers --------------------------------------------------------------

function Y = eval_V(V,vals)
% Evaluate the functional V (polyopvar of intop kernels, or constant tensopvar) at the
% states vals.(name), polynomials in s. The dummy variables of s are numbered over the
% factors of the monomial, and named s_j if the monomial has more than 1 factor.
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        % Constant term: tensopvar with a sopvar without variables
        Y = Y + op.ops{1}.params{1};
        continue
    end
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        if d>1
            xj.varname(strcmp(xj.varname,'s')) = {['s_',num2str(j)]};
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end
