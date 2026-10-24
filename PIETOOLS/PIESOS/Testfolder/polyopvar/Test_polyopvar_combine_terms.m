function Test_polyopvar_combine_terms()
% Test_polyopvar_combine_terms.m -- one local function for each case:
%   case_basic              verify combine_terms.m merges
%                           duplicate-monomial rows (summing their
%                           coefficients) while preserving the represented
%                           polynomial's value, and shrinks the monomial
%                           basis. The
%   case_zero_class         combine_terms of a polynomial functional (intop
%                           coefficients) in which all operators acting on
%                           the monomials of some degrees are zero (empty),
%                           e.g. a row of
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_polyopvar_combine_terms_zero_class into
%   this test, as cases;

cases = {@case_basic, @case_zero_class};
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
    error('Test_polyopvar_combine_terms: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_polyopvar_combine_terms: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- verify combine_terms.m merges
% duplicate-monomial rows (summing their coefficients) while preserving
% the represented polynomial's value, and shrinks the monomial basis. The
% monomial row is duplicated by hand to force a nontrivial merge; oracle:
% the value of the duplicated polynomial is 2x that of the original.

    tol = 1e-8;

    for trial = 1:20
        rng(4000+trial);

        pvar s1 s2
        x = polyopvar('x',[s1,s2],[0,1]);
        z = x.*x;                 % single-row Z = [1,1]
        if size(z.Z,1)~=1
            error("Test_polyopvar_combine_terms failed: setup assumption violated (trial=%d).",trial)
        end

        % duplicate the monomial row and its coefficient block
        z2 = z;
        z2.Z = [z.Z;z.Z];
        z2.C.ops = [z.C.ops, z.C.ops];
        z2.C.depmat.in = [z.C.depmat.in; z.C.depmat.in];
        z2.C.dims.in = [z.C.dims.in; z.C.dims.in];
        z2.chkval = false;
        z2.C.chkval = false;

        z3 = combine_terms(z2);
        if size(z3.Z,1)~=1
            error("Test_polyopvar_combine_terms failed: duplicate rows not merged (trial=%d).",trial)
        end
        if ~z3.chkval
            error("Test_polyopvar_combine_terms failed: chkval not set (trial=%d).",trial)
        end

        xv = rand_poly([1,1],{'s1','s2'},2,5);
        Y_orig = apply_poly(z,{xv});
        Y_dup  = apply_poly(z2,{xv});
        Y_merged = apply_poly(z3,{xv});

        err_dup = cleanpoly(Y_dup-2*Y_orig,1e-10);
        if max(max(abs(err_dup.C)))>tol
            error("Test_polyopvar_combine_terms failed: manual duplication setup did not read as 2x (trial=%d).",trial)
        end
        % merging must not change the value: z3 equals z2 (2*x^2), not the original x^2
        err_merged = cleanpoly(Y_merged-Y_dup,1e-10);
        if max(max(abs(err_merged.C)))>tol
            error("Test_polyopvar_combine_terms failed: combine_terms changed the polynomial's value (trial=%d).",trial)
        end
    end

    disp('Test_polyopvar_combine_terms passed.')
end

%% --------------------------------------------------------------------------
function case_zero_class()
% case_zero_class -- combine_terms of a
% polynomial functional (intop coefficients) in which all operators
% acting on the monomials of some degrees are zero (empty), e.g. a row of
% the right-hand side of a nonlinear PIE with states x1, x2, in which the
% monomials in x2 do not appear (combine_terms must not apply itself to
% the double 0 left for such a class of monomials). The result must keep
% the other monomials with the same value, and zero (empty) operators for
% the monomials with only zero operators. Oracle: the value of f and of
% combine_terms(f) by eval_terms below (integrating each monomial
% directly).

    tol = 1e-8;
    svars = {'s'};
    dom = struct('in',[0,1],'out',zeros(0,2));
    ncases = 0;

    % Monomials a, a o a, b o b, a o b, with different sets of them zero
    Z = [1,0; 1,1; 2,2; 1,2];
    zero_sets = {2,[2,3],[3,4],[2,3,4],[1,2,3]};
    for iz = 1:numel(zero_sets)
        rng(100+iz);
        C = [];
        ops = cell(1,size(Z,1));
        for k=1:size(Z,1)
            d = nnz(Z(k,:));
            ops{k} = rand_intop([1,1],svars,dom.in,d,2);
            Ck = tensopmat(ops{k});
            if k==1,    C = Ck;     else,   C = [C,Ck];     end %#ok<AGROW>
        end
        f = polyopvar();
        f.Z = Z;                    f.C = C;                    f.varname = {'a';'b'};
        f.dims.in = [1;1];          f.dims.out = 1;
        f.vars.in = svars;          f.vars.out = cell(1,0);
        f.dom.in = dom.in;          f.dom.out = zeros(0,2);
        f.depmat.in = true(2,1);    f.depmat.out = zeros(1,0);
        f = initialize(f);
        % zero the operators of some monomials
        cops = f.C.ops;
        for k=zero_sets{iz}
            cops{1,k} = [];
        end
        f.C.ops = cops;

        g = combine_terms(f);

        % same value before and after
        pvar s
        a = rand_poly([1,1],s,2,5);         b = rand_poly([1,1],s,2,5);
        Yf = eval_terms(f,{'a','b'},{a,b});
        Yg = eval_terms(g,{'a','b'},{a,b});
        if abs(Yf)<1e-8
            error("case_zero_class (set %d): the value is zero, the test is vacuous.",iz)
        elseif abs(Yf-Yg)>tol*max(1,abs(Yf))
            error("case_zero_class failed (set %d): combine_terms changed the value %g to %g.",iz,Yf,Yg)
        end
        % monomials of the zeroed classes keep zero operators
        nz = cellfun(@(o) ~isempty(o),g.C.ops);
        if nnz(nz)>numel(Z(:,1))-numel(zero_sets{iz})
            error("case_zero_class failed (set %d): more nonzero operators than expected.",iz)
        end
        ncases = ncases+1;
    end

    fprintf('case_zero_class passed (%d cases).\n',ncases)
end

%% ---- helpers --------------------------------------------------------------

function Y = eval_terms(f,names,vals)
% Value of f at the states vals{k} (polynomials in s) of states names{k}: sum over the monomials
% of the operators applied to the products of the states, each in its own dummy variable
% (s_1, s_2, ... for a monomial of degree >1, as in apply_op)
Zf = full(f.Z);
Y = 0;
for k=1:size(Zf,1)
    op = f.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    v = 1;
    for j=1:d
        xj = vals{strcmp(names,f.varname{idcs(j)})};
        if d>1
            xj.varname(strcmp(xj.varname,'s')) = {['s_',num2str(j)]};
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
Y = double(polynomial(Y));
end
