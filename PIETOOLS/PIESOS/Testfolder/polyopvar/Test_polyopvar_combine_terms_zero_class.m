% Test_polyopvar_combine_terms_zero_class.m -- combine_terms of a
% polynomial functional (intop coefficients) in which all operators
% acting on the monomials of some degrees are zero (empty), e.g. a row of
% the right-hand side of a nonlinear PIE with states x1, x2, in which the
% monomials in x2 do not appear (combine_terms must not apply itself to
% the double 0 left for such a class of monomials). The result must keep
% the other monomials with the same value, and zero (empty) operators for
% the monomials with only zero operators. Oracle: the value of f and of
% combine_terms(f) by eval_terms below (integrating each monomial
% directly).

% DJ (CC), 10/15/2026: Initial coding;

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
        error("Test_polyopvar_combine_terms_zero_class (set %d): the value is zero, the test is vacuous.",iz)
    elseif abs(Yf-Yg)>tol*max(1,abs(Yf))
        error("Test_polyopvar_combine_terms_zero_class failed (set %d): combine_terms changed the value %g to %g.",iz,Yf,Yg)
    end
    % monomials of the zeroed classes keep zero operators
    nz = cellfun(@(o) ~isempty(o),g.C.ops);
    if nnz(nz)>numel(Z(:,1))-numel(zero_sets{iz})
        error("Test_polyopvar_combine_terms_zero_class failed (set %d): more nonzero operators than expected.",iz)
    end
    ncases = ncases+1;
end

fprintf('Test_polyopvar_combine_terms_zero_class passed (%d cases).\n',ncases)


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
