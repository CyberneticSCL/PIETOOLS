% Test_polyopvar_subs_fctr_function.m -- substitution of a factor in a
% polynomial function (tensopvar coefficient with output variables), as
% opposed to Test_polyopvar_subs_fctr.m (functionals, structural check
% only).
%
% Supported cases (see subs_factor in @polyopvar/subs.m):
%   (a) f and the substituted polynomial are both single-factor;
%   (b) f has d factors (pointwise [1,2] or Kronecker [2,2] output), and
%       each monomial of val is degree one, so the factors correspond
%       one-to-one.
% Everything else must raise an informative error.
%
% Verification is semantic and independent of subs: fnew(x,y) is compared
% with Kop applied to the factor list in which factor "fctr" has been
% replaced by the value val(y), evaluated by eval_polyopvar.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;

for Ntot = 1:2
    for d = 1:3
        otypes = {[1,2]};
        if d>1,     otypes = {[1,2],[2,2]};     end
        for io = 1:numel(otypes)
            ot = otypes{io};
            for fctr = 1:d
                for dval = 1:2          % number of monomials in val
                    rng(100*Ntot+10*d+fctr+io);
                    [f,val,Kop] = build_case(Ntot,d,ot,dval);
                    fnew = subs(f,fctr,val);
                    x = rand_pol(svar_all(1:Ntot));
                    y = rand_pol(svar_all(1:Ntot));
                    vy = eval_polyopvar(val,{y});
                    X = repmat({x},1,d);       X{fctr} = vy;
                    Yref = apply_op(Kop,X);
                    xv = cell(numel(fnew.varname),1);
                    for k=1:numel(xv)
                        if strcmp(fnew.varname{k},'x'),     xv{k} = x;
                        else,                               xv{k} = y;      end
                    end
                    Ynew = eval_polyopvar(fnew,xv);
                    % Guard against a vacuous comparison of two zero polynomials
                    Yc = polynomial(Yref).coeff;
                    if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
                        error("Test_polyopvar_subs_fctr_function: reference value is zero (Ntot=%d,d=%d,fctr=%d); test would be vacuous.",Ntot,d,fctr)
                    end
                    err = polynomial(Ynew-Yref);
                    if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
                        error("Test_polyopvar_subs_fctr_function failed (Ntot=%d,d=%d,otype=[%d,%d],fctr=%d,dval=%d): residual %.3e.",...
                            Ntot,d,ot(1),ot(2),fctr,dval,max(abs(full(err.coeff(:)))))
                    end
                    ncases = ncases+1;
                end
            end
        end
    end
end

% % Unsupported cases must error informatively
% (1) single-factor f, val containing a degree-2 monomial
rng(7);
[f,val] = build_case(1,1,[1,2],1);
val2 = build_val_deg2(1);
expect_error(@() subs(f,1,val2),'supported');
% (2) two-factor pointwise f, val with a degree-2 monomial
[f,~] = build_case(1,2,[1,2],1);
expect_error(@() subs(f,1,val2),'supported');
% (3) multi-row f is supported, see Test_polyopvar_subs_multirow

fprintf('Test_polyopvar_subs_fctr_function passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [f,val,Kop] = build_case(N,d,ot,dval)
% Function f(x)(s) = Kop*(x o ... o x), d factors, and val(y) with dval
% degree-one monomials in the state y.
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.C = tensopmat(Kop);       f.varname = {'x'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = svars;          f.vars.out = svars;
f.dom.in = dom.in;          f.dom.out = dom.out;
f.depmat.in = true(1,N);    f.depmat.out = Kop.deps.out;
f = initialize(f);
dep1 = struct('in',true,'out',true);
dep1.in = true(1,N);        dep1.out = true(1,N);
Bop = cell(dval,1);
for k=1:dval
    Bop{k} = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
end
C = Bop{1};
for k=2:dval,   C = [C,Bop{k}];     end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function val = build_val_deg2(N)
% val(y)(s) = (B1 y)(s).*(B2 y)(s), pointwise degree-2 monomial
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(2,N),'out',true(2,N));
Bop = rand_tensopvar(ones(2,2),vars,dom,dep,1,[1,2],2);
val = polyopvar();
val.Z = [1,1];              val.C = tensopmat(Bop);     val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_subs_fctr_function failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("Test_polyopvar_subs_fctr_function failed: expected an error, but none was thrown.")
end
