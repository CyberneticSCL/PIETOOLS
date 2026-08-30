%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPI_EQ_NDOPVAR verifies 'lpi_eq_ndopvar'.
%
% The main check is against the trusted path: an 'ndopvar' is converted to
% an equivalent 'dopvar' with 'ndopvar2dopvar', the stock 'lpi_eq' is applied
% to that, and the two resulting sets of equality constraints are compared.
% Both define an affine set {r : A*r = b}, assembled from prog.expr.At and
% prog.expr.b, and the test verifies that the two affine sets coincide by
% comparing the row spaces of the augmented matrices [A,b].
%
% This is a representation independent comparison: it does not assume the
% two routines generate the constraints in the same order, nor the same
% number of them, only that they cut out the same set of decision variable
% values. Each case is run twice, once on the sparse operator returned by
% 'rand_ndopvar' and once with the coefficients filled in densely, so that
% the constraints actually have close to full rank.
%
% The 'symmetric' option is checked separately on an operator that really is
% self-adjoint, obtained by running 'poslpivar' and converting the result
% with 'dopvar2ndopvar'. For such an operator the reduced set of constraints
% must cut out exactly the same affine set as the full one. That check does
% not reuse any of the adjoint bookkeeping in 'lpi_eq_ndopvar', so it is not
% circular.
%
% The script additionally checks the handling of empty and zero parameters,
% and the input validation.
%
% MMP, 08/28/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(0);

% Each case is {label, dim, deg, dom, ndec}
cases = {
    {'1D, 1x1, deg 1',      [1,1], 1,     [0,1],       4}
    {'1D, 2x2, deg 1',      [2,2], 1,     [0,1],       6}
    {'1D, 2x3, deg 2',      [2,3], 2,     [0,1],       5}
    {'1D, 1x1, shifted',    [1,1], 1,     [-1,2],      4}
    {'2D, 1x1, deg 1',      [1,1], [1;1], [0,1;0,1],   5}
    {'2D, 2x2, deg 1',      [2,2], [1;1], [0,1;0,1],   4}
    {'2D, 1x2, deg [2;1]',  [1,2], [2;1], [0,1;0,1],   6}
    };

fprintf('=== agreement with lpi_eq on the equivalent dopvar ===\n');
npass = 0;
for ic = 1:numel(cases)
    lbl = cases{ic}{1};
    dim = cases{ic}{2};
    deg = cases{ic}{3};
    dom = cases{ic}{4};
    ndec = cases{ic}{5};
    N = size(dom,1);

    % % % Declare an LPI program with N spatial variables and ndec decision vars
    var1 = polynomial(arrayfun(@(k) sprintf('s%d',k),(1:N)','UniformOutput',false));
    var2 = polynomial(arrayfun(@(k) sprintf('s%d_dum',k),(1:N)','UniformOutput',false));
    dvarname = arrayfun(@(k) sprintf('r%d',k),(1:ndec)','UniformOutput',false);

    prog = lpiprogram(var1,var2,dom);
    prog = lpidecvar(prog,dvarname);

    % % % Build a random operator on those decision variables
    P = rand_ndopvar(dim,deg,dom,var1,var2,dvarname);
    if ~isa(P,'ndopvar')
        error('test_lpi_eq_ndopvar: expected rand_ndopvar to return an ndopvar.');
    end

    [nr_s,rk_s] = check_against_trusted(prog,P,[lbl,' (sparse)']);
    [nr_d,rk_d] = check_against_trusted(prog,densify(P),[lbl,' (dense)']);

    fprintf('  passed: %-20s  sparse %3d rows rank %3d | dense %4d rows rank %4d\n',...
            lbl,nr_s,rk_s,nr_d,rk_d);
    npass = npass+2;
end

% % % The symmetric option, checked on a genuinely self-adjoint operator
fprintf('\n=== symmetric option on a self-adjoint operator ===\n');
for n2 = [1,2]
    for dg = [1,2]
        v1 = polynomial({'s1'});
        v2 = polynomial({'s1_dum'});
        prog = lpiprogram(v1,v2,[0,1]);
        % poslpivar returns a self-adjoint positive operator; excluding the
        % first block leaves a purely L2 to L2 operator, which is what the
        % ndopvar class represents.
        [prog,Pop] = poslpivar(prog,[0,n2],dg,struct('exclude',[1,0,0,0],'psatz',0,'sep',0));
        Pn = dopvar2ndopvar(Pop);
        if ~isa(Pn,'ndopvar')
            error('test_lpi_eq_ndopvar: expected dopvar2ndopvar to return an ndopvar.');
        end

        prog_full = lpi_eq_ndopvar(prog,Pn);
        prog_sym = lpi_eq_ndopvar(prog,Pn,'symmetric');
        [Af,bf] = collect_constraints(prog_full,prog);
        [As,bs] = collect_constraints(prog_sym,prog);

        rf = rank_aug(Af,bf);
        rs = rank_aug(As,bs);
        rj = rank_aug([Af;As],[bf;bs]);
        if ~(rf==rs && rj==rf)
            error(['test_lpi_eq_ndopvar: for a self-adjoint operator the symmetric ',...
                   'option must cut out the same affine set (n2=%d, deg=%d): ',...
                   'rank(full)=%d, rank(sym)=%d, rank(joint)=%d.'],n2,dg,rf,rs,rj);
        end
        fprintf(['  passed: n2=%d deg=%d   full %3d rows rank %3d == symmetric ',...
                 '%3d rows rank %3d\n'],n2,dg,size(Af,1),rf,size(As,1),rs);
        npass = npass+1;
    end
end

% % % Empty and zero parameters must generate no constraints
fprintf('\n=== empty and zero parameters ===\n');
v1 = polynomial({'s1'});
v2 = polynomial({'s1_dum'});
prog = lpiprogram(v1,v2,[0,1]);
prog = lpidecvar(prog,{'r1','r2'});
P = densify(rand_ndopvar([1,1],1,[0,1],v1,v2,{'r1','r2'}));

Pz = P;
for ii = 1:numel(Pz.C)
    Pz.C{ii} = 0*Pz.C{ii};
end
prog_z = lpi_eq_ndopvar(prog,Pz);
if prog_z.expr.num ~= prog.expr.num
    error('test_lpi_eq_ndopvar: an all-zero operator should generate no constraints.');
end

Pe = P;
Pe.C{2} = [];
prog_e = lpi_eq_ndopvar(prog,Pe);
prog_f = lpi_eq_ndopvar(prog,P);
if prog_e.expr.num >= prog_f.expr.num
    error('test_lpi_eq_ndopvar: an empty parameter should generate no constraints.');
end
fprintf('  passed: zero operator adds %d constraint blocks, one empty parameter adds %d of %d\n',...
        prog_z.expr.num-prog.expr.num,prog_e.expr.num-prog.expr.num,...
        prog_f.expr.num-prog.expr.num);
npass = npass+1;

% % % Input validation
fprintf('\n=== input validation ===\n');
Pn0 = rand_ndopvar([1,1],1,[0,1],v1,v2,{});
Pbad = rand_ndopvar([1,1],1,[0,1],v1,v2,{'nope'});
probes = { {'nopvar input',      @() lpi_eq_ndopvar(prog,Pn0)}
           {'polynomial input',  @() lpi_eq_ndopvar(prog,polynomial(1))}
           {'double input',      @() lpi_eq_ndopvar(prog,1)}
           {'unregistered dvar', @() lpi_eq_ndopvar(prog,Pbad)}
           {'bad option',        @() lpi_eq_ndopvar(prog,P,'skew')} };
for k = 1:numel(probes)
    raised = false;
    try
        probes{k}{2}();
    catch ME
        raised = true;
        fprintf('  errors   : %-18s %s\n',probes{k}{1},strtrim(ME.message));
    end
    if ~raised
        error('test_lpi_eq_ndopvar: expected an error for %s.',probes{k}{1});
    end
end

fprintf('\nlpi_eq_ndopvar test passed (%d checks).\n',npass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nrows,rk] = check_against_trusted(prog,P,lbl)
% Compare the constraints generated for the ndopvar P against those the stock
% lpi_eq generates for the equivalent dopvar.

prog_new = lpi_eq_ndopvar(prog,P);
[A1,b1] = collect_constraints(prog_new,prog);

Pd = ndopvar2dopvar(P);
prog_ref = lpi_eq(prog,Pd);
[A2,b2] = collect_constraints(prog_ref,prog);

r1 = rank_aug(A1,b1);
r2 = rank_aug(A2,b2);
rboth = rank_aug([A1;A2],[b1;b2]);
if ~(r1==r2 && rboth==r1)
    error(['test_lpi_eq_ndopvar: constraint sets differ for case %s: ',...
           'rank(new)=%d, rank(ref)=%d, rank(joint)=%d.'],lbl,r1,r2,rboth);
end

nrows = size(A1,1);
rk = r1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = densify(P)
% Replace the sparse random coefficients by dense ones of the same size, so
% that the resulting constraints have close to full rank.

for ii = 1:numel(P.C)
    if isempty(P.C{ii})
        continue
    end
    P.C{ii} = sparse(randn(size(P.C{ii})));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b] = collect_constraints(prog,prog0)
% Assemble the affine system A*r = b from the equality constraints that were
% added to prog0 to obtain prog. SOSTOOLS stores each constraint block as
% At (ndec x nrows) and b (nrows x 1), with the constraint At'*r = b.

A = [];     b = [];
for k = prog0.expr.num+1:prog.expr.num
    if ~strcmpi(prog.expr.type{k},'eq')
        error('Expected only equality constraints.');
    end
    A = [A; full(prog.expr.At{k}).'];                                       %#ok<AGROW>
    b = [b; full(prog.expr.b{k})];                                          %#ok<AGROW>
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rank_aug(A,b)
% Rank of the augmented matrix [A,b], with rows normalized so that the rank
% tolerance is not dominated by row scaling.

M = [A,b];
if isempty(M)
    r = 0;
    return
end
nrm = vecnorm(M,2,2);
nrm(nrm<eps) = 1;
M = M./nrm;
r = rank(M,1e-9*max(size(M)));

end
