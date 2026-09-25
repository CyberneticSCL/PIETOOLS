%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSCOPVAR_VS_POSLPIVAR cross-checks 'poscopvar' against 'poslpivar'
% on the MIXED space R^n x L_2^m[a,b].
%
% This is the case 'possopvar' does not cover, so 'test_possopvar_vs_poslpivar'
% has to exclude the R^n block to make the comparison. Here nothing is
% excluded: the two functions parameterize the same family of positive 4-PI
% operators, and all four components are compared at once.
%
%   poslpivar   Pop = Zop'*(g*Q)*Zop with
%               (Zop*[u1;u2])(s) = [u1; Z1(s)u2(s);
%                                   int_a^s Z2(s,th)u2(th)dth;
%                                   int_s^b Z3(s,th)u2(th)dth]
%   poscopvar   Pop = sum_kl sum_ij (Z_{alpha_i^k})* Q (Z_{alpha_j^l}) over
%               the two spaces k,l in {R^n, L_2^m[a,b]}
%
% The correspondence is: the single basis operator of the R^n space is the
% first row of Zop, which is the identity and so has degree 0 in the
% integration variable; and over the L_2 space alpha=1 is Z1, alpha=2 is Z2
% and alpha=3 is Z3, as in the single-space test.
%
% What this exercises that no single-space test can: the one-sided integrals.
% A block pair with the variable on one side only - the R^n-to-L_2 and
% L_2-to-R^n blocks - has no semiseparable integral to do, and its theta is
% eliminated by 'copquadvar's own antiderivative maps rather than by
% 'int_semisep'. The R^n-to-R^n block has no spatial variable on either side
% and is a plain definite integral. Those three paths are new and are what
% this test covers.
%
% As in the single-space test, no correspondence between decision variable
% names or monomial orderings is assumed: each family is linear in its
% decision variables, so each is characterized by the linear span of its
% generators, and the generators are compared by sampling every component on
% a grid and comparing column spans. A permutation or rescaling of Q
% therefore cannot make a wrong operator look right, and cannot make a right
% one look wrong.
%
% MMP, 09/21/2026: Initial coding
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  mopquadvar -> copquadvar, posmopvar -> poscopvar,
%                  posmopvar_generators -> poscopvar_generators,
%                  test_posmopvar_vs_poslpivar -> test_poscopvar_vs_poslpivar.
%                  File was 'test_posmopvar_vs_poslpivar.m'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(11);

% Each case is {label, n (the R^n dimension), m, dom, d1, [d_int,d_mult,d_joint]}
cases = {
    {'n=1, m=1, d1=1, d2=[1 1 2]', 1, 1, [0,1],  1, [1,1,2]}
    {'n=1, m=1, d1=2, d2=[2 1 3]', 1, 1, [0,1],  2, [2,1,3]}
    {'n=2, m=1, d1=1, d2=[1 1 2]', 2, 1, [0,1],  1, [1,1,2]}
    {'n=1, m=2, d1=1, d2=[1 1 2]', 1, 2, [0,1],  1, [1,1,2]}
    {'n=2, m=2, d1=1, d2=[1 1 1]', 2, 2, [0,1],  1, [1,1,1]}
    {'n=1, m=1, shifted domain',   1, 1, [-1,2], 1, [1,1,2]}
    };

% Sample grid on which the kernels are compared.
ngrid = 9;

npass = 0;
for ic = 1:numel(cases)
    [lbl,n,m,dom,d1] = deal(cases{ic}{1:5});
    d2 = cases{ic}{6};

    vars = {'s1'};
    vars_dum = {'s1_dum'};
    pts = grid_points(dom,ngrid);

    % % % poslpivar: the full 4-PI operator, nothing excluded
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    opts_pl = struct('exclude',[0,0,0,0],'psatz',0,'sep',0);
    [~,Pop_pl] = poslpivar(prog,[n,m],{d1,d2,d2},opts_pl);
    Kpl = poslpivar_generators(Pop_pl,vars,vars_dum);

    % % % poscopvar: the R^n basis is the identity, hence degree 0 in the
    % integration variable; the three L_2 bases match Z1, Z2 and Z3
    deg_pm = { struct('int',0), ...
               { struct('int',d1,'mult',0), ...
                 struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
                 struct('int',d2(1),'mult',d2(2),'joint',d2(3)) } };
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    [~,Pop_pm] = poscopvar(prog,[n;m],{{},vars},dom,deg_pm);
    Kpm = poscopvar_generators(Pop_pm,vars);

    % % % Sample both families of generator kernels on the grid
    Phi_pl = sample_generators(Kpl,n,m,vars,vars_dum,pts);
    Phi_pm = sample_generators(Kpm,n,m,vars,vars_dum,pts);

    % % % Compare the column spans
    [r_pl,r_pm,r_both] = compare_spans(Phi_pl,Phi_pm);
    if r_pl~=r_pm || r_both~=r_pl
        error(['test_poscopvar_vs_poslpivar: the two families differ for case '...
               '''%s'': rank(poslpivar) = %d, rank(poscopvar) = %d, '...
               'rank(joint) = %d.'],lbl,r_pl,r_pm,r_both);
    end

    fprintf('  passed: %-28s dim of family = %3d  (generators: %4d poslpivar, %4d poscopvar)\n',...
            lbl,r_both,size(Phi_pl,2),size(Phi_pm,2));
    npass = npass+1;
end

fprintf('poscopvar vs poslpivar test passed (%d of %d cases).\n',npass,numel(cases));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = poslpivar_generators(Pop,vars,vars_dum)
% Extract the kernel of each decision variable of a dopvar returned by
% 'poslpivar', as a 1 x ngen cell of {P,Q1,Q2,R0,R1,R2} tuples. The pieces
% are renamed so that each is a function of the variable the corresponding
% poscopvar block is a function of: Q1 multiplies the INPUT, so its argument
% is the dummy variable, while Q2 and R0 are functions of the output.

comp = {Pop.P, Pop.Q1, Pop.Q2, Pop.R.R0, Pop.R.R1, Pop.R.R2};
tovar = {[], vars_dum, vars, vars, vars, vars};

% Any decision variable that affects the operator appears in one of the six
% pieces. Q11 appears only in P and Q12 only in Q1, so collecting from the
% R kernels alone - which is enough when the R^n block is excluded - would
% miss exactly the generators this test exists to check.
dvars = {};
for i = 1:numel(comp)
    dvars = [dvars; cellstr(string(comp{i}.dvarname(:)))];                  %#ok<AGROW>
end
dvars = unique(dvars);
ngen = numel(dvars);

% The operator is homogeneous in the decision variables, so the generator
% for variable k is obtained by evaluating at the k-th unit vector.
for i = 1:numel(comp)
    Ci = polynomial(pi_dpvar_eval(comp{i},dvars,zeros(ngen,1)));
    if ~isempty(Ci.coefficient) && max(abs(full(Ci.coefficient(:))))>1e-12
        error('Expected the components of poslpivar to have no constant term.');
    end
end

K = cell(1,ngen);
for k = 1:ngen
    ek = zeros(ngen,1);     ek(k) = 1;
    Kk = cell(1,numel(comp));
    for i = 1:numel(comp)
        Ci = pi_dpvar_eval(comp{i},dvars,ek);
        if i==5 || i==6
            % Two-variable kernels: primary to s, dummy to s_dum.
            Ci = rename_pair(Ci,Pop.var1,Pop.var2,vars,vars_dum);
        elseif ~isempty(tovar{i})
            Ci = subs(Ci,Pop.var1,polynomial(tovar{i}(1)));
        end
        Kk{i} = Ci;
    end
    K{k} = Kk;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = poscopvar_generators(Pop,vars)
% Extract the kernel of each decision variable of the 2 x 2 'cdopvar', as a
% 1 x ngen cell of {P,Q1,Q2,R0,R1,R2} tuples in the same layout as
% 'poslpivar_generators'.
%
% Block (1,1) has no spatial variable, (1,2) only an input one, (2,1) only an
% output one and (2,2) both, so the four have 1, 1, 1 and 3 parameter cells
% respectively - the cell array is indexed over the SHARED variables of the
% pair, of which only (2,2) has any.

ngen = numel(Pop.Zd);
for i = 1:2
    for j = 1:2
        Aij = Pop.C{i,j}.params.A;
        for q = 1:numel(Aij)
            if ~isempty(Aij{q}) && max(abs(full(Aij{q}(:))))>1e-12
                error('Expected the blocks of poscopvar to have no constant term.');
            end
        end
    end
end

K = cell(1,ngen);
for k = 1:ngen
    ek = zeros(ngen,1);     ek(k) = 1;
    K11 = pi_sdopvar_kernels(Pop.C{1,1},ek);
    K12 = pi_sdopvar_kernels(Pop.C{1,2},ek);
    K21 = pi_sdopvar_kernels(Pop.C{2,1},ek);
    K22 = pi_sdopvar_kernels(Pop.C{2,2},ek);
    if numel(K11)~=1 || numel(K12)~=1 || numel(K21)~=1 || numel(K22)~=3
        error('Unexpected number of parameter cells in a poscopvar block.');
    end
    % 'pi_sdopvar_kernels' names the input variable by appending '_dum', so
    % block (1,2) already comes back as a function of s_dum and (2,1) of s.
    K{k} = {K11{1},K12{1},K21{1},K22{1},K22{2},K22{3}};
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Phi = sample_generators(K,n,m,vars,vars_dum,pts)
% Sample every component of every generator on the grid, stacking the
% results into the columns of Phi. Every piece is sampled on the same
% (s,s_dum) grid; the pieces that depend on only one of the two, or on
% neither, are simply constant along the rest, which costs nothing and keeps
% the two families' columns aligned by construction.

allvars = [vars,vars_dum];
sz = {[n,n],[n,m],[m,n],[m,m],[m,m],[m,m]};
ngen = numel(K);
Phi = [];
for k = 1:ngen
    col = [];
    for i = 1:6
        Ki = K{k}{i};
        if i==4
            % The multiplier kernel is only defined up to moving degrees
            % between s and s_dum, since delta(s-s_dum) identifies the two.
            % Evaluating on the diagonal removes that ambiguity.
            Ki = subs(Ki,polynomial(vars_dum(1)),polynomial(vars(1)));
        end
        V = pi_poly_grid(Ki,allvars,pts);
        if size(V,2)~=prod(sz{i})
            error('Unexpected kernel dimensions in component %d.',i);
        end
        col = [col;V(:)];                                                   %#ok<AGROW>
    end
    if isempty(Phi)
        Phi = zeros(numel(col),ngen);
    end
    Phi(:,k) = col;                                                         %#ok<AGROW>
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r1,r2,rboth] = compare_spans(A,B)
% Ranks of the two column spans and of their union. Columns are normalized
% first, so that the rank tolerance is not dominated by scaling.

A = normalize_cols(A);
B = normalize_cols(B);
tol = 1e-8*max(size([A,B]));
r1 = rank(A,tol*norm(A));
r2 = rank(B,tol*norm(B));
rboth = rank([A,B],tol*norm([A,B]));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = normalize_cols(A)

nrm = vecnorm(A);
nrm(nrm<eps) = 1;
A = A./nrm;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = rename_pair(P,var1,var2,vars,vars_dum)
% Rename the primary and dummy variables of poslpivar to the given names.

out = P;
if ~isequal(char(var1),char(polynomial(vars(1))))
    out = subs(out,var1,polynomial(vars(1)));
end
if ~isequal(char(var2),char(polynomial(vars_dum(1))))
    out = subs(out,var2,polynomial(vars_dum(1)));
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pts = grid_points(dom,ngrid)
% Tensor grid of sample points in (s,s_dum).

g = linspace(dom(1,1),dom(1,2),ngrid)';
[S,T] = ndgrid(g,g);
pts = [S(:),T(:)];

end
