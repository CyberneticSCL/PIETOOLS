%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSSOPVAR_VS_POSLPIVAR cross-checks 'possopvar' against 'poslpivar'.
%
% For a single spatial variable, the two functions parameterize the same
% family of positive PI operators on L_2^m[a,b]:
%
%   poslpivar   Pop = Zop'*(g*Q)*Zop with Zop built from Z1(s), Z2(eta,s)
%               and Z3(eta,s), excluding the R^n block;
%   possopvar   Pop = sum_ij Z_{alpha_i}* Q_ij Z_{alpha_j} with alpha ranging
%               over {multiplier, lower integral, upper integral}.
%
% Since alpha=1 corresponds to Z1, alpha=2 to Z2 and alpha=3 to Z3, matching
% the degrees should make the two families identical. This is checked without
% relying on any correspondence between decision variable names or monomial
% orderings: each family is linear in its decision variables, so each is
% characterized by the linear span of its generator kernels. The test forms,
% for each family, the matrix whose columns are the generator kernels sampled
% on a grid, and verifies that the two column spans coincide.
%
% Note that this is an independent check of the elimination of the
% integration variable: 'poslpivar' performs it with explicit calls to 'int'
% on dpvar objects, whereas 'possopvar' performs it through the semiseparable
% coefficient tables of INT_SEMISEP.
%
% MP, 08/22/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
rng(7);

% Each case is {label, m, dom, d1, [d_int,d_mult,d_joint]}
cases = {
    {'m=1, d1=1, d2=[1 1 2]', 1, [0,1],  1, [1,1,2]}
    {'m=1, d1=2, d2=[2 1 3]', 1, [0,1],  2, [2,1,3]}
    {'m=1, d1=1, d2=[1 1 1]', 1, [0,1],  1, [1,1,1]}
    {'m=2, d1=1, d2=[1 1 2]', 2, [0,1],  1, [1,1,2]}
    {'m=1, shifted domain',   1, [-1,2], 1, [1,1,2]}
    };

% Sample grid on which the kernels are compared.
ngrid = 9;

npass = 0;
for ic = 1:numel(cases)
    [lbl,m,dom,d1] = deal(cases{ic}{1:4});
    d2 = cases{ic}{5};

    vars = {'s1'};
    vars_dum = {'s1_dum'};
    pts = grid_points(dom,ngrid);

    % % % poslpivar: exclude the R^n block, keep Z1, Z2 and Z3
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    opts_pl = struct('exclude',[1,0,0,0],'psatz',0,'sep',0);
    [~,Pop_pl] = poslpivar(prog,[0,m],{d1,d2,d2},opts_pl);
    Kpl = poslpivar_generators(Pop_pl,vars,vars_dum);

    % % % possopvar: alpha=1 matches Z1, alpha=2 matches Z2, alpha=3 matches Z3
    deg_ps = { struct('int',d1,'mult',0), ...
               struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
               struct('int',d2(1),'mult',d2(2),'joint',d2(3)) };
    prog = lpiprogram(polynomial(vars(:)),[],dom);
    [~,Pop_ps] = possopvar(prog,m,vars,dom,deg_ps);
    Kps = possopvar_generators(Pop_ps);

    % % % Sample both families of generator kernels on the grid
    Phi_pl = sample_generators(Kpl,m,vars,vars_dum,pts);
    Phi_ps = sample_generators(Kps,m,vars,vars_dum,pts);

    % % % Compare the column spans
    [r_pl,r_ps,r_both] = compare_spans(Phi_pl,Phi_ps);
    if r_pl~=r_ps || r_both~=r_pl
        error(['test_possopvar_vs_poslpivar: the two families differ for case ''%s'': '...
               'rank(poslpivar) = %d, rank(possopvar) = %d, rank(joint) = %d.'],...
               lbl,r_pl,r_ps,r_both);
    end

    fprintf('  passed: %-24s dim of family = %3d  (generators: %4d poslpivar, %4d possopvar)\n',...
            lbl,r_both,size(Phi_pl,2),size(Phi_ps,2));
    npass = npass+1;
end

fprintf('possopvar vs poslpivar test passed (%d of %d cases).\n',npass,numel(cases));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = poslpivar_generators(Pop,vars,vars_dum)
% Extract the kernel of each decision variable of a dopvar returned by
% 'poslpivar', as a 1 x ngen cell of {R0,R1,R2} triples.

R = {Pop.R.R0,Pop.R.R1,Pop.R.R2};

% Any decision variable that affects the operator appears in one of the
% three kernels.
dvars = {};
for i = 1:3
    dvars = [dvars; cellstr(string(R{i}.dvarname(:)))];                     %#ok<AGROW>
end
dvars = unique(dvars);
ngen = numel(dvars);

% The kernels are homogeneous in the decision variables, so the generator
% for variable k is obtained by evaluating at the k-th unit vector.
for i = 1:3
    R0i = polynomial(pi_dpvar_eval(R{i},dvars,zeros(ngen,1)));
    if ~isempty(R0i.coefficient) && max(abs(full(R0i.coefficient(:))))>1e-12
        error('Expected the kernels of poslpivar to have no constant term.');
    end
end

% Rename the variables of poslpivar to the naming used elsewhere.
K = cell(1,ngen);
for k = 1:ngen
    ek = zeros(ngen,1);     ek(k) = 1;
    Kk = cell(1,3);
    for i = 1:3
        Kk{i} = rename(pi_dpvar_eval(R{i},dvars,ek),Pop.var1,Pop.var2,vars,vars_dum);
    end
    K{k} = Kk;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function K = possopvar_generators(Pop)
% Extract the kernel of each decision variable of an sdopvar, as a 1 x ngen
% cell of {K_delta,K_lower,K_upper} triples.

ngen = numel(Pop.Zd);
if any(any(abs(cell2mat(reshape(Pop.params.A,1,[])))>1e-12))
    error('Expected the kernels of possopvar to have no constant term.');
end

K = cell(1,ngen);
for k = 1:ngen
    ek = zeros(ngen,1);     ek(k) = 1;
    Kc = pi_sdopvar_kernels(Pop,ek);
    K{k} = {Kc{1},Kc{2},Kc{3}};
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Phi = sample_generators(K,m,vars,vars_dum,pts)
% Sample every generator kernel on the grid, stacking the results into the
% columns of Phi.

allvars = [vars,vars_dum];
ngen = numel(K);
Phi = [];
for k = 1:ngen
    col = [];
    for i = 1:3
        Ki = K{k}{i};
        if i==1
            % The multiplier kernel is only defined up to moving degrees
            % between s and s_dum, since delta(s-s_dum) identifies the two.
            % Evaluating on the diagonal removes that ambiguity.
            Ki = subs(Ki,polynomial(vars_dum(1)),polynomial(vars(1)));
        end
        V = pi_poly_grid(Ki,allvars,pts);
        if size(V,2)~=m*m
            error('Unexpected kernel dimensions.');
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
function out = rename(P,var1,var2,vars,vars_dum)
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
