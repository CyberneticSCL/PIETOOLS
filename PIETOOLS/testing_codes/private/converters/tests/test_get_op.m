function results = test_get_op(groups,ntrials)
% TEST_GET_OP Tests of the extraction of PI operators from PDE terms
% (get_op), against the definition of the terms rather than against the
% routines that build the terms from operators.
%
% Each operator P = get_op(F) is applied, block by block, to random
% polynomial states (apply_sopvar), and P*x is compared with F evaluated at
% those states by polynomial integration, with
%   (P x)(s) = R0(s) x(s) + int_a^s R1(s,th) x(th) dth + int_s^b R2(s,th) x(th) dth
% along each variable s shared by the state and the output, a full
% integral along each variable of the state only, and a multiplier along
% each variable of the output only. Kernels have unequal degrees in s and
% the dummy variable th, and domains are not [0,1]. Groups:
%   'S'  semantics: multipliers, partial and full integrals, R and vector
%        components, several rows, inputs w and u, 2D and 3D, the forms of
%        the output space;
%   'E'  unsupported terms, each of which must throw an error whose
%        message names the problem;
%   'N'  new spatial variables: a variable of no state or input is a
%        variable on [0,1] along which the operator is a multiplier;
%   'R'  random operators: states in R or L2 over random subsets of
%        s1,s2,s3 in random order, random output rows, random block types
%        and kernels; checks the output spaces, the values, and
%        <P x, z> = <x, P' z> on the domains of P;
%   'C'  cross-check with the converter: for x_t = F(x), with F random
%        terms into the spaces of the states, the PIE has T = I and A equal
%        to the operator of F, up to the order of the components.
%
% Run as: results = test_get_op;  or  results = test_get_op({'R'},300);
%
% INPUT
% - groups:  (optional) cellstr, the groups to run, from 'S', 'E', 'N',
%            'R', 'C'; all by default;
% - ntrials: (optional) positive integer, the number of random trials of
%            group 'R' (default 100); group 'C' runs a quarter of these;
%
% OUTPUT
% - results: 'struct' array with fields name, passed and message, one
%            element per case or trial.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_get_op
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
% Initial coding DJ - 11/23/2026

if nargin<1 || isempty(groups),  groups = {'S','E','N','R','C'};  end
if nargin<2,  ntrials = 100;  end
pvar s s_dum s1 s2 s3 s1_dum s2_dum t
a = -1;  b = 2;  c = 0.5;  d = 3;  e = -2;  f = 1;
I = @(p,v,L,U) subs(int(p,v),v,U) - subs(int(p,v),v,L);   % definite integral, polynomial limits
R = cell(0,1);
rng(7);
clear('stateNameGenerator');

% States and polynomial values used by the groups S, E and N
x   = pde_var(1,s,[a,b]);    xp  = rpoly({s},1,3);    xd = subs(xp,s,s_dum);
X2  = pde_var(2,s,[a,b]);    X2p = rpoly({s},2,3);    X2d = subs(X2p,s,s_dum);
x0  = pde_var(1);            x0p = randn;
x02 = pde_var(2);            x02p = randn(2,1);
x2  = pde_var(1,[s1;s2],[a,b; c,d]);   x2p = rpoly({s1,s2},1,2);
x1  = pde_var(1,s1,[a,b]);   x1p = rpoly({s1},1,3);
K1 = s^2*s_dum + 3*s_dum^3 - 1;   K2 = 2*s - s^3*s_dum^2;   K0 = 1 + s^2;
K  = s1*s2_dum + s1_dum^2*s2 + 1;  % 2D kernel, in both dummy variables
y22 = pde_var([s1;s2],[0,1;0,1]);   y22p = rpoly({s1,s2},1,2);

if ismember('S',groups)
fprintf('\n===== S: semantics =====\n');
R{end+1} = chk('S1 multiplier K0(s)*x',      K0*x,                          {},{x},{xp},{K0*xp},{{'s'}});
R{end+1} = chk('S2 lower integral',          int(K1*x(s_dum),s_dum,[a,s]),  {},{x},{xp},{I(K1*xd,s_dum,a,s)},{{'s'}});
R{end+1} = chk('S3 upper integral',          int(K2*x(s_dum),s_dum,[s,b]),  {},{x},{xp},{I(K2*xd,s_dum,s,b)},{{'s'}});
R{end+1} = chk('S4 3-PI operator',           K0*x + int(K1*x(s_dum),s_dum,[a,s]) + int(K2*x(s_dum),s_dum,[s,b]), {},{x},{xp}, ...
               {K0*xp + I(K1*xd,s_dum,a,s) + I(K2*xd,s_dum,s,b)},{{'s'}});
R{end+1} = chk('S5 full integral into R',    int((1+s^2)*x,s,[a,b]),        {},{x},{xp},{I((1+s^2)*xp,s,a,b)},{{}});
R{end+1} = chk('S6 int(s*x,s,[a,s,b]): kernel s', int(s*x,s,[a,s,b]),       {},{x},{xp},{I(s*xp,s,a,b)+0*s},{{'s'}});
R{end+1} = chk('S7 s*int(x,s,[a,b]): s outside', s*int(x,s,[a,b]),          s, {x},{xp},{s*I(xp,s,a,b)},{{'s'}});
R{end+1} = chk('S8 full integral, kernel K(s,th)', int(K1*x(s_dum),s_dum,[a,b]), {},{x},{xp},{I(K1*xd,s_dum,a,b)},{{'s'}});
R{end+1} = chk('S9 scaled and negated terms', -3*int(x(s_dum),s_dum,[s,b]) - x, {},{x},{xp},{-3*I(xd,s_dum,s,b) - xp},{{'s'}});
M1 = [s, 1; s_dum^2, 2*s*s_dum];   M0 = [1, s; 0, 2];
R{end+1} = chk('S10 2x2 lower integral + multiplier', M0*X2 + int(M1*X2(s_dum),s_dum,[a,s]), {},{X2},{X2p}, ...
               {M0*X2p + I(M1*X2d,s_dum,a,s)},{{'s'}});
R{end+1} = chk('S11 R -> L2, output space x',  (1+s)*x0,                    x, {x0},{x0p},{(1+s)*x0p},{{'s'}});
R{end+1} = chk('S12 R^2 -> L2, output space x',[s,2]*x02,                   x, {x02},{x02p},{[s,2]*x02p},{{'s'}});
R{end+1} = chk('S13 (x0,x) -> L2',             x0 + int(x(s_dum),s_dum,[a,s]), {},{x0,x},{x0p,xp},{x0p + I(xd,s_dum,a,s)},{{'s'}});
R{end+1} = chk('S14 (x0,x) -> R',              2*x0 + int(s*x,s,[a,b]),     {},{x0,x},{x0p,xp},{2*x0p + I(s*xp,s,a,b)},{{}});
R{end+1} = chk('S15 rows R x L2',              [x0 + int(x,s,[a,b]); s*x0 + K0*x], {},{x0,x},{x0p,xp}, ...
               {x0p + I(xp,s,a,b); s*x0p + K0*xp},{{},{'s'}});
R{end+1} = chk('S16 rows, one output space per row', [int(x,s,[a,b]); 3*x0], {{};s},{x0,x},{x0p,xp},{I(xp,s,a,b); 3*x0p+0*s},{{},{'s'}});
R{end+1} = chk('S17 zero term keeps a column', x0 + 0*x,                    {},{x0,x},{x0p,xp},{x0p+0*s},{{'s'}});
w = pde_var('in',1);   u = pde_var('control',1,s,[a,b]);   wp = randn;   up = rpoly({s},1,2);
R{end+1} = chk('S18 state, input w and control u', K0*x + s*w + int(K1*u(s_dum),s_dum,[a,s]), {},{x,w,u},{xp,wp,up}, ...
               {K0*xp + s*wp + I(K1*subs(up,s,s_dum),s_dum,a,s)},{{'s'}});
Ka = subs(K,s2_dum,s2);   Kb = subs(K,s1_dum,s1);
x2d1 = subs(x2p,s1,s1_dum);   x2d2 = subs(x2p,s2,s2_dum);   x2d12 = subs(x2d1,s2,s2_dum);
R{end+1} = chk('S19 2D multiplier',            (s1+s2^2)*x2,                {},{x2},{x2p},{(s1+s2^2)*x2p},{{'s1','s2'}});
R{end+1} = chk('S20 2D lower in s1, multiplier in s2', int(Ka*x2(s1_dum,s2),s1_dum,[a,s1]), {},{x2},{x2p},{I(Ka*x2d1,s1_dum,a,s1)},{{'s1','s2'}});
R{end+1} = chk('S21 2D upper in s2, multiplier in s1', int(Kb*x2(s1,s2_dum),s2_dum,[s2,d]), {},{x2},{x2p},{I(Kb*x2d2,s2_dum,s2,d)},{{'s1','s2'}});
R{end+1} = chk('S22 2D lower in s1, upper in s2', int(K*x2(s1_dum,s2_dum),[s1_dum;s2_dum],[a,s1; s2,d]), {},{x2},{x2p}, ...
               {I(I(K*x2d12,s1_dum,a,s1),s2_dum,s2,d)},{{'s1','s2'}});
R{end+1} = chk('S23 2D full in s2 -> L2[s1]',  int(s2*x2,s2,[c,d]),         {},{x2},{x2p},{I(s2*x2p,s2,c,d)},{{'s1'}});
R{end+1} = chk('S24 2D full in both -> R',     int(s1*s2*x2,[s1;s2],[a,b; c,d]), {},{x2},{x2p},{I(I(s1*s2*x2p,s1,a,b),s2,c,d)},{{}});
R{end+1} = chk('S25 2D full in s2, upper in s1', int(Ka*x2(s1_dum,s2),[s1_dum;s2],[s1,b; c,d]), {},{x2},{x2p}, ...
               {I(I(Ka*x2d1,s1_dum,s1,b),s2,c,d)},{{'s1'}});
z12 = pde_var('out',1,[s1;s2],[a,b; c,d]);
R{end+1} = chk('S26 L2[s1] -> L2[s1,s2]',      (s1+s2)*x1,                  z12,{x1},{x1p},{(s1+s2)*x1p},{{'s1','s2'}});
R{end+1} = chk('S27 (x1,x2) -> L2[s1,s2]',     s2*x1 + x2,                  {},{x1,x2},{x1p,x2p},{s2*x1p + x2p},{{'s1','s2'}});
R{end+1} = chk('S28 (x1,x2) -> L2[s1]',        x1 + int(x2,s2,[c,d]),       {},{x1,x2},{x1p,x2p},{x1p + I(x2p,s2,c,d)},{{'s1'}});
x3 = pde_var(1,[s1;s2;s3],[a,b; c,d; e,f]);   x3p = rpoly({s1,s2,s3},1,2);   K3 = s1*s1_dum + s3^2;
R{end+1} = chk('S29 3D lower s1, full s2, multiplier s3', int(K3*x3(s1_dum,s2,s3),[s1_dum;s2],[a,s1; c,d]), {},{x3},{x3p}, ...
               {I(I(K3*subs(x3p,s1,s1_dum),s1_dum,a,s1),s2,c,d)},{{'s1','s3'}});
R{end+1} = chk('S30 full in s2 into L2[s1,s2]: constant in s2', int(y22,s2,[0,1]), [s1,s2],{y22},{y22p},{I(y22p,s2,0,1)+0*s2},{{'s1','s2'}});
R{end+1} = chk('S31 0*x is the zero operator', 0*x,                         {},{x},{xp},{0*s},{{'s'}});
R{end+1} = chk('S32 x - x is the zero operator', x - x,                     {},{x},{xp},{0*s},{{'s'}});
% The forms of the output space give one operator.
F = int(x2,s2,[c,d]);
P1 = get_op(F,[s1,s2]);  P2 = get_op(F,{'s1','s2'});  P3 = get_op(F,{[s1,s2]});  P4 = get_op(F,z12);
mx = max([opdiff(P1,P2),opdiff(P1,P3),opdiff(P1,P4)]);
R{end+1} = report('S33 output space as polynomial, cellstr, cell, pde_struct: equal',mx<1e-12,sprintf('max |diff| %.1e',mx));
end

if ismember('E',groups)
fprintf('\n===== E: unsupported terms =====\n');
y = pde_var(s,[0,1]);   xr = pde_var();
R{end+1} = errm('E1 derivative',                  @() get_op(diff(x,s)),                      'Derivatives');
R{end+1} = errm('E2 evaluation at a boundary',    @() get_op(subs(x,s,a)),                    'Evaluation at a position');
R{end+1} = errm('E3 product of states',           @() get_op(x*x),                            'not supported for nonlinear terms');
R{end+1} = errm('E4 delay',                       @() get_op(subs(x,t,t-1)),                  'Delayed');
R{end+1} = errm('E5 output variable',             @() get_op(pde_var('out',1)),               'Outputs');
R{end+1} = errm('E6 time derivative',             @() get_op(diff(y,t)),                      'temporal derivatives');
R{end+1} = errm('E7 equation',                    @() get_op(y==0),                           'free terms');
R{end+1} = errm('E8 partial integral, s not in the output space', @() get_op(int(x(s_dum),s_dum,[a,s]),{}), 'requires ''s''');
R{end+1} = errm('E9 state on s into a space without s', @() get_op(x,{}),                     'depends on variable ''s''');
R{end+1} = errm('E10 2D state into L2[s1], no integral over s2', @() get_op(y22,s1),          'depends on variable ''s2''');
R{end+1} = errm('E11 x in R, s*x into R',         @() get_op(s*xr,xr),                        'depend on ''s''');
R{end+1} = errm('E12 variable r of no state, output space y', @() get_op(pvar('r')*y,y),      'depend on ''r''');
R{end+1} = errm('E13 dummy variable in a multiplier', @() get_op(s_dum*y,s),                  'depend on ''s_dum''');
R{end+1} = errm('E14 dummy variable left after int over s', @() get_op(int(s_dum*y,s,[0,s])), 'outside of an integral');
R{end+1} = errm('E15 dummy variable left in a 2D kernel', @() get_op(int(K*x2(s1_dum,s2),s1_dum,[a,s1])), 'outside of an integral');
R{end+1} = errm('E16 y(s_dum) without integral, output s', @() get_op(y(s_dum),s),            'without an integral');
R{end+1} = errm('E17 y(s_dum) without integral',  @() get_op(y(s_dum)),                       'outside of an integral');
R{end+1} = errm('E18 time in a coefficient',      @() get_op(t*y),                            'is the time');
R{end+1} = errm('E19 time in the output space',   @() get_op(x,[s,t]),                        'is the time');
R{end+1} = errm('E20 integral over a subinterval',@() get_op(int(y(s_dum),s_dum,[0,0.5])),    'limit of the domain');
R{end+1} = errm('E21 lower limit inside the domain', @() get_op(int(y(s_dum),s_dum,[0.2,s])), 'limit of the domain');
R{end+1} = errm('E22 limit 2*s',                  @() get_op(int(y(s_dum),s_dum,[0,2*s])),    'single variables');
R{end+1} = errm('E23 one variable on two domains',@() get_op(x + pde_var(1,s,[0,1])),         'different domains');
R{end+1} = errm('E24 output space: wrong number of lists', @() get_op([x;x],{s}),             'one list of variables per free equation');
R{end+1} = errm('E25 output space: numeric',      @() get_op(y,3),                            'list of variables');
R{end+1} = errm('E26 output space: pde_struct of other size', @() get_op(y,pde_var(2,s,[0,1])), 'size of PDE variable');
R{end+1} = errm('E27 output space: pde_struct on other domain', @() get_op(y,pde_var(1,s,[0,2])), 'different domains');
% Nonlinear terms, alone or in a sum with linear terms
R{end+1} = errm('E28 nonlinear: x^2',             @() get_op(x^2),                            'not supported for nonlinear terms');
R{end+1} = errm('E29 nonlinear: R state times L2 state', @() get_op(x0*x),                    'not supported for nonlinear terms');
R{end+1} = errm('E30 nonlinear term in a linear sum', @() get_op(K0*x + s*x0*x),              'not supported for nonlinear terms');
R{end+1} = errm('E31 nonlinear: 2D product',      @() get_op(y22*y22),                        'not supported for nonlinear terms');
% Nonlinear terms that the parser refuses to build
wn = pde_var('in',1);
R{end+1} = errm('E32 state times input (refused when built)', @() get_op(x*wn),               'not supported');
R{end+1} = errm('E33 integral of x*x (refused when built)', @() get_op(int(x*x,s,[a,b])),     'not supported');
R{end+1} = errm('E34 x times x_t (refused when built)', @() get_op(x*diff(x,t)),              'not supported');
end

if ismember('N',groups)
fprintf('\n===== N: new spatial variables on [0,1] =====\n');
r = pvar('r');   xn = pde_var();   xnp = randn;   y = pde_var(s,[0,1]);   yp = rpoly({s},1,2);
Q1 = get_op(s*xn);
R{end+1} = report('N1 get_op(s*x), x in R: L2[s] on [0,1]',isequal(Q1.vars,{'s'}) && isequal(Q1.dom,[0,1]),'');
R{end+1} = chk('N2 get_op(s*x): (Q1 x)(s) = s*x', s*xn,                     {},{xn},{xnp},{s*xnp},{{'s'}});
QQ = Q1'*Q1;   mx = pmax(apply_sopvar(QQ.C{1,1},1) - 1/3);
R{end+1} = report('N3 Q1''*Q1 = int_0^1 s^2 ds = 1/3',mx<1e-12,sprintf('|diff| %.1e',mx));
R{end+1} = chk('N4 R^2 -> L2, output space by list', [s,2]*x02,             s, {x02},{x02p},{[s,2]*x02p},{{'s'}});
P = get_op([s,2]*x02,s);
R{end+1} = report('N5 ... on [0,1]',isequal(P.dom,[0,1]),'');
R{end+1} = chk('N6 R -> L2[s3] by list',       xn,                          s3,{xn},{xnp},{xnp+0*s3},{{'s3'}});
R{end+1} = chk('N7 s1*s2*x: new pair L2[s1,s2]', s1*s2*xn,                  {},{xn},{xnp},{s1*s2*xnp},{{'s1','s2'}});
P = get_op(s1*s2*xn);
R{end+1} = report('N8 ... on [0,1]^2',isequal(P.dom,[0,1;0,1]),'');
R{end+1} = chk('N9 r*y: new r next to the variable s of y', r*y,            {},{y},{yp},{r*yp},{{'r','s'}});
R{end+1} = chk('N10 r*int_0^s y: new r, partial integral in s', r*int(y(s_dum),s_dum,[0,s]), {},{y},{yp}, ...
               {r*I(subs(yp,s,s_dum),s_dum,0,s)},{{'r','s'}});
y02 = pde_var(1,s,[0,2]);   y02p = rpoly({s},1,2);
R{end+1} = chk('N11 [s*x; y]: s of y on [0,2]',  [s*xn; y02],               {},{xn,y02},{xnp,y02p},{s*xnp; y02p},{{'s'},{'s'}});
P = get_op([s*xn; y02]);
R{end+1} = report('N12 ... s keeps the domain [0,2] of y',isequal(P.dom,[0,2]),'');
end

if ismember('R',groups)
fprintf('\n===== R: random operators (%d trials) =====\n',ntrials);
for trial = 1:ntrials
    R{end+1} = random_trial(trial); %#ok<AGROW>
end
end

if ismember('C',groups)
nc = max(1,round(ntrials/4));
fprintf('\n===== C: cross-check with the converter (%d trials) =====\n',nc);
for trial = 1:nc
    R{end+1} = convert_trial(trial); %#ok<AGROW>
end
end

results = [R{:}];
fprintf('\n%d of %d cases passed.\n',sum([results.passed]),numel(results));
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function r = chk(name,terms,sout,vars_in,polys_in,expected,exp_vars)
% R = CHK(NAME,TERMS,SOUT,VARS_IN,POLYS_IN,EXPECTED,EXP_VARS) extracts the
% operators of TERMS, applies them to the polynomial values of the states
% and inputs, and compares the result and the output spaces with the
% expected ones.
% INPUT
% - name:     'char', the name of the case;
% - terms:    'pde_struct', free terms;
% - sout:     the output space passed to get_op, or {} for the default;
% - vars_in:  cell of 'pde_struct', the PDE variables of the terms;
% - polys_in: cell of 'polynomial' column vectors, their values;
% - expected: cell of 'polynomial' column vectors, the value of each row;
% - exp_vars: cell of cellstr, the variables of the output space of each row;
% OUTPUT
% - r:        'struct', the result.
try
    if iscell(sout) && isempty(sout),  [Px,Pw,Pu] = get_op(terms);
    else,                              [Px,Pw,Pu] = get_op(terms,sout);
    end
    M = numel(terms.free);   out = cell(M,1);   sp_ok = true;
    for i = 1:M,  out{i} = polynomial(zeros(terms.free{i}.size,1));  end
    Ps = {Px,Pw,Pu};   objs = {'x','w','u'};
    for kk = 1:3
        P = Ps{kk};   obj = objs{kk};
        if isempty(P),  continue;  end
        for j = 1:numel(terms.(obj))
            k = find(cellfun(@(v) ~isempty(v.(obj)) && v.(obj){1}.ID==terms.(obj){j}.ID,vars_in),1);
            for i = 1:M
                if ~isempty(P.C{i,j}),  out{i} = out{i} + apply_sopvar(P.C{i,j},polys_in{k});  end
            end
        end
        for i = 1:M
            sp_ok = sp_ok && isempty(setxor(P.vars(P.space_out(i,:)),exp_vars{i}));
        end
    end
    mx = 0;
    for i = 1:M,  mx = max(mx,pmax(out{i} - expected{i}));  end
    r = report(name,sp_ok && mx<1e-9,sprintf('max |diff| %.1e%s',mx,repmat(', OUTPUT SPACE WRONG',1,~sp_ok)));
catch ME
    r = report(name,false,sprintf('%s [%s:%d]',ME.message,ME.stack(1).name,ME.stack(1).line));
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function r = errm(name,f,phrase)
% R = ERRM(NAME,F,PHRASE) passes if F throws an error whose message
% contains PHRASE.
% INPUT
% - name:   'char', the name of the case;
% - f:      function handle with no input;
% - phrase: 'char', the text the message must contain;
% OUTPUT
% - r:      'struct', the result.
try
    f();
    r = report(name,false,'no error');
catch ME
    r = report(name,contains(ME.message,phrase),ME.message);
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function r = random_trial(trial)
% R = RANDOM_TRIAL(TRIAL) declares random states and output rows, builds
% random terms with their values computed by polynomial integration, and
% checks the output spaces, the values and the adjoint of the extracted
% operator. The trial is reproducible from its number.
% INPUT
% - trial: positive integer, the number of the trial (the seed);
% OUTPUT
% - r:     'struct', the result.
rng(1000+trial);
V = {'s1','s2','s3'};   pv = {pvar('s1'),pvar('s2'),pvar('s3')};
dom = [-1,2; 0,3; 0.5,1.5];
nd = randi(3);   use_new = nd<3 && rand<0.35;   % s3 as a variable of no state
% States: R or L2 over a random subset of s1..s_nd, in random order
nx = randi(3);   X = cell(nx,1);   Xp = cell(nx,1);   Xv = cell(nx,1);
for j = 1:nx
    sub = find(rand(1,nd)<0.6);   sub = sub(randperm(numel(sub)));   dj = randi(2);
    if isempty(sub),  X{j} = pde_var(dj);
    else,             X{j} = pde_var(dj,vertcat(pv{sub}),dom(sub,:));
    end
    Xv{j} = sub;   Xp{j} = rpoly(pv(sub),dj,2);
end
% Output rows, and their terms with expected values
M = randi(3);   O = cell(M,1);   rows = cell(M,1);   expect = cell(M,1);
for i = 1:M
    O{i} = find(rand(1,nd)<0.5);
    if use_new && rand<0.6,  O{i} = [O{i},3];  end
    dOi = randi(2);   expect{i} = polynomial(zeros(dOi,1));
    cols = find(rand(1,nx)<0.7);   if isempty(cols),  cols = randi(nx);  end
    for j = cols
        [trm,ex] = rand_block(X{j},Xp{j},Xv{j},O{i},dOi,dom);
        if isempty(rows{i}),  rows{i} = trm;  else,  rows{i} = rows{i} + trm;  end
        expect{i} = expect{i} + ex;
    end
end
terms = vertcat(rows{:});
name = sprintf('R%d: %d variables, %d states, %d rows%s',trial,nd,nx,M,repmat(', new s3',1,use_new));
try
    P = get_op(terms,cellfun(@(o) V(o),O,'UniformOutput',false));
    % Output spaces
    sp_ok = true;
    for i = 1:M,  sp_ok = sp_ok && isempty(setxor(P.vars(P.space_out(i,:)),V(O{i})));  end
    % Values of P*x, column by column
    idx = cellfun(@(c) find(cellfun(@(v) v.x{1}.ID==c.ID,X),1),terms.x);
    out = cell(M,1);   mx = 0;
    for i = 1:M
        out{i} = polynomial(zeros(P.dim_out(i),1));
        for jj = 1:numel(idx)
            if ~isempty(P.C{i,jj}),  out{i} = out{i} + apply_sopvar(P.C{i,jj},Xp{idx(jj)});  end
        end
        mx = max(mx,pmax(out{i} - expect{i}));
    end
    % Adjoint: <P x, z> = <x, P' z>, on the domains of P
    Pt = P';   lhs = 0;   rhs = 0;   Z = cell(M,1);
    for i = 1:M
        Z{i} = rpoly(cellfun(@pvar,V(O{i}),'UniformOutput',false),P.dim_out(i),2);
        lhs = lhs + ip(out{i},Z{i},P,P.space_out(i,:));
    end
    for jj = 1:numel(idx)
        yj = polynomial(zeros(P.dim_in(jj),1));
        for i = 1:M
            if ~isempty(Pt.C{jj,i}),  yj = yj + apply_sopvar(Pt.C{jj,i},Z{i});  end
        end
        rhs = rhs + ip(Xp{idx(jj)},yj,P,P.space_in(jj,:));
    end
    adj = abs(lhs-rhs)/max(1,abs(lhs));
    r = report(name,sp_ok && mx<1e-8 && adj<1e-8,sprintf('values %.1e, adjoint %.1e%s',mx,adj,repmat(', OUTPUT SPACE WRONG',1,~sp_ok)));
catch ME
    r = report(name,false,sprintf('%s [%s:%d]',ME.message,ME.stack(1).name,ME.stack(1).line));
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function r = convert_trial(trial)
% R = CONVERT_TRIAL(TRIAL) declares random states in R or L2 over s1, s2,
% and the PDE x_t = F(x) with F random terms into the spaces of the
% states; the PIE then has T = I and A the operator of F, which the
% converter builds by its own code. The trial is reproducible from its
% number.
% INPUT
% - trial: positive integer, the number of the trial (the seed);
% OUTPUT
% - r:     'struct', the result.
rng(2000+trial);
clear('stateNameGenerator');
pv = {pvar('s1'),pvar('s2')};   dom = [0,1; -1,1];   t = pvar('t');
nd = randi(2);   nx = randi(3);   X = cell(nx,1);   Xv = cell(nx,1);
for j = 1:nx
    sub = find(rand(1,nd)<0.7);   dj = randi(2);
    if isempty(sub),  X{j} = pde_var(dj);  else,  X{j} = pde_var(dj,vertcat(pv{sub}),dom(sub,:));  end
    Xv{j} = sub;
end
F = cell(nx,1);   eqs = cell(nx,1);
for i = 1:nx
    for j = 1:nx
        if rand<0.4 && j~=i,  continue;  end
        trm = rand_block(X{j},[],Xv{j},Xv{i},X{i}.x{1}.size,dom);
        if isempty(F{i}),  F{i} = trm;  else,  F{i} = F{i} + trm;  end
    end
    eqs{i} = diff(X{i},t)==F{i};
end
name = sprintf('C%d: %d variables, %d states',trial,nd,nx);
try
    Pg = get_op(vertcat(F{:}),vertcat(X{:}));
    Ig = get_op(vertcat(X{:}));
    [~,PIE] = evalc('convert(initialize(vertcat(eqs{:}),true),''pie'')');
    % The converter orders the components by their spaces: match each PIE
    % component to the first unmatched one with the same space and size.
    perm = zeros(1,nx);   used = false(1,nx);
    for k = 1:nx
        j = find(~used & all(Pg.space_in==PIE.A.space_in(k,:),2)' & Pg.dim_in(:)'==PIE.A.dim_in(k),1);
        perm(k) = j;   used(j) = true;
    end
    dA = opdiff(PIE.A,permop(Pg,perm));   dT = opdiff(PIE.T,permop(Ig,perm));
    r = report(name,dA<1e-9 && dT<1e-9,sprintf('|A - get_op(F)| %.1e, |T - I| %.1e',dA,dT));
catch ME
    r = report(name,false,sprintf('%s [%s:%d]',ME.message,ME.stack(1).name,ME.stack(1).line));
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [trm,ex] = rand_block(Xj,Xpj,Xvj,Oi,dOi,dom)
% [TRM,EX] = RAND_BLOCK(XJ,XPJ,XVJ,OI,DOI,DOM) returns a random PI term of
% the state XJ into the output space with variables OI, and its value EX
% at the polynomial XPJ computed by polynomial integration: along each
% shared variable a random multiplier, lower, upper or full integral,
% along each variable of the state only a full integral, and a random
% kernel in the variables of the output and of integration.
% INPUT
% - Xj:   'pde_struct', a PDE state variable;
% - Xpj:  'polynomial' column vector, its value ([] for no value);
% - Xvj:  integer array, the indices in s1,s2,s3 of its variables, in order;
% - Oi:   integer array, the indices of the variables of the output space;
% - dOi:  positive integer, the size of the output;
% - dom:  3 x 2 array, the domains of s1, s2, s3;
% OUTPUT
% - trm:  'pde_struct', the term;
% - ex:   'polynomial' dOi x 1, its value ([] if Xpj is empty).
pv = {pvar('s1'),pvar('s2'),pvar('s3')};   pd = {pvar('s1_dum'),pvar('s2_dum'),pvar('s3_dum')};
typ = zeros(1,numel(Xvj));                 % 1 multiplier, 2 lower, 3 upper, 4 full, 5 full along own variable
args = cell(1,numel(Xvj));   kvars = pv(Oi);   xs = Xpj;
for k = 1:numel(Xvj)
    v = Xvj(k);
    if ismember(v,Oi),  typ(k) = randi(4);  else,  typ(k) = 5;  end
    if typ(k)==1 || typ(k)==5
        args{k} = pv{v};
    else
        args{k} = pd{v};   kvars{end+1} = pd{v}; %#ok<AGROW>
        if ~isempty(xs),  xs = subs(xs,pv{v},pd{v});  end
    end
    if typ(k)==5,  kvars{end+1} = pv{v};  end %#ok<AGROW>
end
K = rpoly(kvars,dOi,2,Xj.x{1}.size);
if isempty(Xvj),  xe = Xj;  else,  xe = Xj(args{:});  end
trm = K*xe;   ex = [];
if ~isempty(xs),  ex = K*xs;  end
for k = 1:numel(Xvj)
    v = Xvj(k);   L = dom(v,1);   U = dom(v,2);
    switch typ(k)
        case 2,  trm = int(trm,pd{v},[L,pv{v}]);   lim = {pd{v},L,pv{v}};
        case 3,  trm = int(trm,pd{v},[pv{v},U]);   lim = {pd{v},pv{v},U};
        case 4,  trm = int(trm,pd{v},[L,U]);       lim = {pd{v},L,U};
        case 5,  trm = int(trm,pv{v},[L,U]);       lim = {pv{v},L,U};
        otherwise,  continue
    end
    if ~isempty(ex)
        G = int(ex,lim{1});   ex = subs(G,lim{1},lim{3}) - subs(G,lim{1},lim{2});
    end
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function val = ip(p,q,P,mask)
% VAL = IP(P,Q,P,MASK) returns the L2 inner product of the polynomial
% vectors P and Q over the variables MASK of the registry of the operator,
% on its domains.
% INPUT
% - p, q:  'polynomial' column vectors of equal size;
% - P:     'opvar' object, whose registry gives the variables and domains;
% - mask:  logical row, the variables of the registry to integrate over;
% OUTPUT
% - val:   'double' scalar.
g = sum(polynomial(p).*polynomial(q));
for k = find(mask)
    v = pvar(P.vars{k});   G = int(g,v);
    g = subs(G,v,P.dom(k,2)) - subs(G,v,P.dom(k,1));
end
val = double(g);
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function p = rpoly(vars,m,deg,n)
% P = RPOLY(VARS,M,DEG,N) returns an M x N polynomial matrix of total
% degree at most DEG in the variables VARS, with random coefficients.
% INPUT
% - vars:  cell of 'polynomial' pvars (may be empty);
% - m, n:  positive integers, the size (n defaults to 1);
% - deg:   nonnegative integer, the degree;
% OUTPUT
% - p:     m x n 'polynomial'.
if nargin<4,  n = 1;  end
p = polynomial(randn(m,n));   mons = {polynomial(1)};
for dd = 1:deg
    new = {};
    for k = 1:numel(vars)
        for q = 1:numel(mons),  new{end+1} = mons{q}*vars{k};  end %#ok<AGROW>
    end
    mons = [mons,new]; %#ok<AGROW>
end
for q = 2:numel(mons)
    if rand<0.5,  p = p + randn(m,n)*mons{q};  end
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function mx = opdiff(A,B)
% MX = OPDIFF(A,B) returns the largest coefficient of A - B, or Inf if the
% spaces, domains or dimensions of A and B differ.
% INPUT
% - A, B:  'opvar' objects;
% OUTPUT
% - mx:    'double' scalar.
if ~isequal(A.vars,B.vars) || ~isequal(A.dom,B.dom) || ~isequal(A.space_out,B.space_out) ...
        || ~isequal(A.space_in,B.space_in) || ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    mx = inf;   return
end
D = A-B;   mx = 0;
for k = 1:numel(D.C)
    if isempty(D.C{k}),  continue;  end
    for j = 1:numel(D.C{k}.params)
        pj = D.C{k}.params{j};
        if ~isempty(pj),  mx = max(mx,full(max(abs(pj(:)))));  end
    end
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Q = permop(P,perm)
% Q = PERMOP(P,PERM) returns P with its block rows and columns in the
% order PERM.
% INPUT
% - P:     square 'opvar' object (as a block grid);
% - perm:  permutation of 1:size(P.C,1);
% OUTPUT
% - Q:     'opvar' object.
meta = struct('vars',{P.vars},'dom',P.dom,'space_out',P.space_out(perm,:), ...
              'space_in',P.space_in(perm,:),'dim_out',P.dim_out(perm),'dim_in',P.dim_in(perm));
Q = opvar(P.C(perm,perm),meta);
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function mx = pmax(p)
% MX = PMAX(P) returns the largest absolute coefficient of the polynomial P.
% INPUT
% - p:   'polynomial' or 'double';
% OUTPUT
% - mx:  'double' scalar.
p = polynomial(p);
if isempty(p.coefficient),  mx = 0;  else,  mx = full(max(abs(p.coefficient(:))));  end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function r = report(name,passed,msg)
% R = REPORT(NAME,PASSED,MSG) prints and returns the result of a case.
% INPUT
% - name:    'char', the name of the case;
% - passed:  logical scalar;
% - msg:     'char', details;
% OUTPUT
% - r:       'struct' with fields name, passed and message.
if passed,  tag = 'PASS';  else,  tag = 'FAIL';  end
fprintf('%s  %-52s %s\n',tag,name,msg);
r = struct('name',name,'passed',logical(passed),'message',msg);
end
