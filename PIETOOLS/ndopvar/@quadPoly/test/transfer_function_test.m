clear;

pvar t s
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);   
z = pde_var('output',1);
w = pde_var('input',1);          u = pde_var('control',1);
% Declare the sytem equations
PDE = [diff(x,t) == diff(x,s) + s*w;
       z == int(x,s,[0,1]);
       subs(x,s,0)==0;];

display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;      A = PIE.A;
B1 = PIE.B1;    C1 = PIE.C1;
Astar = A'; C1star = C1';
Tstar = T';
B1star = B1';
%%
syms s1 t1;
lamVals = linspace(2,10,10);
Minvt = cell(1,length(lamVals));
RM1t = Minvt; RM2t = Minvt; lamT_A_inv_lamT_A = RM1t;
for i=1:numel(lamVals)
    lam = lamVals(i);
    M = quadPoly.Sym2quadPoly(-1+0*s1,{'s1'},{});
    FG1 = quadPoly.Sym2quadPoly(lam*(-1)+0*s1,{'s1'},{'t1'});
    FG2 = quadPoly.Sym2quadPoly(lam*(0)+0*s1,{'s1'},{'t1'});
    [Minv,RM1,RM2,~] = sopvar.inv_1D(M,FG1,FG2,[0,1]);
    Minvt{i} = Minv;
    RM1t{i} = RM1;
    RM2t{i} = RM2;
    
    lamT_A = lam*T-A;
    inv_lamT_A = T;
    inv_lamT_A.R = struct('R0',quadPoly2polynomial(Minv),...
                            'R1',quadPoly2polynomial(RM1),...
                             'R2',quadPoly2polynomial(RM2));

    lamT_A_inv_lamT_A{i} = (lam*T-A)*inv_lamT_A;
end