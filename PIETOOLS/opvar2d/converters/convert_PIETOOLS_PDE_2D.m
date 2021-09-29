%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_2D.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE = convert_PIETOOLS_PDE_2D(PDE)
% Convert 2D PDE to associated PIE
% PDE must be a data structure with fields:
% - dom: 2x2 array of domain [a1,b1;a2,b2] of spatial variables.
% - n: Structure with field:
%       > n_pde: 3x3 DIAGONAL matrix with elements n_pde(i,j) specifying
%                how many state elements occur which are differentiable up 
%                to order (i-1) in x and (j-1) in y;
%       > nx: Integer providing number of ODE state components.
% - ODE: Structure specifying ODE X_{t} = A*X + Bxr*r with fields:
%       > A: Matrix A in ODE system;
%       > Bxr: Matrix Bxr in ODE system, mapping PDE output r to ODE;
%       > Cv: Matrix Cv in ODE system, mapping X to PDE input v=Cv*X.
% - PDE: Structure specifying PDE defined as
% \dot{[u_0] = A{1,1} [u_0] + A{2,1} [u_1x] + A{1,2} [u_1y] + A{2,2} [u_1xy] + A{3,1} [u_2xx] + A{1,3} [u_2yy] + A{3,2} [u_2xxy] + A{2,3} [u_2xyy] + A{3,3} [u_2xxyy]
%      [u_1]          [u_1]          [u_2x]          [u_2y]          [u_2xy]
%      [u_2]}         [u_2]
%        Must have fields A, Bpv and Drb, where Bpv describes ODE->PDE
%        interconnection, \dot{u} = ... + Bpv*v, and Drb describes PDE->ODE
%        interconnection, r = Drb*u_bf.
% - BC: Structure with field Ebb, which must be an opvar2d object
%       specifying BCs as 0 = Ebb*ubf
%
% Here:
% u_bf = [u_1(a1,a2); u_1(b1,a2); u_1(a1,b2); u_1(b1,b2); 
%         u_2(a1,a2); u_2(b1,a2); u_2(a1,b2); u_2(b1,b2); u_2x(a1,a2); u_2x(b1,a2); u_2x(a1,b2); u_2x(b1,b2);
%         u_2y(a1,a2); u_2y(b1,a2); u_2y(a1,b2); u_2y(b1,b2); u_2xy(a1,a2); u_2xy(b1,a2); u_2xy(a1,b2); u_2xy(b1,b2);
%         u_1x(x,a2); u_1x(x,b2); 
%         u_2xx(x,a2); u_2xx(x,b2); u_2xxy(x,a2); u_2xxy(x,b2);
%         u_1y(a1,y); u_1y(b1,y);
%         u_2yy(a1,y); u_2yy(b1,y); u_2xyy(a1,y); u_2xyy(b1,y)]
%
%
%
% This system gets transformed into PIE:
%
% T*\dot{[X]        = A*[X]
%        [\hat{u}]}     [\hat{u}]
%
% where \hat{u} = [u_0; u_1xy; u_2xxyy] is the fundamental state, and
% where T and A are opvar2d, collectied in structure PIE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Initialization

% Initialize PDE in case this has not been done
PDE = initialize_PIETOOLS_PDE_2D(PDE);

% Extract dimensions of PDE state components
nx = PDE.n.nx;
n00 = PDE.n.n_pde(1,1);
n11 = PDE.n.n_pde(2,2);
n22 = PDE.n.n_pde(3,3);

% Extract lower and  upper boundaries of spatial domain
dom = PDE.dom;
a = dom(:,1);   b = dom(:,2);

% Initialize PIE variables, as well as dummy variables eta
var1 = PDE.vars(:,1);   var2 = PDE.vars(:,2);
ss1 = var1(1);          tt1 = var2(1);
ss2 = var1(2);          tt2 = var2(2);
pvar eta1 eta2;

% Extract PDE defining matrix Aij
Am = PDE.PDE.A;
for m=1:numel(Am)
    Am{m} = cell2mat(Am{m}.coeff);
end

% Extract BCs defining operator Ebb
Ebb = PDE.BC.Ebb;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Describe full bc and full state in terms of core bc and fundamental state:
%
% [u_bf] = [K11] [u_bc] + [K12] [\hat{u}_11]
%                               [\hat{u}_22]
%
% [u_11] = [K21] [u_bc] + [K22] [\hat{u}_11]
% [u_22]                        [\hat{u}_22]

% Map from core bc to full bc
K11 = opvar2d();
K11.I = [a,b];      K11.var1 = var1;  K11.var2 = var2;

R00_1 = [eye(n11);eye(n11);eye(n11);eye(n11)];
R00_2 = kron([1,0,0,0;1,b(1)-a(1),0,0;1,0,b(2)-a(2),0;1,b(1)-a(1),b(2)-a(2),(b(2)-a(2))*(b(1)-a(1));
              0,1,0,0;0,1,0,0;0,1,0,b(2)-a(2);0,1,0,b(2)-a(2);
              0,0,1,0;0,0,1,b(1)-a(1);0,0,1,0;0,0,1,b(1)-a(1);
              0,0,0,1;0,0,0,1;0,0,0,1;0,0,0,1],eye(n22));
          
K11.R00 = [R00_1 , zeros(size(R00_1,1),size(R00_2,2));
            zeros(size(R00_2,1),size(R00_1,2)) , R00_2];

        
R0x_1 = [zeros(n11);eye(n11);zeros(n11);eye(n11)];        
R0x_2 = polynomial(zeros(16*n22,2*n22));
R0x_2(n22+1:2*n22,1:n22) = (b(1)-ss1)*eye(n22);
R0x_2(3*n22+1:4*n22,1:n22) = (b(1)-ss1)*eye(n22);
R0x_2(3*n22+1:4*n22,n22+1:2*n22) = (b(2)-a(2))*(b(1)-ss1)*eye(n22);
R0x_2(5*n22+1:6*n22,1:n22) = eye(n22);
R0x_2(7*n22+1:8*n22,1:n22) = eye(n22);
R0x_2(7*n22+1:8*n22,n22+1:2*n22) = (b(2)-a(2))*eye(n22);
R0x_2(9*n22+1:10*n22,n22+1:2*n22) = (b(1)-ss1)*eye(n22);
R0x_2(11*n22+1:12*n22,n22+1:2*n22) = (b(1)-ss1)*eye(n22);
R0x_2(13*n22+1:14*n22,n22+1:2*n22) = eye(n22);
R0x_2(15*n22+1:16*n22,n22+1:2*n22) = eye(n22);

K11.R0x = [R0x_1 , zeros(size(R0x_1,1),size(R0x_2,2));
            zeros(size(R0x_2,1),size(R0x_1,2)) , R0x_2];
    
        
R0y_1 = [zeros(n11);zeros(n11);eye(n11);eye(n11)];
R0y_2 = polynomial(zeros(16*n22,2*n22));
R0y_2(2*n22+1:3*n22,1:n22) = (b(2)-ss2)*eye(n22);
R0y_2(3*n22+1:4*n22,1:n22) = (b(2)-ss2)*eye(n22);
R0y_2(3*n22+1:4*n22,n22+1:2*n22) = (b(1)-a(1))*(b(2)-ss2)*eye(n22);
R0y_2(6*n22+1:7*n22,n22+1:2*n22) = (b(2)-ss2)*eye(n22);
R0y_2(7*n22+1:8*n22,n22+1:2*n22) = (b(2)-ss2)*eye(n22);
R0y_2(10*n22+1:11*n22,1:n22) = eye(n22);
R0y_2(11*n22+1:12*n22,1:n22) = eye(n22);
R0y_2(11*n22+1:12*n22,n22+1:2*n22) = (b(1)-a(1))*eye(n22);
R0y_2(14*n22+1:15*n22,n22+1:2*n22) = eye(n22);
R0y_2(15*n22+1:16*n22,n22+1:2*n22) = eye(n22);

K11.R0y = [R0y_1 , zeros(size(R0y_1,1),size(R0y_2,2));
            zeros(size(R0y_2,1),size(R0y_1,2)) , R0y_2];

        
Rxx_1 = [eye(n11);eye(n11)];
Rxx_2 = [eye(n22),zeros(n22);
         eye(n22),(b(2)-a(2))*eye(n22);
         zeros(n22),eye(n22);
         zeros(n22),eye(n22)];
     
K11.Rxx{1} = [Rxx_1 , zeros(size(Rxx_1,1),size(Rxx_2,2));
               zeros(size(Rxx_2,1),size(Rxx_1,2)) , Rxx_2];


Ryy_1 = [eye(n11);eye(n11)];
Ryy_2 = [eye(n22),zeros(n22);
         eye(n22),(b(1)-a(1))*eye(n22);
         zeros(n22),eye(n22);
         zeros(n22),eye(n22)];

K11.Ryy{1} = [Ryy_1 , zeros(size(Ryy_1,1),size(Ryy_2,2));
               zeros(size(Ryy_2,1),size(Ryy_1,2)) , Ryy_2];

K11.dim = K11.dim;


% Map from fundamental state to full bc
K12 = opvar2d();    
K12.I = [a,b];      K12.var1 = var1;  K12.var2 = var2;

R02_1 = [zeros(n11);zeros(n11);zeros(n11);eye(n11)];
R02_2 = polynomial(zeros(16*n22,n22));
R02_2(3*n22+1:4*n22,:) = (b(2)-ss2)*(b(1)-ss1)*eye(n22);
R02_2(7*n22+1:8*n22,:) = (b(2)-ss2)*eye(n22);
R02_2(11*n22+1:12*n22,:) = (b(1)-ss1)*eye(n22);
R02_2(15*n22+1:16*n22,:) = eye(n22);

K12.R02 = [R02_1 , zeros(size(R02_1,1),size(R02_2,2));
          zeros(size(R02_2,1),size(R02_1,2)) , R02_2];

         
Rx2_1 = [zeros(n11);eye(n11)];
Rx2_2 = [zeros(n22);(b(2)-ss2)*eye(n22);zeros(n22);eye(n22)];

K12.Rx2{1} = [Rx2_1 , zeros(size(Rx2_1,1),size(Rx2_2,2));
             zeros(size(Rx2_2,1),size(Rx2_1,2)) , Rx2_2];


Ry2_1 = [zeros(n11);eye(n11)];
Ry2_2 = [zeros(n22);(b(1)-ss1)*eye(n22);zeros(n22);eye(n22)];

K12.Ry2{1} = [Ry2_1 , zeros(size(Ry2_1,1),size(Ry2_2,2));
             zeros(size(Ry2_2,1),size(Ry2_1,2)) , Ry2_2];

K12.dim = K12.dim;


% Map from core bc to full state
K21 = opvar2d();
K21.I = [a,b];      K21.var1 = var1;  K21.var2 = var2;

R20_1 = eye(n11);
R20_2 = [eye(n22),(ss1-a(1))*eye(n22),(ss2-a(2))*eye(n22),(ss1-a(1))*(ss2-a(2))*eye(n22)];

K21.R20 = [R20_1 , zeros(size(R20_1,1),size(R20_2,2));
          zeros(size(R20_2,1),size(R20_1,2)) , R20_2];


R2x_1 = eye(n11);
R2x_2 = [(ss1-tt1)*eye(n22), (ss2-a(2))*(ss1-tt1)*eye(n22)];

K21.R2x{2} = [R2x_1 , zeros(size(R2x_1,1),size(R2x_2,2));
             zeros(size(R2x_2,1),size(R2x_1,2)) , R2x_2];
         

R2y_1 = eye(n11);
R2y_2 = [(ss2-tt2)*eye(n22), (ss1-a(1))*(ss2-tt2)*eye(n22)];

K21.R2y{2} = [R2y_1 , zeros(size(R2y_1,1),size(R2y_2,2));
             zeros(size(R2y_2,1),size(R2y_1,2)) , R2y_2];


K21.dim = K21.dim;


% Map from fundamental state to full state
K22 = opvar2d();
K22.I = [a,b];      K22.var1 = var1;  K22.var2 = var2;

R22_1 = eye(n11);
R22_2 = (ss2-tt2)*(ss1-tt1)*eye(n22);

K22.R22{2,2} = [R22_1 , zeros(size(R22_1,1),size(R22_2,2));
                zeros(size(R22_2,1),size(R22_1,2)) , R22_2];

            
K22.dim = K22.dim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Describe BCs in terms of u_bc and fundamental state:
%
% [0] = [B] [u_bf] = [B o K11] [u_bc] + [B o K12] [\hat{u}_11] = [R] [u_bc] +  [F] [\hat{u}_11]
%                                                 [\hat{u}_22]                     [\hat{u}_22]

R = opvar2d();
R.I = [a,b];    R.var1 = var1;  R.var2 = var2;

R.R00 = Ebb.R00*K11.R00 + int(Ebb.R0x*K11.Rx0,ss1,a(1),b(1)) + int(Ebb.R0y*K11.Ry0,ss2,a(2),b(2));
R.R0x = Ebb.R00*K11.R0x + Ebb.R0x*K11.Rxx{1};
R.R0y = Ebb.R00*K11.R0y + Ebb.R0y*K11.Ryy{1};

R.Rx0 = Ebb.Rx0*K11.R00;
R.Rxy = Ebb.Rx0*K11.R0y + Ebb.Rxy*K11.Ryy{1};
R.Ry0 = Ebb.Ry0*K11.R00;
R.Ryx = Ebb.Ry0*K11.R0x + Ebb.Ryx*K11.Rxx{1};

R.Rxx{1} = Ebb.Rxx{1}*K11.Rxx{1};
R.Rxx{2} = Ebb.Rxx{2}*K11.Rxx{1} + Ebb.Rx0*subs(K11.R0x,ss1,tt1);
R.Rxx{3} = R.Rxx{2};        %<-- Note: we require B.Rxx{2}=B.Rxx{3}

R.Ryy{1} = Ebb.Ryy{1}*K11.Ryy{1};
R.Ryy{2} = Ebb.Ryy{2}*K11.Ryy{1} + Ebb.Ry0*subs(K11.R0y,ss2,tt2);
R.Ryy{3} = R.Ryy{2};        %<-- Note: we require B.Ryy{2}=B.Ryy{3}

R.dim = R.dim;

%%%%%%%%%%

F = opvar2d();
F.I = [a,b];    F.var1 = var1;  F.var2 = var2;

F.R02 = Ebb.R00*K12.R02 + Ebb.R0x*K12.Rx2{1} + Ebb.R0y*K12.Ry2{1};

F.Rx2{1} = Ebb.Rxx{1}*K12.Rx2{1};
F.Rx2{2} = Ebb.Rx0*subs(K12.R02,ss1,tt1) + Ebb.Rxx{2}*subs(K12.Rx2{1},ss1,tt1) + Ebb.Rxy*subs(K12.Ry2{1},ss1,tt1);
F.Rx2{3} = F.Rx2{2};    %<-- Note: we require B.Rxx{2}=B.Rxx{3}

F.Ry2{1} = Ebb.Ryy{1}*K12.Ry2{1};
F.Ry2{2} = Ebb.Ry0*subs(K12.R02,ss2,tt2) + Ebb.Ryy{2}*subs(K12.Ry2{1},ss2,tt2) + Ebb.Ryx*subs(K12.Rx2{1},ss2,tt2);
F.Ry2{3} = F.Ry2{2};    %<-- Note: we require B.Rxx{2}=B.Rxx{3}

F.dim = F.dim;


% % Describe core bc in terms of fundamental state:
%
% [u_bc] = [G] [\hat{u}_11] = - [R^{-1} o F] [\hat{u}_11] =  - [(B o K11)^{-1} o (B o K12)] [\hat{u}_11]
%              [\hat{u}_22]                  [\hat{u}_22]                                   [\hat{u}_22]

Rhat = inv_opvar2d(R);

%%%%%%%%%%

G = opvar2d();
G.I = [a,b];    G.var1 = var1;  G.var2 = var2;

G.R02 = -Rhat.R00*F.R02 - Rhat.R0x*F.Rx2{1} - int(var_swap(Rhat.R0x*F.Rx2{2},ss1,tt1),tt1,a(1),b(1)) ...
                        - Rhat.R0y*F.Ry2{1} - int(var_swap(Rhat.R0y*F.Ry2{2},ss2,tt2),tt2,a(2),b(2));
                   
G.Rx2{1} = -Rhat.Rxx{1}*F.Rx2{1};
G.Rx2{2} = -Rhat.Rx0*subs(F.R02,ss1,tt1) - Rhat.Rxx{1}*F.Rx2{2} - Rhat.Rxx{2}*subs(F.Rx2{1},ss1,tt1) ...
         - int(subs(Rhat.Rxx{2},tt1,eta1)*subs(F.Rx2{2},ss1,eta1),eta1,a(1),b(1)) ...
         - Rhat.Rxy*subs(F.Ry2{1},ss1,tt1) - int(var_swap(Rhat.Rxy*subs(F.Ry2{2},ss1,tt1),ss2,tt2),tt2,a(2),b(2));
G.Rx2{3} = G.Rx2{2};

G.Ry2{1} = -Rhat.Ryy{1}*F.Ry2{1};
G.Ry2{2} = -Rhat.Ry0*subs(F.R02,ss2,tt2) - Rhat.Ryy{1}*F.Ry2{2} - Rhat.Ryy{2}*subs(F.Ry2{1},ss2,tt2) ...
         - int(subs(Rhat.Ryy{2},tt2,eta2)*subs(F.Ry2{2},ss2,eta2),eta2,a(2),b(2)) ...
         - Rhat.Ryx*subs(F.Rx2{1},ss2,tt2) - int(var_swap(Rhat.Ryx*subs(F.Rx2{2},ss2,tt2),ss1,tt1),tt1,a(1),b(1));
G.Ry2{3} = G.Ry2{2};


% % Describe full state in terms of fundamental state:
%
% [u_11] = [K21] [u_bc] + [K22] [\hat{u}_11] = [K21 o G] [\hat{u}_11] + [K22] [\hat{u}_11] = [H + K22] [\hat{u}_11]
% [u_22]                        [\hat{u}_22]             [\hat{u}_22]         [\hat{u}_22]             [\hat{u}_22]

H = opvar2d();
H.I = [a,b];    H.var1 = var1;  H.var2 = var2;

H.R22{2,2} = K21.R20*subs(G.R02,[ss1;ss2],[tt1;tt2]) ...
           + K21.R2x{2}*subs(G.Rx2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(G.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + K21.R2y{2}*subs(G.Ry2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(G.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{3,2} = K21.R20*subs(G.R02,[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(G.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + K21.R2y{2}*subs(G.Ry2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(G.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{2,3} = K21.R20*subs(G.R02,[ss1;ss2],[tt1;tt2]) ...
           + K21.R2x{2}*subs(G.Rx2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(G.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(G.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{3,3} = K21.R20*subs(G.R02,[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(G.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(G.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);

H.dim = H.dim;


% % Construct the T operator
%
% [u] = [T] [\hat{u}] = [I , 0      ] [\hat{u}_0               ]
%                       [0 , H + K22] [[\hat{u}_11; \hat{u}_22]]

T11 = H + K22;

T0 = eye(n00);
Tpde = poly_opvar2d(blkdiag(T0,T11));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Describe PDE states through PI states
%
% \dot{[u_0 ] = A_0 [u_0] + sum_{i,j=1}^2 A_{ij} ([u_11] 
%      [u_11]       [u_1]                         [u_22])_(i,j)
%      [u_22]}      [u_2] 
% \dot{u} = A_{ij} [u]_{(i,j)}
%
% ([u_1]         = [U_{ij}, 0     ] [\hat{u}_{11}]
%  [u_2])_(i,j)    [0     , V_{ij}] [\hat{u}_{22}]


M10 = T11;   
M10.R22{1,2} = subs(M10.R22{2,2} - M10.R22{3,2},T11.var2(1),T11.var1(1));
M10.R22{1,3} = subs(M10.R22{2,3} - M10.R22{3,3},T11.var2(1),T11.var1(1));
M10.R22{2,2} = diff(polynomial(M10.R22{2,2}),T11.var1(1));
M10.R22{3,2} = diff(polynomial(M10.R22{3,2}),T11.var1(1));
M10.R22{2,3} = diff(polynomial(M10.R22{2,3}),T11.var1(1));
M10.R22{3,3} = diff(polynomial(M10.R22{3,3}),T11.var1(1));

M01 = T11;
M01.R22{2,1} = subs(M01.R22{2,2} - M01.R22{2,3},T11.var2(2),T11.var1(2));
M01.R22{3,1} = subs(M01.R22{3,2} - M01.R22{3,3},T11.var2(2),T11.var1(2));
M01.R22{2,2} = diff(polynomial(M01.R22{2,2}),T11.var1(2));
M01.R22{3,2} = diff(polynomial(M01.R22{3,2}),T11.var1(2));
M01.R22{2,3} = diff(polynomial(M01.R22{2,3}),T11.var1(2));
M01.R22{3,3} = diff(polynomial(M01.R22{3,3}),T11.var1(2));

M11 = op_slice(M10,(n11+1:n11+n22),(n11+1:n11+n22));
M11.R22{2,2} = diff(polynomial(M11.R22{2,2}),T11.var1(2));
M11.R22{3,2} = diff(polynomial(M11.R22{3,2}),T11.var1(2));
M11.R22{2,3} = diff(polynomial(M11.R22{2,3}),T11.var1(2));
M11.R22{3,3} = diff(polynomial(M11.R22{3,3}),T11.var1(2));

   
M20 = op_slice(M10,(n11+1:n11+n22),(n00+n11+1:n00+n11+n22)); 
M20.R22{1,2} = subs(M20.R22{2,2} - M20.R22{3,2},T11.var2(1),T11.var1(1));
M20.R22{1,3} = subs(M20.R22{2,3} - M20.R22{3,3},T11.var2(1),T11.var1(1));
M20.R22{2,2} = diff(polynomial(M20.R22{2,2}),T11.var1(1));
M20.R22{3,2} = diff(polynomial(M20.R22{3,2}),T11.var1(1));
M20.R22{2,3} = diff(polynomial(M20.R22{2,3}),T11.var1(1));
M20.R22{3,3} = diff(polynomial(M20.R22{3,3}),T11.var1(1));

M02 = op_slice(M01,(n11+1:n11+n22),(n11+1:n11+n22));  
M02.R22{2,1} = subs(M02.R22{2,2} - M02.R22{2,3},T11.var2(2),T11.var1(2));
M02.R22{3,1} = subs(M02.R22{3,2} - M02.R22{3,3},T11.var2(2),T11.var1(2));
M02.R22{2,2} = diff(polynomial(M02.R22{2,2}),T11.var1(2));
M02.R22{3,2} = diff(polynomial(M02.R22{3,2}),T11.var1(2));
M02.R22{2,3} = diff(polynomial(M02.R22{2,3}),T11.var1(2));
M02.R22{3,3} = diff(polynomial(M02.R22{3,3}),T11.var1(2));

M21 = M11;   
M21.R22{1,2} = subs(M11.R22{2,2} - M11.R22{3,2},T11.var2(1),T11.var1(1));
M21.R22{1,3} = subs(M11.R22{2,3} - M11.R22{3,3},T11.var2(1),T11.var1(1));
M21.R22{2,2} = diff(polynomial(M11.R22{2,2}),T11.var1(1));
M21.R22{3,2} = diff(polynomial(M11.R22{3,2}),T11.var1(1));
M21.R22{2,3} = diff(polynomial(M11.R22{2,3}),T11.var1(1));
M21.R22{3,3} = diff(polynomial(M11.R22{3,3}),T11.var1(1));

M12 = M11;   
M12.R22{2,1} = subs(M11.R22{2,2} - M11.R22{2,3},T11.var2(2),T11.var1(2));
M12.R22{3,1} = subs(M11.R22{3,2} - M11.R22{3,3},T11.var2(2),T11.var1(2));
M12.R22{2,2} = diff(polynomial(M11.R22{2,2}),T11.var1(2));
M12.R22{3,2} = diff(polynomial(M11.R22{3,2}),T11.var1(2));
M12.R22{2,3} = diff(polynomial(M11.R22{2,3}),T11.var1(2));
M12.R22{3,3} = diff(polynomial(M11.R22{3,3}),T11.var1(2));

opvar2d I_11 O_10 O_20 I_22;
I_11.I = dom; O_10.I = dom; O_20.I = dom; I_22.I = dom;
I_11.dim = [0,0;0,0;0,0;n11,n00+n11]; 
O_10.dim = [0,0;0,0;0,0;n11+n22,n00];
O_20.dim = [0,0;0,0;0,0;n22,n00+n11];  
I_22.dim = [0,0;0,0;0,0;n22,n22];

I_11.R22{1,1} = [zeros(n11,n00) , eye(n11)];
I_22.R22{1,1} = eye(n22);

M10 = [O_10 , M10];
M01 = [O_10 , M01];
M11 = blkdiag(I_11,M11); % Exploit the fact that u_11xy = hat{u_11} to construct M11
   
M20 = [O_20 , M20];
M02 = [O_20 , M02];
M21 = [O_20 , M21];
M12 = [O_20 , M12];

M22 = [O_20 , I_22]; % Exploit the fact that u_22xxyy = hat{u_22} to construct M22


% % Describe the right-hand side of the PDE using PI operator
%
% T \dot{\hat{u}} = A \hat{u}

A00 = Am{1,1}*Tpde;
A10 = Am{2,1}*M10;
A01 = Am{1,2}*M01;
A11 = Am{2,2}*M11;
A20 = Am{3,1}*M20;
A02 = Am{1,3}*M02;
A21 = Am{3,2}*M21;
A12 = Am{2,3}*M12;
A22 = Am{3,3}*M22;

Apde = A00 + A10 + A01 + A11 + A20 + A02 + A21 + A12 + A22;
Apde = poly_opvar2d(Apde);
Apde.dim = Apde.dim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Include ODE dynamics
if nx~=0

Iode = opvar2d([nx,nx;0,0;0,0;0,0]);
Iode.I = dom;   Iode.var1 = var1;  Iode.var2 = var2;
Iode.R00 = eye(nx);

Aode = opvar2d([nx,nx;0,0;0,0;0,0]);
Aode.I = dom;   Aode.var1 = var1;  Aode.var2 = var2;
Aode.R00 = PDE.ODE.A;

T = poly_opvar2d(blkdiag(Iode,Tpde));
A = poly_opvar2d(blkdiag(Aode,Apde));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Include ODE-PDE interconnections
if nx~=0 && ~all(PDE.n.n_pde==0,'all')

% Contribution from ODE to PDE state
A.R20 = PDE.PDE.Bpv*PDE.ODE.Cv;

% Contribution from boundary to ODE state
n_pde_sum = [0,cumsum(PDE.n.n_pde(:))'];
R02 = A.R02;
for m=1:numel(PDE.PDE.Drb)
    tmp = PDE.PDE.Drb{m};
    
    % Extract elements of A associated to desired state component
    Rstate = tmp.Rstate;
    Rloc = Rstate*[1;3]+1;
    indcs = n_pde_sum(Rloc)+1:n_pde_sum(Rloc+1);
    if ~isempty(indcs)
        Tr = Tpde(:,indcs);

        deg = tmp.D';   % Order of derivatives to take
        loc = tmp.delta+1;  % location to evaluate: upper or lower boundary
        dl1 = dom(1,loc(1));    dl2 = dom(2,loc(2));
        % Differentiate the A operator to desired degree
        Tdif = diff_opvar2d(Tr,var1.^deg,'pure');
        % Evaluate the A operator at the desired boundary
        Tdel = delta_opvar2d(Tdif,var1,[dl1;dl2],'pure');
        % Add the PDE to ODE contribution to A
        R02(:,indcs) = R02(:,indcs) + PDE.ODE.Bxr*tmp.coeff*Tdel.R02;
    end
end
A.R02 = R02;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get rid of terms below tolerance, and define the PIE structure

PIE.T = zremove(T,1e-14);
PIE.A = zremove(A,1e-14);

end