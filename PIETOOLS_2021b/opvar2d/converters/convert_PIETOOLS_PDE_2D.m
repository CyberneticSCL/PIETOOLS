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
% \dot{[X_0] = A{1,1} [X_0] + A{2,1} [X_1x] + A{1,2} [X_1y] + A{2,2} [X_1xy] + A{3,1} [X_2xx] + A{1,3} [X_2yy] + A{3,2} [X_2xxy] + A{2,3} [X_2xyy] + A{3,3} [X_2xxyy]
%      [X_1]          [X_1]          [X_2x]          [X_2y]          [X_2xy]
%      [X_2]}         [X_2]
%        Must have fields A, Bpv and Drb, where Bpv describes ODE->PDE
%        interconnection, \dot{u} = ... + Bpv*v, and Drb describes PDE->ODE
%        interconnection, r = Drb*u_bf.
% - BC: Structure with field Ebb, which must be an opvar2d object
%       specifying BCs as 0 = Ebb*ubf
%
% Here:
% X_bf = [X_1(a1,a2); X_1(b1,a2); X_1(a1,b2); X_1(b1,b2); 
%         X_2(a1,a2); x_2(b1,a2); X_2(a1,b2); X_2(b1,b2); X_2x(a1,a2); X_2x(b1,a2); X_2x(a1,b2); X_2x(b1,b2);
%         X_2y(a1,a2); X_2y(b1,a2); X_2y(a1,b2); X_2y(b1,b2); X_2xy(a1,a2); X_2xy(b1,a2); X_2xy(a1,b2); X_2xy(b1,b2);
%         X_1x(x,a2); X_1x(x,b2); 
%         X_2xx(x,a2); X_2xx(x,b2); X_2xxy(x,a2); X_2xxy(x,b2);
%         X_1y(a1,y); X_1y(b1,y);
%         X_2yy(a1,y); X_2yy(b1,y); X_2xyy(a1,y); X_2xyy(b1,y)]
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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 12/30/2021
% 02/18/2022 - DJ: Update to include inputs and outputs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Initialization

% Initialize PDE in case this has not been done
PDE = initialize_PIETOOLS_PDE_2D(PDE);

% Extract dimensions of PDE state components
nx = PDE.n.nx;
nz = PDE.n.nz;      nw = PDE.n.nw;      nv = PDE.n.nv;
ny = PDE.n.ny;      nu = PDE.n.nu;      nr = PDE.n.nr;
n_pde = PDE.n.n_pde;
if any(any(n_pde-diag(diag(n_pde))))
    error('Only PDEs with equal differentiability in both spatial variables are currently supported')
end
n00 = n_pde(1,1);   n11 = n_pde(2,2);   n22 = n_pde(3,3);
np = sum(sum(n_pde));

% Extract lower and  upper boundaries of spatial domain
dom = PDE.dom;
a = dom(:,1);   b = dom(:,2);

% Initialize PIE variables, as well as dummy variables eta
vars = PDE.vars;
var1 = vars(:,1);       var2 = vars(:,2);
ss1 = var1(1);          tt1 = var2(1);
ss2 = var1(2);          tt2 = var2(2);
pvar eta1 eta2;


% % Extract the number of signals

% ODE states
if numel(nx)~=1
    error('Number of ODE signals should be specified as a scalar integer')
else
    nx_op = [nx;0;0;0];     % ODE signals do not vary in time
end

% PDE state
np_op = [0;0;0;np];     % We assume only PDEs in 2 spatial variables

% Exogenous inputs
if numel(nw)==1
    nw_op = [nw;0;0;0]; % Assume w does not vary in space
elseif numel(nw)==4
    nw_op = nw(:);
    nw = sum(nw_op);
else
    error('Number of exogenous inputs should be specified as scalar integer or 4x1 array');
end

% Actuator inputs
if numel(nu)==1
    nu_op = [nu;0;0;0]; % Assume u does not vary in space
elseif numel(nu)==4
    nu_op = nu(:);
    nu = sum(nu_op);
else
    error('Number of actuator inputs should be specified as scalar integer or 4x1 array');
end

% Regulated outputs
if numel(nz)==1
    nz_op = [nz;0;0;0]; % Assume z does not vary in space
elseif numel(nz)==4
    nz_op = nz(:);
    nz = sum(nz_op);
else
    error('Number of regulated outputs should be specified as scalar integer or 4x1 array');
end

% Regulated outputs
if numel(ny)==1
    ny_op = [ny;0;0;0]; % Assume y does not vary in space
elseif numel(ny)==4
    ny_op = ny(:);
    ny = sum(ny_op);
else
    error('Number of regulated outputs should be specified as scalar integer or 4x1 array');
end

% PDE-2-ODE signals
if numel(nv)==1
    nv_op = [nv;0;0;0]; % Assume v does not vary in space
elseif numel(nv)==4
    nv_op = nv(:);
    nv = sum(nv_op);
else
    error('Number of PDE to ODE signals should be specified as scalar integer or 4x1 array');
end

% ODE-2-PDE signals
if numel(nr)==1
    nr_op = [nr;0;0;0]; % Assume r does not vary in space
elseif numel(nr)==4
    nr_op = nr(:);
    nr = sum(nr_op);
else
    error('Number of PDE to ODE signals should be specified as scalar integer or 4x1 array');
end

% Define the number of boundary conditions
nBC_op = [n11+4*n22 ; n11+2*n22; n11+2*n22; 0];
nBF_op = [4*n11+16*n22; 2*n11+4*n22; 2*n11+4*n22; 0];
np12_op = [0;0;0;n11+n22];
np0_op = [0;0;0;n00];



% Extract BCs defining operator Ebb
Ebx = opvar2d(PDE.BC.Ebv*PDE.ODE.Cv,[nBC_op,nx_op],dom,vars);
Ebb = PDE.BC.Ebb;
Ebw = opvar2d(PDE.BC.Ebv*PDE.ODE.Dvw,[nBC_op,nw_op],dom,vars);
Ebu = opvar2d(PDE.BC.Ebv*PDE.ODE.Dvu,[nBC_op,nu_op],dom,vars);
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Describe full bc and full state in terms of core bc and fundamental state:
%
% [X_bf] = [K11] [X_bc] + [K12] [\hat{X}_11]
%                               [\hat{X}_22]
%
% [X_11] = [K21] [X_bc] + [K22] [\hat{X}_11]
% [X_22]                        [\hat{X}_22]

% Map from core bc to full bc
K11 = opvar2d([],[nBF_op,nBC_op],dom,vars);
% K11 = opvar2d();
% K11.I = [a,b];      K11.var1 = var1;  K11.var2 = var2;

R00_1 = [eye(n11);eye(n11);eye(n11);eye(n11)];
R00_2 = kron([1,0,0,0;1,b(1)-a(1),0,0;1,0,b(2)-a(2),0;1,b(1)-a(1),b(2)-a(2),(b(2)-a(2))*(b(1)-a(1));
              0,1,0,0;0,1,0,0;0,1,0,b(2)-a(2);0,1,0,b(2)-a(2);
              0,0,1,0;0,0,1,b(1)-a(1);0,0,1,0;0,0,1,b(1)-a(1);
              0,0,0,1;0,0,0,1;0,0,0,1;0,0,0,1],eye(n22));
          
K11.R00 = [R00_1 , zeros(size(R00_1,1),size(R00_2,2));
            zeros(size(R00_2,1),size(R00_1,2)) , R00_2];

        
R0x_1 = [zeros(n11);eye(n11);zeros(n11);eye(n11)];        
R0x_2 = polynomial(zeros(16*n22,2*n22));
if n22>0
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
end

K11.R0x = [R0x_1 , zeros(size(R0x_1,1),size(R0x_2,2));
            zeros(size(R0x_2,1),size(R0x_1,2)) , R0x_2];
    
        
R0y_1 = [zeros(n11);zeros(n11);eye(n11);eye(n11)];
R0y_2 = polynomial(zeros(16*n22,2*n22));
if n22>0
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
end

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
K12 = opvar2d([],[nBF_op,np12_op],dom,vars);
% K12 = opvar2d();    
% K12.I = [a,b];      K12.var1 = var1;  K12.var2 = var2;

R02_1 = [zeros(n11);zeros(n11);zeros(n11);eye(n11)];
R02_2 = polynomial(zeros(16*n22,n22));
if n22>0
    R02_2(3*n22+1:4*n22,:) = (b(2)-ss2)*(b(1)-ss1)*eye(n22);
    R02_2(7*n22+1:8*n22,:) = (b(2)-ss2)*eye(n22);
    R02_2(11*n22+1:12*n22,:) = (b(1)-ss1)*eye(n22);
    R02_2(15*n22+1:16*n22,:) = eye(n22);
end

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
K21 = opvar2d([],[np12_op,nBC_op],dom,vars);
% K21 = opvar2d();
% K21.I = [a,b];      K21.var1 = var1;  K21.var2 = var2;

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
K22 = opvar2d([],[np12_op,np12_op],dom,vars);
% K22 = opvar2d();
% K22.I = [a,b];      K22.var1 = var1;  K22.var2 = var2;

R22_1 = eye(n11);
R22_2 = (ss2-tt2)*(ss1-tt1)*eye(n22);

K22.R22{2,2} = [R22_1 , zeros(size(R22_1,1),size(R22_2,2));
                zeros(size(R22_2,1),size(R22_1,2)) , R22_2];

            
K22.dim = K22.dim;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Describe BCs in terms of X_bc and fundamental state:
%
% 0 = Ebb*X_bf + Ebx*x + Ebw*w + Ebu*u 
%   = (Ebb*K11)*X_bc + (Ebb*K12)*[\hat{X}_11;\hat{X}_22] + Ebx*x + Ebw*w + Ebu*u 
%   = R*X_bc + F*[\hat{X}_11;\hat{X}_22] + Ebx*x + Ebw*w + Ebu*u

R = opvar2d([],[nBC_op,nBC_op],dom,vars);

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

F = opvar2d([],[nBC_op,np12_op],dom,vars);

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
% 0 = Ebb*X_bf + Ebx*x + Ebw*w + Ebu*u 
%   = (Ebb*K11)*X_bc + (Ebb*K12)*[\hat{X}_11;\hat{X}_22] + Ebx*x + Ebw*w + Ebu*u 
% X_bc  = Gbb*[\hat{X}_11;\hat{X}_22] + Gbx*x + Gbw*w + Gbu*u
%       = (-R^{-1}*F)*[\hat{X}_11;\hat{X}_22] + (-R^{-1}*Ebx)*x + (-R^{-1}*Ebw)*w + (-R^{-1}*Ebu)*u
%

Rhat = inv_opvar2d(R);

%%%%%%%%%%

Gbb = opvar2d([],[nBC_op,np12_op],dom,vars);

Gbb.R02 = -Rhat.R00*F.R02 - Rhat.R0x*F.Rx2{1} - int(var_swap(Rhat.R0x*F.Rx2{2},ss1,tt1),tt1,a(1),b(1)) ...
                        - Rhat.R0y*F.Ry2{1} - int(var_swap(Rhat.R0y*F.Ry2{2},ss2,tt2),tt2,a(2),b(2));
                   
Gbb.Rx2{1} = -Rhat.Rxx{1}*F.Rx2{1};
Gbb.Rx2{2} = -Rhat.Rx0*subs(F.R02,ss1,tt1) - Rhat.Rxx{1}*F.Rx2{2} - Rhat.Rxx{2}*subs(F.Rx2{1},ss1,tt1) ...
         - int(subs(Rhat.Rxx{2},tt1,eta1)*subs(F.Rx2{2},ss1,eta1),eta1,a(1),b(1)) ...
         - Rhat.Rxy*subs(F.Ry2{1},ss1,tt1) - int(var_swap(Rhat.Rxy*subs(F.Ry2{2},ss1,tt1),ss2,tt2),tt2,a(2),b(2));
Gbb.Rx2{3} = Gbb.Rx2{2};

Gbb.Ry2{1} = -Rhat.Ryy{1}*F.Ry2{1};
Gbb.Ry2{2} = -Rhat.Ry0*subs(F.R02,ss2,tt2) - Rhat.Ryy{1}*F.Ry2{2} - Rhat.Ryy{2}*subs(F.Ry2{1},ss2,tt2) ...
         - int(subs(Rhat.Ryy{2},tt2,eta2)*subs(F.Ry2{2},ss2,eta2),eta2,a(2),b(2)) ...
         - Rhat.Ryx*subs(F.Rx2{1},ss2,tt2) - int(var_swap(Rhat.Ryx*subs(F.Rx2{2},ss2,tt2),ss1,tt1),tt1,a(1),b(1));
Gbb.Ry2{3} = Gbb.Ry2{2};


Gbx = -Rhat*Ebx;
Gbw = -Rhat*Ebw;
Gbu = -Rhat*Ebu;

% % Describe full state in terms of fundamental state:
%
% [X_11] = [K21] [X_bc] + [K22] [\hat{X}_11] = [K21 o G] [\hat{X}_11] + [K22] [\hat{X}_11] = [H + K22] [\hat{X}_11]
% [X_22]                        [\hat{X}_22]             [\hat{X}_22]         [\hat{X}_22]             [\hat{X}_22]

% [X_11; X_22] = K21*X_bc + K22*[\hat{X}_11;\hat{X}_22]
%              = (K21*Gbb)*[\hat{X}_11;\hat{X}_22] + (K21*Gbx)*x + (K21*Gbw)*w + (K21*Gbu)*u + K22*[\hat{X}_11;\hat{X}_22]
%              = (H + K22)*[\hat{X}_11;\hat{X}_22] + T_x2p*x + Tw*w + Tu*u

H = opvar2d([],[np12_op,np12_op],dom,vars);

H.R22{2,2} = K21.R20*subs(Gbb.R02,[ss1;ss2],[tt1;tt2]) ...
           + K21.R2x{2}*subs(Gbb.Rx2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(Gbb.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + K21.R2y{2}*subs(Gbb.Ry2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(Gbb.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{3,2} = K21.R20*subs(Gbb.R02,[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(Gbb.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + K21.R2y{2}*subs(Gbb.Ry2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(Gbb.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{2,3} = K21.R20*subs(Gbb.R02,[ss1;ss2],[tt1;tt2]) ...
           + K21.R2x{2}*subs(Gbb.Rx2{1},[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(Gbb.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(Gbb.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);
H.R22{3,3} = K21.R20*subs(Gbb.R02,[ss1;ss2],[tt1;tt2]) ...
           + int(subs(K21.R2x{2},tt1,eta1)*subs(Gbb.Rx2{2},[ss1;ss2],[eta1;tt2]),eta1,a(1),ss1) ...
           + int(subs(K21.R2y{2},tt2,eta2)*subs(Gbb.Ry2{2},[ss1;ss2],[tt1;eta2]),eta2,a(2),ss2);

H.dim = H.dim;


% % Construct the T operator
%
% X = T_pde*\hat{X} + T_o2p*x + Tw_pde*w + Tu_pde*u
% T_pde*\hat{X} = [I , 0 ; 0, H + K22] * [\hat{X}_00; [\hat{X}_11; \hat{X}_22]]

T11 = H + K22;
Tx12 = K21*Gbx;
Tw12 = K21*Gbw;
Tu12 = K21*Gbu;

T0 = opvar2d(eye(n00),[np0_op,np0_op],dom,vars);
T_pde = poly_opvar2d(blkdiag(T0,T11));
iszero_Tpde = T_pde==0;

Tx0 = opvar2d(zeros([n00,nx]),[np0_op,nx_op],dom,vars);
T_o2p = poly_opvar2d([Tx0;Tx12]);
iszero_To2p = T_o2p==0;

Tw0 = opvar2d(zeros([n00,nw]),[np0_op,nw_op],dom,vars);
Tw_pde = poly_opvar2d([Tw0;Tw12]);
iszero_Tw = Tw_pde==0;

Tu0 = opvar2d(zeros([n00,nu]),[np0_op,nu_op],dom,vars);
Tu_pde = poly_opvar2d([Tu0;Tu12]);
iszero_Tu = Tu_pde==0;


% Define the operators T, Tw and Tu such that
% [x;X] = T*[x;\hat{X}] + Tw*w + Tu*u

% T = [I,0; T_o2p, T_pde];
I_ode = opvar2d(eye(nx),[nx_op,nx_op],dom,vars);
O_p2o = opvar2d(zeros(nx,np),[nx_op,np_op],dom,vars);
T = poly_opvar2d([I_ode, O_p2o; T_o2p, T_pde]);

% Tw = [0;Tw_pde]
Ow = opvar2d(zeros(nx,nw),[nx_op,nw_op],dom,vars);
Tw = poly_opvar2d([Ow; Tw_pde]);

% Tw = [0;Tu_pde]
Ou = opvar2d(zeros(nx,nu),[nx_op,nu_op],dom,vars);
Tu = poly_opvar2d([Ou; Tu_pde]);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Define ODE-PDE interconnection dynamics
% Express system in the form
%
% [\dot{x}] = [ A_ode , A_p2o , Bw_ode, Bu_ode] [ x ]
% [\dot{X}]   [ A_o2p , A_pde , Bw_pde, Bu_pde] [ X ]
% [     z ]   [ Cz_ode, Cz_pde, Dzw   , Dzu   ] [ w ]
% [     y ]   [ Cy_ode, Cy_pde, Dyw   , Dyu   ] [ u ]
%
% Where all the operators are PI operators

% % Initialize first the maps using just ODE input signals
% ODE to output maps:
A_ode = opvar2d(PDE.ODE.A,[nx_op,nx_op],dom,vars);  
Cz_ode = opvar2d(PDE.ODE.Cz,[nz_op,nx_op],dom,vars);
Cy_ode = opvar2d(PDE.ODE.Cy,[ny_op,nx_op],dom,vars);

% w to output maps
Bw_ode = opvar2d(PDE.ODE.Bw,[nx_op,nw_op],dom,vars);
Dzw = opvar2d(PDE.ODE.Dzw,[nz_op,nw_op],dom,vars);
Dyw = opvar2d(PDE.ODE.Dyw,[ny_op,nw_op],dom,vars);

% u to output maps
Bu_ode = opvar2d(PDE.ODE.Bu,[nx_op,nu_op],dom,vars);
Dzu = opvar2d(PDE.ODE.Dzu,[nz_op,nu_op],dom,vars);
Dyu = opvar2d(PDE.ODE.Dyu,[ny_op,nu_op],dom,vars);


% % For the PDE dynamics, we have to consider the signals r and v
% In the PDE, we have
%
% \dot{X} = \sum_ij A_pde{i,j} * d/dx^i d/dy^j * Sij * X 
%            + \sum_ijkl Bpb{i,j,k,l} delta_k(x) delta_l(y) * d/dx^i d/dy^j * Sij * X
%               + Bpv * v
%       r = \sum_ij Crp{i,j} * d/dx^i d/dy^j * Sij * X 
%            + \sum_ijkl Drb{i,j,k,l} Delta_k(x) Delta_l(y) * d/dx^i d/dy^j * Sij * X
%               + Drv * v % <-- We assume Drv=0
%
% Which becomes:
%
% T*[x;\dot{\hat{X}}] + Tw*\dot{w} + Tu*\dot{u}
%       = \sum_ij A_pde(i,j) * d/dx^i d/dy^j * Sij * [T*[x;\hat{X}] + Tw*w + Tu*u]
%            + \sum_ijkl Bpb(i,j,k,l) * delta_k(x) delta_l(y) * d/dx^i d/dy^j * Sij * [T*[x;\hat{X}] + Tw*w + Tu*u]
%               + Bpv * [Cv*x + Dvw*w + Dvu*u]
%     r = \sum_ij Crp{i,j} * d/dx^i d/dy^j * Sij * [T*[x;\hat{X}] + Tw*w + Tu*u]
%          + \sum_ijkl Drb{i,j,k,l} Delta_k(x) Delta_l(y) * d/dx^i d/dy^j * Sij * [T*[x;\hat{X}] + Tw*w + Tu*u]
%

% We initialize some immediate contributions
A_o2p = opvar2d(PDE.PDE.Bpv*PDE.ODE.Cv,[np_op,nx_op],dom,vars);
Bw_pde = opvar2d(PDE.PDE.Bpv*PDE.ODE.Dvw,[np_op,nw_op],dom,vars);
Bu_pde = opvar2d(PDE.PDE.Bpv*PDE.ODE.Dvu,[np_op,nu_op],dom,vars);

Br = opvar2d(PDE.ODE.Br,[nx_op,nr_op],dom,vars);
Dzr = opvar2d(PDE.ODE.Dzr,[nz_op,nr_op],dom,vars);
Dyr = opvar2d(PDE.ODE.Dyr,[ny_op,nr_op],dom,vars);


% % For the remaining contributions, we will need to compute the composition
% of PI operators T, Tw, Tu with the differential operators d/dx^i d/dy^j,
% and the delta operators Drb{i,j,k,l} Delta_k(x) Delta_l(y).

% Initialize contributions as zero
A_pde = opvar2d([],[np_op,np_op],dom,vars);     % Map from PDE state to PDE state
Crp = opvar2d([],[nr_op,np_op],dom,vars);       % Map from PDE state to PDE-to-ODE signal
Crx = opvar2d([],[nr_op,nx_op],dom,vars);       % Map from ODE state to PDE-to-ODE signal
Drw = opvar2d([],[nr_op,nw_op],dom,vars);       % Map from input w to PDE-to-ODE signal
Dru = opvar2d([],[nr_op,nu_op],dom,vars);       % Map from input u to PDE-to-ODE signal


% % For the derivatives, we will need to be able to extract only those
% variable differentiable up to the desired order. To achieve this,
% we build a cell such that state component j corresponds to indcs 
% "state_indcs{j}" of the full PDE state X.
state_indx = cumsum(PDE.n.n_pde(:))';
state_indcs = [1,state_indx(1:end-1)+1;state_indx];
state_indcs = mat2cell(state_indcs,2,ones(1,length(state_indcs)));
state_indcs = cellfun(@(m) (m(1):m(2))',state_indcs,'UniformOutput',false);
state_indcs = reshape(state_indcs,numel(PDE.n.n_pde),1);

% In addition, the following arrays will be used to determine the
% appropriate state component and associated variables
sz_Np = [1,cumprod(size(n_pde))];
n_pde_sum = [0,cumsum(n_pde(:))'];

% We define a table of degrees, so that state component j is differentiable
% up to degree deg_table(j,:)
%N = size(PDE.n.n_pde)-1;
%Ncell = num2cell(N);
%degvals_cell = cellfun(@(m) (0:m)',Ncell,'uni',0);
%degvals_grid = degvals_cell;
%[degvals_grid{:}] = ndgrid(degvals_cell{:});
%deg_table = cell2mat(cellfun(@(m)m(:),degvals_grid,'uni',0));


% % Now, we sum over each term in the PDE, and add the contribution to the 
% appropriate operators
% dST_pde_cell = cell(size(Apde_m));
% dST_o2p_cell = cell(size(Apde_m));
% dSTw_cell = cell(size(Apde_m));
% dSTu_cell = cell(size(Apde_m));
for m=1:numel(PDE.PDE.A)
    Apde_m = PDE.PDE.A{m};
    Crp_m = PDE.PDE.Crp{m};
    
    ddeg = Apde_m.D';               % Degree of derivative considered in term m
    Rstates = Apde_m.Rstate(:);     % Which state components does App_m map?
    Rlocs = cell2mat(Rstates)*sz_Np(1:end-1)' + 1; % Linear indices associated to state components
    Rindcs = cell2mat(state_indcs(Rlocs));  % Linear indices of variables of state components
    npR = sum(n_pde(Rlocs));        % Number of state variables associated to state components Rstates
    npR_op = [0;0;0;npR];           % Associated column dimension of the opvar2d
    
    if ~isempty(Rindcs)
        % Compute coefficients mapping derivative "ddeg" of states "Rstates"
        % to the PDE state and output signal r
        
        % First, check that the coefficients are not zero
        Apde_m = cell2mat_poly(Apde_m.coeff);
        iszero_Apde = all(all(isequal(Apde_m,0)));
        Crp_m = cell2mat_poly(Crp_m.coeff);
        iszero_Crp = all(all(isequal(Crp_m,0)));
        
        % If coefficients are zero, there's no need to process this term
        if ~iszero_Apde || ~iszero_Crp
            % Convert polynomials to opvar2d objects
            Apde_m = opvar2d(Apde_m,[np_op,npR_op],dom,vars);
            Crp_m = opvar2d(Crp_m,[nr_op,npR_op],dom,vars);
             
            % % Process PDE-PDE contributions
            if ~iszero_Tpde
                % Extract rows of T operators associated to components Rstates
                ST_pde = T_pde(Rindcs,:);
                % Compute composition of PI operator T with differential operator
                % (d/dx^ddeg(1) * d/dy^ddeg(2))
                dST_pde = diff_opvar2d(ST_pde,var1.^ddeg,'pure');
                % Add contribution of each input to PDE dynamics
                if ~iszero_Apde
                    A_pde = A_pde + Apde_m*dST_pde;
                end
                % Add contribution of internal PDE state to output signal r
                if ~iszero_Crp
                    Crp = Crp + Crp_m*dST_pde;
                end
            end
            
            % % Repeat for contributions from ODE state
            if nx~=0 && ~iszero_To2p
                ST_o2p = T_o2p(Rindcs,:);
                dST_o2p = diff_opvar2d(ST_o2p,var1.^ddeg,'pure');
                if ~iszero_Apde
                    A_o2p = A_o2p + Apde_m*dST_o2p;
                end
                if ~iszero_Crp
                    Crx = Crx + Crp_m*dST_o2p;
                end
            end
            
            % % Repeat for contributions from input w
            if nw~=0 && ~iszero_Tw
                STw = Tw_pde(Rindcs,:);
                dSTw = diff_opvar2d(STw,var1.^ddeg,'pure');
                if ~iszero_Apde
                    Bw_pde = Bw_pde + Apde_m*dSTw;
                end
                if ~iszero_Crp
                    Drw = Drw + Crp_m*dSTw;
                end
            end
            
            % % Repeat for contributions from input u
            if nu~=0 && ~iszero_Tu
                STu = Tu_pde(Rindcs,:);
                dSTu = diff_opvar2d(STu,var1.^ddeg,'pure');
                if ~iszero_Apde
                    Bu_pde = Bu_pde + Apde_m*dSTu;
                end
                if ~iszero_Crp
                    Dru = Dru + Crp_m*dSTu;
                end
            end
        end
    end
end

% % Next: Add contributions from boundary PDE states
% NOTE: Only corner boundary contributions (no edge boundaries) are currently supported
for m=1:numel(PDE.PDE.Drb)
    Bpb_m = PDE.PDE.Bpb{m};
    Drb_m = PDE.PDE.Drb{m};
     
    % Extract elements of A associated to desired state component
    ddeg = Bpb_m.D';                    % Degree of derivative considered in term m
    loc = Bpb_m.delta+1;                % Location to evaluate: upper or lower boundary
    dl1 = dom(1,loc(1));    dl2 = dom(2,loc(2));    % Associated boundary values
    
    Rstate = Bpb_m.Rstate;
    Rloc = Rstate*sz_Np(1:end-1)' + 1;  % Linear index associated to state component
    Rindcs = n_pde_sum(Rloc)+1:n_pde_sum(Rloc+1);    % Linear indices of variables of state component
    nbR = sum(n_pde(Rloc));             % Number of state variables associated to state components Rstates
    nbR_op = [nbR;0;0;0];               % Associated column dimension of the opvar2d
    
    if ~isempty(Rindcs)
        % Compute coefficients mapping derivative "ddeg" of states "Rstate"
        % at boundary "[dl1,dl2]" to the PDE state and output signal r
        
        % First, check that the coefficients are not zero
        Bpb_m = cell2mat_poly(Bpb_m.coeff);
        iszero_Bpb = all(all(isequal(Bpb_m,0)));
        Drb_m = Drb_m.coeff;
        iszero_Drb = all(all(isequal(Drb_m,0)));
        
        % If coefficients are zero, there's no need to process this term
        if ~iszero_Bpb || ~iszero_Drb
            % Convert polynomials to opvar2d objects
            Bpb_m = opvar2d(Bpb_m,[np_op,nbR_op],dom,vars);
            Drb_m = opvar2d(Drb_m,[nr_op,nbR_op],dom,vars);
            
            % % Process PDE-BC contributions
            if ~iszero_Tpde
                % Extract rows of T operators associated to component Rstate
                ST_pde = T_pde(Rindcs,:);
                % Compute composition of PI operators T with differential operator
                % (d/dx^ddeg(1) * d/dy^ddeg(2))
                T_pde_dif = diff_opvar2d(ST_pde,var1.^ddeg,'pure');
                % Take composition of PI operators T with delta operator
                % (Delta(x=dl1) * Delta(y=dl2))
                T_pde_del = delta_opvar2d(T_pde_dif,var1,[dl1;dl2],'pure');
                % Add boundary contribution from each input to PDE state evolution
                if ~iszero_Bpb
                    A_pde = A_pde + Bpb_m*T_pde_del;
                end
                % Add boundary contribution from each input to r
                if ~iszero_Drb
                    Crp = Crp + Drb_m*T_pde_del;
                end
            end
            
            % % Repeat for contributions from ODE state
            if nx~=0 && ~iszero_To2p
                ST_o2p = T_o2p(Rindcs,:);
                T_o2p_dif = diff_opvar2d(ST_o2p,var1.^ddeg,'pure');
                T_o2p_del = delta_opvar2d(T_o2p_dif,var1,[dl1;dl2],'pure');
                if ~iszero_Bpb
                    A_o2p = A_o2p + Bpb_m*T_o2p_del;
                end
                if ~iszero_Drb
                    Crx = Crx + Drb_m*T_o2p_del;
                end
            end
            
            % % Repeat for contributions from input w
            if nw~=0 && ~iszero_Tw
                STw = Tw_pde(Rindcs,:);
                Tw_dif = diff_opvar2d(STw,var1.^ddeg,'pure');
                Tw_del = delta_opvar2d(Tw_dif,var1,[dl1;dl2],'pure');
                if ~iszero_Bpb
                    Bw_pde = Bw_pde + Bpb_m*Tw_del;
                end
                if ~iszero_Drb
                    Drw = Drw + Drb_m*Tw_del;
                end
            end
            
            % % Repeat for contributions from input u
            if nu~=0 && ~iszero_Tu
                STu = Tu_pde(Rindcs,:);
                Tu_dif = diff_opvar2d(STu,var1.^ddeg,'pure');
                Tu_del = delta_opvar2d(Tu_dif,var1,[dl1;dl2],'pure');
                if iszero_Bpb
                    Bu_pde = Bu_pde + Bpb_m*Tu_del;
                end
                if ~iszero_Drb
                    Dru = Dru + Drb_m*Tu_del;
                end
            end
            
        end
    end
end
% Now that we have collected all terms (all derivatives, all delta
% operators) into respective PI operators, we can define the remaining
% input-state-output maps

% Include contribution of PDE output r to ODE state
A_ode = A_ode + Br*Crx;
A_p2o = Br*Crp;
A = poly_opvar2d([A_ode,A_p2o; A_o2p, A_pde]);

Bw_ode = Bw_ode + Br*Drw;
Bw = poly_opvar2d([Bw_ode;Bw_pde]);

Bu_ode = Bu_ode + Br*Dru;
Bu = poly_opvar2d([Bu_ode;Bu_pde]);

% Include contribution of PDE output r to outputs z and y
Cz_ode = Cz_ode + Dzr*Crx;
Cz_pde = Dzr*Crp;
Cz = poly_opvar2d([Cz_ode,Cz_pde]);

Cy_ode = Cy_ode + Dyr*Crx;
Cy_pde = Dyr*Crp;
Cy = poly_opvar2d([Cy_ode,Cy_pde]);

Dzw = poly_opvar2d(Dzw + Dzr*Drw);
Dyw = poly_opvar2d(Dyw + Dyr*Drw);
Dzu = poly_opvar2d(Dzu + Dzr*Dru);
Dyu = poly_opvar2d(Dyu + Dyr*Dru);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finally, get rid of terms below tolerance, and define the PIE structure

tol = 1e-14;
PIE.T = clean_opvar(T,tol);
PIE.Tw = clean_opvar(Tw,tol);
PIE.Tu = clean_opvar(Tu,tol);
PIE.A = clean_opvar(A,tol);
PIE.Bw = clean_opvar(Bw,tol);
PIE.Bu = clean_opvar(Bu,tol);
PIE.Cz = clean_opvar(Cz,tol);
PIE.Cy = clean_opvar(Cy,tol);
PIE.Dzw = clean_opvar(Dzw,tol);
PIE.Dzu = clean_opvar(Dzu,tol);
PIE.Dyw = clean_opvar(Dyw,tol);
PIE.Dyu = clean_opvar(Dyu,tol);

PIE.dom = dom;
PIE.vars = vars;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MATout = cell2mat_poly(CELLin)
% A simple (inefficient!) for-loop to concatenate polynomial objects in a
% 2D cell. Does not check if dimensions of the elements match.

% Use the more efficient built-in function if possible
try MATout = cell2mat(CELLin);
    return
catch
    [Nr,Nc] = size(CELLin);
    MATout = [];
    for i=1:Nr
        row_i = [];
        for j=1:Nc
            row_i = [row_i,CELLin{i,j}];
        end
        MATout = [MATout;row_i];
    end
end
end