function PIEcl = closedLoopPIE(PIE, K, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% closedLoopPIE.m     PIETOOLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given a PIE of the form
% T \dot x +Tw w + Tu u= A x + B1 w + B2 u
% z = C1 x + D11 w + D12 u
% y = C2 x + D21 w + D22 u
% where the PI operators are stored as fields of PIE object
% this function returns the closed loop pie_struct object that corresponds to the 
% PIE 
% if type='controller' (default)
% (T + Tu*K) \dot x +Tw w = (A + B2*K) x + B1 w 
% z = (C1 + D12*K) x + D11 w 
% y = (C2 + D22*K) x + D21 w 
% 
% if type='observer'
% [T 0][\dot x   ] +[Tw] w +[Tu] u = [A     0     ][x   ] +[B1   ] w +[B2   ] u 
% [0 T][\dot xhat]  [0 ]    [0 ]     [-K*C2 A+K*C2][xhat] +[K*D21]    [K*D22]
% 
% z = [C1 0 ][x   ] +[D11] w +[D12] u
%     [0  C1][xhat]  [0  ]    [0  ]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE: a structure with fields T Tw Tu A B1 B2 C1 C2 D11 D12 D21 D22
% K : an opvar object with dimension equal to [B1.dim(:,2) T.dim(:,1)]
% type: optional string input; 'controller' (default) or 'observer' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUTS:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEcl: a pie_struct object with fields Tcl Tw Acl B1 C1cl D11 C2cl D21

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 5/30/2021
% added observer closed-loop construction, SS - 1/18/2024
% DJ, 12/07/2024: Make sure variables of output match those of input;
% DJ, 01/04/2025: Add support for 2D PIEs;

if nargin==2
    type = 'controller';
elseif nargin==3
    % do nothing
end

if ~isa(PIE,'pie_struct')
    error('First input to closedLoopPIE() must be a pie_struct object');
end
PIE = initialize(PIE);
dim = PIE.dim;                                                              % DJ, 01/04/2025
if dim>2
    error("Construction of PIEs in more than 2 spatial variables is currently not supported.")
end

if ~isa(K, 'opvar') && ~isa(K,'opvar2d')
    error('Second input to closedLoopPIE() must be an opvar or opvar2d class object');
end

if ~strcmp(type,'controller') && ~strcmp(type,'observer')
    error('Third input to closedLoopPIE() must be "controller" or "observer"');
end

fieldList = fieldnames(PIE);

requiredList = {'T'; 'Tw'; 'Tu'; 'A'; 'B1'; 'B2'; 'C1'; 'D11'; 'D12'; 'C2'; 'D21'; 'D22'};

if ~ismember(fieldList, requiredList)
    error('PIE object must have fields T Tw Tu A B1 B2 C1 C2 D11 D12 D21 D22');
end

if ~isa(PIE.T,'opvar') && ~isa(PIE.T,'opvar2d')
    error('PIE.T must be an opvar class object');
end
if ~isa(PIE.Tw,'opvar') && ~isa(PIE.Tw,'opvar2d')
    error('PIE.Tw must be an opvar class object');
end
if ~isa(PIE.Tu,'opvar') && ~isa(PIE.Tu,'opvar2d')
    error('PIE.Tu must be an opvar class object');
end
if ~isa(PIE.A,'opvar') && ~isa(PIE.A,'opvar2d')
    error('PIE.A must be an opvar class object');
end
if ~isa(PIE.B1,'opvar') && ~isa(PIE.B1,'opvar2d')
    error('PIE.B1 must be an opvar class object');
end
if ~isa(PIE.B2,'opvar') && ~isa(PIE.B2,'opvar2d')
    error('PIE.B2 must be an opvar class object');
end
if ~isa(PIE.C1,'opvar') && ~isa(PIE.C1,'opvar2d')
    error('PIE.C1 must be an opvar class object');
end
if ~isa(PIE.C2,'opvar') && ~isa(PIE.C2,'opvar2d')
    error('PIE.C2 must be an opvar class object');
end
if ~isa(PIE.D11,'opvar') && ~isa(PIE.D11,'opvar2d')
    error('PIE.D11 must be an opvar class object');
end
if ~isa(PIE.D12,'opvar') && ~isa(PIE.D12,'opvar2d')
    error('PIE.D12 must be an opvar class object');
end
if ~isa(PIE.D21,'opvar') && ~isa(PIE.D21,'opvar2d')
    error('PIE.D21 must be an opvar class object');
end
if ~isa(PIE.D22,'opvar') && ~isa(PIE.D22,'opvar2d')
    error('PIE.D22 must be an opvar class object');
end


% for convenience extract PI operators from PIE
T = PIE.T;      Tw = PIE.Tw;    Tu = PIE.Tu;
A = PIE.A;      B1 = PIE.B1;    B2 = PIE.B2;
C1 = PIE.C1;    D11 = PIE.D11;  D12 = PIE.D12;
C2 = PIE.C2;    D21 = PIE.D21;  D22 = PIE.D22;

% Build controller closed-loop system by default
if strcmp(type,'controller') 
    if ~all(K.dim(:,2)==Tu.dim(:,1))
        error('Input dimensions of controller gain should match size of the PIE state')
    elseif ~all(K.dim(:,1)==Tu.dim(:,2))
        error('Output dimensions of controller gain should match size of the actuator input')
    end

    pie_struct PIEcl;
    PIEcl.vars = PIE.vars;                                                  % DJ, 12/07/2024
    PIEcl.dom = PIE.dom;
    PIEcl.T = T + Tu*K;     PIEcl.Tw = Tw;
    PIEcl.A = A + B2*K;     PIEcl.B1 = B1;
    PIEcl.C1 = C1 + D12*K;  PIEcl.D11 = D11;
    PIEcl.C2 = C2 + D22*K;  PIEcl.D21 = D21;

    PIEcl = initialize(PIEcl);
    return
end

% % Otherwise, build observer closed-loop system
L = K;
if ~all(L.dim(:,2)==C2.dim(:,1))
    error('Input dimensions of observer gain should match size of the observed output')
elseif ~all(L.dim(:,1)==C2.dim(:,2))
    error('Output dimensions of observer gain should match size of the PIE state')
end

nx = T.dim(:,2); nw = Tw.dim(:,2); nu = Tu.dim(:,2);
nz = C1.dim(:,1);

if dim<=1
    opvar Tcl Twcl Tucl Acl B1cl B2cl C1cl C2cl D11cl D12cl D21cl D22cl;
else
    opvar2d Tcl Twcl Tucl Acl B1cl B2cl C1cl C2cl D11cl D12cl D21cl D22cl;
end
C2cl.dim = [zeros(length(nx),1),2*nx]; 
D21cl.dim = [zeros(length(nw),1),nw]; 
D22cl.dim = [zeros(length(nu),1),nu];
C1cl.dim = [2*nz,2*nx]; D11cl.dim = [2*nz,nw]; D12cl.dim = [2*nz,nu];
Tcl.dim = [2*nx,2*nx]; Twcl.dim = [2*nx,nw]; Tucl.dim = [2*nx,nu];
Acl.dim = [2*nx,2*nx]; B1cl.dim = [2*nx,nw]; B2cl.dim = [2*nx,nu];
Tcl.var1 = T.var1;  Twcl.var1 = T.var1;     Tucl.var1 = T.var1;             % DJ, 12/07/2024
Acl.var1 = T.var1;  B1cl.var1 = T.var1;     B2cl.var1 = T.var1;
C1cl.var1 = T.var1; D11cl.var1 = T.var1;    D12cl.var1 = T.var1;
C2cl.var1 = T.var1; D21cl.var1 = T.var1;    D22cl.var1 = T.var1;
Tcl.var2 = T.var2;  Twcl.var2 = T.var2;     Tucl.var2 = T.var2;
Acl.var2 = T.var2;  B1cl.var2 = T.var2;     B2cl.var2 = T.var2;
C1cl.var2 = T.var2; D11cl.var2 = T.var2;    D12cl.var2 = T.var2;
C2cl.var2 = T.var2; D21cl.var2 = T.var2;    D22cl.var2 = T.var2;
Tcl.I = T.I;        Twcl.I = T.I;           Tucl.I = T.I;             
Acl.I = T.I;        B1cl.I = T.I;           B2cl.I = T.I;
C1cl.I = T.I;       D11cl.I = T.I;          D12cl.I = T.I;
C2cl.I = T.I;       D21cl.I = T.I;          D22cl.I = T.I;

% construct parameters of z-output
if dim<=1                                                                   % DJ, 01/04/2025
    C1cl.P = [C1.P zeros(nz(1),nx(1)); zeros(nz(1),nx(1)) C1.P]; 
    C1cl.Q1 = [C1.Q1 zeros(nz(1),nx(2)); zeros(nz(1),nx(2)) C1.Q1];
    D11cl.P = [D11.P; zeros(nz(1),nw(1))]; 
    D12cl.P = [D12.P; zeros(nz(1),nu(1))]; 

    % construct dynamics parameters
    LC2 = L*C2;
    Acl.P = [A.P zeros(nx(1),nx(1)); -LC2.P   A.P+LC2.P]; 
    Acl.Q1 = [A.Q1 zeros(nx(1),nx(2)); -LC2.Q1   A.Q1+LC2.Q1];
    Acl.Q2 = [A.Q2 zeros(nx(2),nx(1)); -LC2.Q2   A.Q2+LC2.Q2];
    Acl.R.R0 = [A.R.R0 zeros(nx(2),nx(2)); -LC2.R.R0   A.R.R0+LC2.R.R0];
    Acl.R.R1 = [A.R.R1 zeros(nx(2),nx(2)); -LC2.R.R1   A.R.R1+LC2.R.R1];
    Acl.R.R2 = [A.R.R2 zeros(nx(2),nx(2)); -LC2.R.R2   A.R.R2+LC2.R.R2];
    LD21 = L*D21;
    B1cl.P = [B1.P; LD21.P]; 
    B1cl.Q2 = [B1.Q2; LD21.Q2];
    LD22 = L*D22;
    B2cl.P = [B2.P; LD22.P]; 
    B2cl.Q2 = [B2.Q2; LD22.Q2];
    
    Tcl.P = [T.P zeros(nx(1),nx(1));zeros(nx(1),nx(1))   T.P]; 
    Tcl.Q1 = [T.Q1 zeros(nx(1),nx(2)); zeros(nx(1),nx(2))   T.Q1];
    Tcl.Q2 = [T.Q2 zeros(nx(2),nx(1)); zeros(nx(2),nx(1))   T.Q2];
    Tcl.R.R0 = [T.R.R0 zeros(nx(2),nx(2));zeros(nx(2),nx(2))   T.R.R0];
    Tcl.R.R1 = [T.R.R1 zeros(nx(2),nx(2));zeros(nx(2),nx(2))   T.R.R1];
    Tcl.R.R2 = [T.R.R2 zeros(nx(2),nx(2));zeros(nx(2),nx(2))   T.R.R2];
    
    Twcl.P = [Tw.P;zeros(nx(1),nw(1))]; 
    Twcl.Q2 = [Tw.Q2; zeros(nx(2),nw(1))];
    
    Tucl.P = [Tu.P;zeros(nx(1),nu(1))]; 
    Tucl.Q2 = [Tu.Q2; zeros(nx(2),nu(1))];
else
    opvar2d Oxx Ozx Oxw Oxu Ozw Ozu;
    Oxx.dim = [nx,nx];  Oxx.var1 = T.var1;  Oxx.var2 = T.var2;  Oxx.I = T.I;
    Ozx.dim = [nz,nx];  Ozx.var1 = T.var1;  Ozx.var2 = T.var2;  Ozx.I = T.I;
    Oxw.dim = [nx,nw];  Oxw.var1 = T.var1;  Oxw.var2 = T.var2;  Oxw.I = T.I;
    Oxu.dim = [nx,nu];  Oxu.var1 = T.var1;  Oxu.var2 = T.var2;  Oxu.I = T.I;
    Ozw.dim = [nz,nw];  Ozw.var1 = T.var1;  Ozw.var2 = T.var2;  Ozw.I = T.I;
    Ozu.dim = [nz,nu];  Ozu.var1 = T.var1;  Ozu.var2 = T.var2;  Ozu.I = T.I;
    
    % Construct dynamics parameters
    LC2 = L*C2;                 LD21 = L*D21;           LD22 = L*D22;
    Tcl = [T,Oxx; Oxx,T];       Twcl = [Tw;Oxw];        Tucl = [Tu;Oxu];
    Acl = [A,Oxx; -LC2,A+LC2];  B1cl = [B1;LD21];       B2cl = [B2;LD22];

    % Construct output parameters
    C1cl = [C1,Ozx; Ozx,C1];    D11cl = [D11;Ozw];      D12cl = [D12;Ozu];
end


% assign values to pie_struct object container
pie_struct PIEcl;
PIEcl.vars = PIE.vars;                                                      % DJ, 12/07/2024
PIEcl.dom = PIE.dom;
PIEcl.T=Tcl;    PIEcl.Tw=Twcl;      PIEcl.Tu=Tucl;
PIEcl.A=Acl;    PIEcl.B1=B1cl;      PIEcl.B2=B2cl;
PIEcl.C1=C1cl;  PIEcl.D11=D11cl;    PIEcl.D12=D12cl;
PIEcl.C2=C2cl;  PIEcl.D21=D21cl;    PIEcl.D22=D22cl;

PIEcl = initialize(PIEcl);

end