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
% z = (C1 + D12*K)x + D11 w 
% y = (C2 + D22*K) x + D21 w 
% 
% if type='observer'
% [T 0; 0 T][\dot x; \dot xhat] +[Tw; 0] w+[Tu; 0] u = [A 0; -K*C2 A+K*C2][x; xhat] + [B1; K*D21] w +[B2; K*D22] u
% z = [C1 0;0 C1][x; xhat] + [D11; 0] w + [D12; 0] u
% y = [0 0][x; xhat] + [0] w + [0] u 
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

if nargin==2
    type = 'controller';
elseif nargin==3
    % do nothing
end

if ~isa(PIE,'pie_struct')
    error('First input to closedLoopPIE() must be a pie_struct object');
end

if ~isa(K, 'opvar')
    error('Second input to closedLoopPIE() must be an opvar class object');
end

if ~strcmp(type,'controller') && strcmp(type,'observer')
    error('Third input to closedLoopPIE() must be "controller" or "observer"');
end

fieldList = fieldnames(PIE);

requiredList = {'T'; 'Tw'; 'Tu'; 'A'; 'B1'; 'B2'; 'C1'; 'D11'; 'D12'; 'C2'; 'D21'; 'D22'};

if ~ismember(fieldList, requiredList)
    error('PIE object must have fields T Tw Tu A B1 B2 C1 C2 D11 D12 D21 D22');
end

if ~isa(PIE.T,'opvar')
    error('PIE.T must be an opvar class object');
end
if ~isa(PIE.Tw,'opvar')
    error('PIE.Tw must be an opvar class object');
end
if ~isa(PIE.Tu,'opvar')
    error('PIE.Tu must be an opvar class object');
end
if ~isa(PIE.A,'opvar')
    error('PIE.A must be an opvar class object');
end
if ~isa(PIE.B1,'opvar')
    error('PIE.B1 must be an opvar class object');
end
if ~isa(PIE.B2,'opvar')
    error('PIE.B2 must be an opvar class object');
end
if ~isa(PIE.C1,'opvar')
    error('PIE.C1 must be an opvar class object');
end
if ~isa(PIE.C2,'opvar')
    error('PIE.C2 must be an opvar class object');
end
if ~isa(PIE.D11,'opvar')
    error('PIE.D11 must be an opvar class object');
end
if ~isa(PIE.D12,'opvar')
    error('PIE.D12 must be an opvar class object');
end
if ~isa(PIE.D21,'opvar')
    error('PIE.D21 must be an opvar class object');
end
if ~isa(PIE.D22,'opvar')
    error('PIE.D22 must be an opvar class object');
end


% for convenience extract PI operators from PIE
T = PIE.T; Tw = PIE.Tw; Tu = PIE.Tu;
A = PIE.A; B1 = PIE.B1; B2 = PIE.B2;
C1 = PIE.C1; D11 = PIE.D11; D12 = PIE.D12;
C2 = PIE.C2; D21 = PIE.D21; D22 = PIE.D22;

if strcmp(type,'controller') % build controller closed-loop system
PIEcl.T = PIE.T+PIE.Tu*K;
PIEcl.Tw = PIE.Tw;
PIEcl.A = PIE.A+PIE.B2*K; 
PIEcl.B1 = PIE.B1;
PIEcl.C1 = PIE.C1+PIE.D12*K;
PIEcl.D11 = PIE.D11;
PIEcl.C2 = PIE.C2+PIE.D22*K;
PIEcl.D21 = PIE.D21;

opvar dummy;
dummy.I = PIE.T.I; dummy.dim = [PIE.T.dim(:,1),[0;0]];
dummy.var1 = PIE.T.var1; dummy.var2 = PIE.T.var2;

PIEcl.Tu = dummy;
PIEcl.B2 = dummy;
dummy.dim = [PIE.C1.dim(:,1),[0;0]];
PIEcl.D12 = dummy;
dummy.dim = [PIE.C2.dim(:,1),[0;0]];
PIEcl.D22 = dummy;
else % build observer closed-loop system
L = K;

nx = T.dim(:,2); nw = Tw.dim(:,2); nu = Tu.dim(:,2);
nz = C1.dim(:,1); ny = C2.dim(:,1);

opvar Tcl Twcl Tucl Acl B1cl B2cl C1cl C2cl D11cl D12cl D21cl D22cl;
C2cl.dim = [[1;0],2*nx]; D21cl.dim = [[1;0],nw]; D22cl.dim = [[1;0],nu];
C1cl.dim = [2*nz,2*nx]; D11cl.dim = [2*nz,nw]; D12cl.dim = [2*nz,nu];
Tcl.dim = [2*nx,2*nx]; Twcl.dim = [2*nx,nw]; Tucl.dim = [2*nx,nu];
Acl.dim = [2*nx,2*nx]; B1cl.dim = [2*nx,nw]; B2cl.dim = [2*nx,nu];

% construct parameters of z-output
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

% assign values to pie_struct object container
PIEcl.T=Tcl;PIEcl.Tw=Twcl;PIEcl.Tu=Tucl;
PIEcl.A=Acl;PIEcl.B1=B1cl;PIEcl.B2=B2cl;
PIEcl.C1=C1cl;PIEcl.D11=D11cl;PIEcl.D12=D12cl;
PIEcl.C2=C2cl;PIEcl.D21=D21cl;PIEcl.D22=D22cl;
end

PIEcl.dom = PIE.dom;
PIEcl.vars = PIE.vars;

PIEcl = pie_struct(PIEcl);
end