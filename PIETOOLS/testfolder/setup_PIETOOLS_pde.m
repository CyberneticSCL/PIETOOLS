if ~exist('a','var')||~exist('b','var')
    error('Domain not defined');
end
X = [a,b];

if ~exist('nu','var')
    warning('Number of inputs not defined. Defaulting to zero');
    nu=0;
end
if ~exist('nw','var')
    warning('Number of disturbances not defined. Defaulting to zero');
    nw=0;
end
if ~exist('ny','var')
    warning('Number of outputs not defined. Defaulting to zero');
    ny=0;
end
if ~exist('nz','var')
    warning('Number of regulated outputs not defined. Defaulting to zero');
    nz=0;
end

if ~exist('A','var') && ~exist('E0','var')
    disp('Number of ODE states Defaulting to zero');
    no = 0;
elseif ~exist('E0','var')
    no = size(A,1);
else 
    no = size(E0,1);
end

if ~exist('A','var')
    disp('A not defined. Defaulting to zero');
    A = zeros(no);
end

if exist('A0','var')
    np = size(A0,1);
elseif exist('A1','var')
    np = size(A1,1);
elseif exist('A2','var')
    np = size(A2,1);
end

if ~exist('A0','var')
    disp('A0 not defined. Defaulting to zero');
    n1 = 0;
    A0=zeros(np);
end

if ~exist('A1','var')
    disp('A1 not defined. Defaulting to zero');
    n2 = 0;
    A1=zeros(np,n2+n3);
end

if ~exist('A2','var')
    disp('A2 not defined. Defaulting to zero');
    n3 = 0;
    A2=zeros(np,n3);
end

if ~exist('E','var')
    disp('E not defined. Defaulting to zero');
    E = zeros(np,no);
elseif size(E,1)~=np
    disp('E has wrong dimensions. Resetting to zero');
    E = zeros(np,no);
end

if ~exist('E0','var')
    disp('E0 not defiend. Defaulting to zero');
    E0 = zeros(no,2*n2+4*n3);
end

if ~exist('Ea','var')
    Ea=zeros(no,np);
    disp('Ea undefined, defaulting to 0');
end
if ~exist('Eb','var')
    Eb=zeros(no,n2+n3);
    disp('Eb undefined, defaulting to 0');
end


if ~exist('B','var')
    error('Boundary conditions not defined.');
end

if ~exist('Bx0','var')
    disp('Bx0 not defined. Defaulting to zero');
    Bx0 = zeros(n2+2*n3,no);
end

if ~exist('Bw','var')
    disp('Bw not defined. Defaulting to zero');
    Bw = zeros(n2+2*n3,nw);
end

if ~exist('Bu','var')
    disp('Bu not defined. Defaulting to zero');
    Bu = zeros(n2+2*n3,nu);
end

if ~exist('B11','var')
    B11=zeros(no,nw);
    disp('B11 undefined, defaulting to 0')
end

if ~exist('B12','var')
    B12=zeros(no,nu);
    disp('B12 undefined, defaulting to 0')
end

if ~exist('B21','var')
    B21=zeros(np,nw);
    disp('B21 undefined, defaulting to 0')
end

if ~exist('B22','var')
    B22=zeros(np,nu);
    disp('B22 undefined, defaulting to 0')
end

if ~exist('C1','var')
    C1=zeros(nz,no);
    disp('C1 undefined, defaulting to 0')
end
if ~exist('C2','var')
    C2=zeros(ny,no);
    disp('C2 undefined, defaulting to 0')
end
if ~exist('C10','var')
    C10=zeros(nz,2*n2+4*n3);
    disp('C10 undefined, defaulting to 0')
end
if ~exist('C20','var')
    C20=zeros(ny,2*n2+4*n3);
    disp('C20 undefined, defaulting to 0')
end
if ~exist('Ca1','var')
    Ca1=zeros(nz,np);
    disp('Ca1 undefined, defaulting to 0')
end
if ~exist('Ca2','var')
    Ca2=zeros(ny,np);
    disp('Ca2 undefined, defaulting to 0')
end
if ~exist('Cb1','var')
    Cb1=zeros(nz,n2+n3);
    disp('Cb1 undefined, defaulting to 0')
end
if ~exist('Cb2','var')
    Cb2=zeros(ny,n2+n3);
    disp('Cb2 undefined, defaulting to 0')
end


if ~exist('D11','var')
    D11=zeros(nz,nw);
    disp('D11 undefined, defaulting to 0')
end
if ~exist('D12','var')
    D12=zeros(nz,nu);
    disp('D12 undefined, defaulting to 0')
end
if ~exist('D21','var')
    D21=zeros(ny,nw);
    disp('D21 undefined, defaulting to 0')
end
if ~exist('D22','var')
    D22=zeros(ny,nu);
    disp('D22 undefined, defaulting to 0')
end



T = [eye(n2) zeros(n2,n3) zeros(n2,n3);
     eye(n2) zeros(n2,n3) zeros(n2,n3);
     zeros(n3,n2) eye(n3) zeros(n3);
     zeros(n3,n2) eye(n3) (b-a)*eye(n3);
     zeros(n3,n2) zeros(n3) eye(n3);
     zeros(n3,n2) zeros(n3) eye(n3)];
Q = [zeros(n2,n1) zeros(n2) zeros(n2,n3);
     zeros(n2,n1) eye(n2) zeros(n2,n3);
     zeros(n3,n1) zeros(n3,n2) zeros(n3);
     zeros(n3,n1) zeros(n3,n2) (b-theta)*eye(n3);
     zeros(n3,n1) zeros(n3,n2) zeros(n3);
     zeros(n3,n1) zeros(n3,n2) eye(n3)];
K = [zeros(n1,n2) zeros(n1,n3) zeros(n1,n3);
     eye(n2) zeros(n2,n3) zeros(n2,n3);
     zeros(n3,n2) eye(n3) (s-a)*eye(n3)];
L0 = [eye(n1) zeros(n1,n2) zeros(n1,n3);
     zeros(n2,n1) zeros(n2) zeros(n2,n3);
     zeros(n3,n1) zeros(n3,n2) zeros(n3)];
L1 = [zeros(n1) zeros(n1,n2) zeros(n1,n3);
     zeros(n2,n1) eye(n2) zeros(n2,n3);
     zeros(n3,n1) zeros(n3,n2) (s-theta)*eye(n3)];
V = [zeros(n2,n2) zeros(n2,n3) zeros(n2,n3);
     zeros(n3,n2) zeros(n3,n3) eye(n3)];
F0 = [zeros(n2,n1) eye(n2) zeros(n2,n3);
      zeros(n3,n1) zeros(n3,n2) zeros(n3)];
F1 = [zeros(n2,n1) zeros(n2) zeros(n2,n3);
      zeros(n3,n1) zeros(n3,n2) eye(n3)];


%--------------------------------------------------------------------------
% converts all states to fundamental state and assembles A operator
opvar H0 H1 Hbf;
disp('Converting all states to fundamental state');
% H0 maps [z x1 x2s x3ss] to [z x1 x2 x3] 
H0.P = eye(no); H0.Q1 = zeros(no,np); H0.Q2 = K*inv(B*T)*Bx0;
H0.R.R0 = L0;
H0.R.R1 = L1 -K*inv(B*T)*B*Q; H0.R.R2 = -K*inv(B*T)*B*Q;
H0.dim = [no no; np np];
H0.I = [a b];
H0.var1 = s;
H0.var2 = theta;

% H1 maps [z x1 x2s x3ss] to [x2s x3s]
H1.Q2 = V*inv(B*T)*Bx0;
H1.R.R0 = F0; 
H1.R.R1 = F1 - V*inv(B*T)*B*Q; H1.R.R2 = - V*inv(B*T)*B*Q;
H1.dim = [0 no; n2+n3 np];
H1.I = [a b];
H1.var1 = s;
H1.var2 = theta;

% Hbf maps [z x1 x2s x3ss] to [x2(0) x2(L) x3(0) x3(L) x3s(0) x3s(L)]
Hbf.P = T*inv(B*T)*Bx0;
Hbf.Q1 = var_swap(-T*inv(B*T)*B*Q+Q, s, theta);
Hbf.dim = [2*n2+4*n3 no; 0 np];
Hbf.I = [a b];
Hbf.var1 = s;
Hbf.var2 = theta;


% convert operator E to fundamental state operator Ef
opvar EAop EBop;
EAop.Q1 = Ea; 
EAop.dim = [no no;0 np]; %?????????? 
EAop.I = [a b]; EAop.var1 = s; EAop.var2 = theta;
EBop.Q1 = Eb;
EBop.dim = [no 0;0 n2+n3]; 
EBop.I = [a b]; EBop.var1 = s; EBop.var2 = theta;

Af1 = E0*Hbf+EAop*H0+EBop*H1;
Af1.P = Af1.P+A;

opvar A2p;
A2p.R.R0 = [zeros(n1) zeros(n2) A2]; A2p.dim = [0 no; np np];
opvar Ip;
Ip.R.R0 = eye(np); Ip.dim = [0 no; np np];  
Ip.var1 = s; Ip.var2 = theta; Ip.I = [a b]; 

Af2 = A0*Ip*H0+A1*H1+A2p;
Af2.Q2 = Af2.Q2+E;
% Af = [A     E1*H3+EA*H0+EB*H2;
%       E     A0*H0+A1*H2+[zeros(np,n1) zeros(np,n2) A2]];
Af = [Af1; Af2];
clear Af1 Af2 EAop EBop A2p;

%--------------------------------------------------------------------------
opvar Bf1 Bf2 Df11 Df12 Df21 Df22;

%Assemble B operator
Bf1.P = B11;
Bf1.Q2 = B21;
Bf1.dim = [no nw; np 0];
Bf1.I = [a b];Bf1.var1 = s;Bf1.var2 = theta;

Bf2.P = B12;
Bf2.Q2 = B22;
Bf2.dim = [no nu; np 0];
Bf2.I = [a b];Bf2.var1 = s;Bf2.var2 = theta;


opvar CA1 CB1 CA2 CB2;
CA1.Q1 = Ca1; CB1.Q1 = Cb1;
CA1.dim = [nz no;0 np];
CA1.I = [a b];CA1.var1 = s;CA1.var2 = theta;
CB1.dim = [nz 0;0 n2+n3];
CB1.I = [a b];CB1.var1 = s;CB1.var2 = theta;
CA2.Q1 = Ca2; CB2.Q1 = Cb2;
CA2.dim = [ny no;0 np];
CA2.I = [a b];CA2.var1 = s;CA2.var2 = theta;
CB2.dim = [ny 0;0 n2+n3];
CB2.I = [a b];CB2.var1 = s;CB2.var2 = theta;


%Assemble C operators
Cf1 = C10*Hbf + CA1*H0 + CB1*H1; 
Cf1.P = Cf1.P + C1;
Cf2 = C20*Hbf + CA2*H0 + CB2*H1; 
Cf2.P = Cf2.P + C2;

clear CA1 CB1 CA2 CB2;

%Assemble D operators
Df11.P = D11; Df11.dim = [nz nw; 0 0]; Df11.I = [a b];
Df12.P = D12; Df12.dim = [nz nu; 0 0]; Df12.I = [a b];
Df21.P = D21; Df21.dim = [ny nw; 0 0]; Df21.I = [a b];
Df22.P = D22; Df22.dim = [ny nu; 0 0]; Df22.I = [a b];

Eop = H0;
Aop = Af;
C1op = Cf1;
C2op = Cf2;
B1op = Bf1;
B2op = Bf2;
D12op = Df12;
D11op = Df11;
D21op = Df21;
D22op = Df22;
nx1 = no; nx = no; nx2 = np; nxb = np;