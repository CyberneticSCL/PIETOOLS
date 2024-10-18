function [Pinv] = inv_opvar_old(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv_opvar(P) computes the inverse operator of the
% operator
%
% P[ x1    ](s)=[ P.P*x1 +int_I(1)^I(2) P.Q(th)x_2(th)d th                    ]
%  [ x2(s) ]    [ P.Q(s)^T*x1 +P.R.R0(s)x_2(s) + int_I(1)^I(2) P.R.R1(s,th)x_2(th)dth ]
%
% INPUT
%   P: positive definite opvar to invert
%   I = [l1 u1] interval of integration
%
% OUTPUT 
%   Pinv: inverse opvar object. Inverse opvar is a numerical inversion and
%   should be used with care and reservations.
% The inverse operator has the form
% P^-1[ x1    ](s)=[ Ph*x1 +int_I(1)^I(2) Qh(th)x_2(th)d th                    ]
%     [ x2(s) ]    [ Qh(s)^T*x1 +tau*Sh(s)x_2(s) + tau* int_I(1)^I(2) Rh(s,th)x_2(th)dth ]
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv_opvar
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP - 7_1_2020
% Assure parameters are dpvar before converting to poly, DJ - 12/30/2021
%       <-- Not really an optimal fix...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(P,'opvar')||isa(P,'dopvar')
try 
    err = double(P.R.R1-P.R.R2);
    if err~=0
        error('Inverse for opvar can only be calculated for opvar with R1=R2');
    end
catch ME
    error('Inverse for opvar can only be calculated for opvar with R1=R2');
end
end

I = P.I; s = P.var1; theta = P.var2;
l=I(1);u=I(2);
tau=u-l;
r=tau;

if P.dim(2,2) == 0
    opvar Pinv;
    Pinv.I = I;
    Pinv.P = inv(double(P.P)); Pinv.Q1=[]; Pinv.Q2 = []; Pinv.R.R0= []; 
    return
end

NP = polynomial(P.P);
NQ = polynomial(dpvar2poly(dpvar(P.Q1)));   % DJ, 12/30/2021
Sa = polynomial(dpvar2poly(dpvar(P.R.R0)));
Ra = polynomial(dpvar2poly(dpvar(P.R.R1)));
NR = Ra; NS = Sa;




n_dim=size(NP,1);
m_dim=size(NS,1);

% if ~iscell(I)
%     error(['I must be a cell array of intervals'])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step is to construct the separable representation of the
% Polynomials Q and R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For ZQ and NR, we need to convert to the forms
% NQ=H Z(s) and NR= Z(s)^T Gam Z(theta) as per paper
% or for us this is Q=HZ(theta) and NR= Z(theta)^T Gam Z(ksi)
% however, we need to be consistent about our Z

dQ=NQ.degmat;
dmQ=max(max(dQ));
dR=NR.degmat;
dmR=max(max(dR));
dm=max([dmQ dmR]);
ZCth=polynomial(eye(dm+1),[0:dm]',{'s'},[dm+1 1]); % common vector of monomials
bigZCth=[];
for i=1:m_dim
    bigZCth=blkdiag(bigZCth,ZCth);
end

% Now use bigZCth to find H!
nZ=length(ZCth);
bigC = zeros(n_dim,nZ*m_dim);
for i=1:n_dim
    for j=1:m_dim
        [CQij,Ztemp,etemp] = poly2basis(NQ(i,j),ZCth); % examine each element of NQ
        bigC(i,(nZ*(j-1)+1):(nZ*j))=CQij.';
    end
end
H=bigC;


% H*bigZCth-NQ  uncomment to verify the representation (should be 0)

% Now we deal with NR

for i=1:m_dim
    for j=1:m_dim % address the decomposition of each element separately and place in the larger matrix Q
        dmij=NR(i,j).degmat;
        cfij=NR(i,j).coeff;
        CN=zeros(dm+1,dm+1);
        for k=1:length(cfij) % take each coefficient and put it in its proper place in CN, which is then assembled into CbigN
            CN(full(dmij(k,1))+1,full(dmij(k,2))+1)=cfij(k);
        end
        Gam(((dm+1)*(i-1)+1):((dm+1)*i),((dm+1)*(j-1)+1):((dm+1)*j))=CN;
    end
end
%bigZCksi=subs(bigZCth,th,ksi);
%bigZCksi.'*Gam*bigZCth-NR % this should be zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we have the right representation, we can construct the inverse
% using Keqin's formula.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The first step is to find a polynomial approximation to S^{-1}
if isa(NS,'double')
    Sh =inv(NS);
else  
N=100; orderapp=6; % specifies the degree of the polynomial approximation and the number of data points to fit at.
% It is recommended to use at least a 4th order approximation

interval = tau/N; % specifies the degree of the polynomial approximation and the number of data points to fit at.
% It is recommended to use at least a 4th order approximation
ii=0;


th = s; ksi = theta;

for ss=[l:interval:u]
    ii=ii+1;
    NShtemp(:,:,ii)=inv(double(subs(NS,th,ss))); % Calculates the value of the inverse of S at every point in the interval
end

% The following fits a polynomial to every element of S^{-1}
for i=1:m_dim
    for j=1:m_dim
        Data1=squeeze(NShtemp(i,j,:))';
        tempCoeffs =polyfit([l:interval:u],Data1,orderapp); % uses matlab internal polynomial representation 
%         [temp,S,mu]=polyfit([l:interval:u],Data1,orderapp); % uses matlab internal polynomial representation
%        S.normr
%         syms th2   %Converts to a symbolic variable representation
%         Sinv_temp=s2p(poly2sym(temp,th2)); % this approach works, but seems fragile. 
       
        Sinv(i,j)=th.^(orderapp:-1:0)*tempCoeffs';
        clear Data1 temp Sinv_temp
    end
end
% This is the polynomial that is returned. It is not used to calculate the
% following integral, however.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Now we compute the matrix K (does not use the approximate S^{-1})
% This is an anonymous function which is the integrand 
int_fun = @(sss) double(subs(bigZCth,th,sss))*inv(double(subs(NS,th,sss)))*double(subs(bigZCth,th,sss)).';
% Uses Matlab internal numerical integration routine.
K=integral(int_fun,l,u,'ArrayValued',true);

% Now ready to compute the matrices to return

Sh=Sinv;
end
Z_h_th=bigZCth*Sinv; % This inherits its degree from S^{-1}

P_inv=inv(double(NP));

ndtemp=size(K,1);

T=inv(eye(ndtemp)+K*Gam-r*K*H.'*P_inv*H); % A Matrix
H_h=-P_inv*H*T; % A Matrix

%ndtemp2=size(P_inv,1);

Ph=(eye(n_dim) + r*P_inv*H*T*K*H.')*P_inv;
Qh = H_h*Z_h_th; % Same degree as S^{-1}

Gam_h=(r*T.'*H.'*P_inv*H-Gam)*inv(eye(ndtemp)+K*Gam); % A Matrix

Z_h_ksi=subs(Z_h_th,th,ksi);   
Rh = Z_h_th.'*Gam_h*Z_h_ksi;  % Twice the degree of S^{-1}

opvar Pinv;
Pinv.I = I;
Pinv.P = Ph; Pinv.Q1=Qh; Pinv.Q2 = Qh'; Pinv.R.R0= Sh; Pinv.R.R1 = Rh; Pinv.R.R2 = Rh;
end




