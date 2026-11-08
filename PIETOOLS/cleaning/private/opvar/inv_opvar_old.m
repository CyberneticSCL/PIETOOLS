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
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 12/07/2024: Set vars of Pinv equal to vars of P;
% DJ, 03/25/2025: Bugfix in case degmats of parameters are empty;
% DJ, 04/07/2025: Add support for non self-adjoint operators;
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

I = P.I; s = P.var1; s_dum = P.var2;
l=I(1);     u=I(2);

if P.dim(2,2) == 0
    opvar Pinv;
    Pinv.I = I;
    Pinv.var1 = s;  Pinv.var2 = s_dum;                                      % DJ, 12/07/2024
    Pinv.P = inv(double(P.P)); Pinv.Q1=[]; Pinv.Q2 = []; Pinv.R.R0= []; 
    return
end

% Extract the parameters defining the operator
NP = double(P.P);
NQ1 = polynomial(dpvar2poly(dpvar(P.Q1)));   % DJ, 12/30/2021
NQ2 = polynomial(dpvar2poly(dpvar(P.Q2)));                                  % DJ, 04/07/2025
Sa = polynomial(dpvar2poly(dpvar(P.R.R0)));
Ra = polynomial(dpvar2poly(dpvar(P.R.R1)));
NR = Ra; NS = Sa;

% Check that the operator is square
[n1_dim,n2_dim] = size(NP);
[m1_dim,m2_dim] = size(NS);
if n1_dim~=n2_dim || m1_dim~=m2_dim
    error("Inversion of 'opvar' objects with different row and column dimensions is currently not supported.")
end

% if ~iscell(I)
%     error(['I must be a cell array of intervals'])
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The first step is to construct the separable representation of the
% Polynomials NQ1, NQ2, and NR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For NQ1, NQ2 and NR, we need to convert to the forms
%   NQ1 = H1*Z(s), NQ2 = Z(s)^T*H2 and NR= Z(s)^T Gam Z(theta)
% where we need to be consistent about our Z

dQ1 = NQ1.degmat;   dQ2 = NQ2.degmat;
dmQ = max([max(max(dQ1)),max(max(dQ2))]);
dR = NR.degmat;
dmR = max(max(dR));
dm = max([dmQ dmR]);
if isempty(dm) || dm==0                                                     % DJ, 03/25/2025
    Zs = polynomial(1);
    dm = 0;
else
    Zs = polynomial(eye(dm+1),(0:dm)',s.varname,[dm+1 1]); % common vector of monomials
end
bigZsL = [];
for i=1:m1_dim
    bigZsL = blkdiag(bigZsL,Zs);
end
bigZsR = [];
for i=1:m2_dim
    bigZsR = blkdiag(bigZsR,Zs);
end

% Now find H1 s.t. NQ1 = H1*bigZs
nZ=length(Zs);
H1 = zeros(n1_dim,nZ*m2_dim);
for i=1:n1_dim
    for j=1:m2_dim
        [CQ1ij,~,~] = poly2basis(NQ1(i,j),Zs); % examine each element of NQ1
        H1(i,(nZ*(j-1)+1):(nZ*j)) = CQ1ij.';
    end
end
% H1*bigZsR-NQ1  uncomment to verify the representation (should be 0)

% Now find H2 s.t. NQ2 = bigZs'*H2
H2 = zeros(nZ*m1_dim,n2_dim);
for i=1:m1_dim
    for j=1:n2_dim
        [CQ2ij,~,~] = poly2basis(NQ2(i,j),Zs); % examine each element of NQ2
        H2((nZ*(i-1)+1):(nZ*i),j) = CQ2ij;
    end
end
% bigZsL'*H2-NQ2  uncomment to verify the representation (should be 0)

% Now find Gam s.t. NR(s,theta) = bigZCth(s)'*Gam*bigZCth(theta)
Gam = zeros(size(bigZsL,1),size(bigZsR,1));
for i=1:m1_dim
    for j=1:m2_dim % address the decomposition of each element separately and place in the larger matrix Q
        dmij=NR(i,j).degmat;
        if isempty(dmij)                                                    % DJ, 03/25/2025
            dmij = [0,0];
        elseif size(dmij,2)==1
            if strcmp(NR(i,j).varname{1},s.varname)
                dmij = [dmij,zeros(size(dmij))];
            else
                dmij = [zeros(size(dmij)),dmij];
            end
        elseif size(dmij,2)==2
            if strcmp(NR(i,j).varname{1},s_dum.varname{1})
                dmij = fliplr(dmij);
            end
        else
            error("Parameters 'R.R1' and 'R.R2' can depend on at most two variables.")
        end
        cfij = NR(i,j).coeff;
        CN = zeros(dm+1,dm+1);
        for k=1:length(cfij) % take each coefficient and put it in its proper place in CN, which is then assembled into Gam
            CN(full(dmij(k,1))+1,full(dmij(k,2))+1) = cfij(k);
        end
        Gam(((dm+1)*(i-1)+1):((dm+1)*i),((dm+1)*(j-1)+1):((dm+1)*j)) = CN;
    end
end
%bigZsL'*Gam*subs(bigZsR,s,s_dum)-NR % this should be zero

% At this point, we can represent
%   P = [0,0 ] + [I,0 ] [NP,H1 ] [I,0 ]
%       [0,Sa]   [0,ZL] [H2,Gam] [0,ZR],
% where (ZL*v)(s)=bigZsL(s)'*v, and (ZR*u)=int_{l}^{u}bigZsR(th)*u(th)dth.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now that we have the right representation, we can construct the inverse.
% We assume it has the form
%   Pinv = [0,0 ] + [I,0   ] [NP_h,H1_h ] [I,0   ]
%          [0,Sh]   [0,ZL_h] [H2_h,Gam_h] [0,ZR_h],
% where Sh(s) = Sa(s)^{-1}, and
%   (ZL_h*v)(s)=Sh(s)*bigZsL(s)'*v,   
%   (ZR_h*u)=int_{l}^{u}bigZsR(th)*Sh(th)*u(th)dth,
% so that
%   Pinv*P = [0,0] + [0,0   ] [NP,H1 ] [I,0 ]
%            [0,I]   [0,ZL_h] [H2,Gam] [0,ZR]
%               + [I,0   ] [NP_h,H1_h ] [0,0 ]
%                 [0,ZL_h] [H2_h,Gam_h] [0,ZR]
%                   + [I,0   ] [NP_h,H1_h ] [I,0] [NP,H1 ] [I,0 ]
%                     [0,ZL_h] [H2_h,Gam_h] [0,K] [H2,Gam] [0,ZR],
% where K = ZR_h*ZL = int_{l}^{u}bigZsR(th)*Sh(th)*bigZsL(th)'dth.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The first step is to find a polynomial approximation to Sa^{-1},
% and compute K = int_{l}^{u}bigZsR(th)*Sh(th)*bigZsL(th)'dth.
if isdouble(NS)                                                             % DJ, 03/25/2025
    Sh = polynomial(pinv(double(NS)));
    K = double(int(bigZsR*Sh*bigZsL.',s,l,u));
else  
    % If NS is polynomial, it's inverse likely isn't. To approximate,
    % we compute an inverse at different grid points, and interpolate
    % by a polynomial
    N=100;      % number of data points to fit at
    orderapp=6; % degree of the polynomial approximation
    % It is recommended to use at least a 4th order approximation
    
    s_step = (u-l)/N; % dustance between grid points
    ii=0;
    
    for sval=(l:s_step:u)
        ii=ii+1;
        NShtemp(:,:,ii)=inv(double(subs(NS,s,sval))); % Calculates the value of the inverse of S at every point in the interval
    end
    
    % The following fits a polynomial to every element of S^{-1}
    Sinv = polynomial(zeros(m1_dim,m1_dim));
    for i=1:m1_dim
        for j=1:m1_dim
            Data1=squeeze(NShtemp(i,j,:))';
            tempCoeffs =polyfit((l:s_step:u),Data1,orderapp); % uses matlab internal polynomial representation 
    %         [temp,S,mu]=polyfit([l:interval:u],Data1,orderapp); % uses matlab internal polynomial representation
    %        S.normr
    %         syms th2   %Converts to a symbolic variable representation
    %         Sinv_temp=s2p(poly2sym(temp,th2)); % this approach works, but seems fragile. 
           
            Sinv(i,j)=s.^(orderapp:-1:0)*tempCoeffs';
            clear Data1 temp Sinv_temp
        end
    end
    % This is the polynomial that is returned. It is not used to calculate the
    % following integral, however.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Now we compute the matrix K (does not use the approximate S^{-1})
    % This is an anonymous function which is the integrand 
    int_fun = @(sss) double(subs(bigZsL,s,sss))*inv(double(subs(NS,s,sss)))*double(subs(bigZsL,s,sss)).';
    % Uses Matlab internal numerical integration routine.
    K=integral(int_fun,l,u,'ArrayValued',true);
    
    % We're ready to compute the matrices to return
    Sh=Sinv;
end

% Next, we need to find NP_h, H1_h, H2_h, Gam_h such that
%   [0,0   ] + [0,H1_h ] + [NP_h,H1_h ] [NP,H1     ] = [I,0]
%   [H2,Gam]   [0,Gam_h]   [H2_h,Gam_h] [K*H2,K*Gam]   [0,0]
% Thus
%   [0,H1_h ] + [NP_h,H1_h ] [NP,H1     ] = [I,0] - [0,0 ]
%   [0,Gam_h]   [H2_h,Gam_h] [K*H2,K*Gam]   [0,0]   [H2,Gam]
% We can solve this for NP_h, H1_h, H2_h, and Gam_h
AMat = blkdiag(eye(n2_dim,n1_dim),zeros(size(Gam))) - [zeros(n2_dim,n2_dim+size(Gam,2)); H2,Gam];
BMat = blkdiag(zeros(n1_dim,n2_dim),eye(size(Gam,2))) + [NP,H1; K*H2,K*Gam];
bigMat_h = AMat/BMat;                                                       % DJ, 04/07/2025

NP_h = bigMat_h(1:n2_dim,1:n1_dim);   
H1_h = bigMat_h(1:n2_dim,n1_dim+1:end);
H2_h = bigMat_h(n2_dim+1:end,1:n1_dim);
Gam_h = bigMat_h(n2_dim+1:end,n1_dim+1:end);

opvar Pinv;
Pinv.I = I;
Pinv.var1 = s;  Pinv.var2 = s_dum;                                          % DJ, 12/07/2024
Pinv.P = NP_h; 
Pinv.Q1 = H1_h*bigZsR*Sh;                                 % Same degree as S^{-1}
Pinv.Q2 = Sh*bigZsL'*H2_h; 
Pinv.R.R0= Sh; 
Pinv.R.R1 = (Sh*bigZsL')*Gam_h*subs(bigZsR*Sh,s,s_dum);   % Twice degree as S^{-1}
Pinv.R.R2 = Pinv.R.R1;

end