%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a synthesis code for H-infty optimal controller 
% design (w/o control at the boundary) for a 4-PIE System defined
% by the 9 4-PI operator representation
% Eop \dot x = Aop x(t)+B1op w(t)+B2op u(t)
% z(t)=C1op x(t)+D11 w(t)+D12 u(t)
% y(t)=C2op x(t)+D21 w(t)+D22 u(t)
% u(t)=Kop x(t)
%
% Where for now, we assume $z(t)\in \R^{nz}$, $x(t) \in \R^{nx1} \times
% L_2^{nx2}$, w(t)\in \R^{nw}, u(t)\in \R^{nu}, y(t) \in \R^{ny}
% The domain is s \in [a,b]
%
% The following must be defined externally:
% Eop,Aop,B1op,C1op,B2op,C2op,D11,D12,D21,D22,nu,ny,nw,nz,nx1,nx2,n_order1,n_order2,a,b
% s and th must be pvars
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - executive_PIETOOLS_Hinf_control
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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS  - 7_26_2019
%

if ~exist('override1')
    override1=1; % default to no Psatz in the LF
end

if ~exist('override2')
    override2=0; % default to Psatz in the Derivative
end

% setup an SOS program

varlist = [s; theta];
pvar gamma;
prog = sosprogram(varlist,[gamma]);

% domain of spatial variables
% p-compact form of the interval
%g1 = (X(2)-s)*(s-X(1));

% hinf norm variable which needs to be minimized
%prog = sosdecvar(prog, gamma); %this sets gamma as decision var
prog = sossetobj(prog, gamma); %this minimizes gamma, comment for feasibility test

% setup the variables for lyapunov function candidate
disp('Parameterizing Positive function');
[prog, P1op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd1,options1);

if override1~=1
    [prog, P2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

[prog,Zop] = sos_opvar(prog,[nu nx1;0 nx2],X,s,theta,ddZ);

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nxb);  

%Assemble the big operator

Dop = [-gamma*eye(nz)   D11              (C1op*Pop+D12op*Zop)*(Eop');
        D11'                -gamma*eye(nw)  B1op';
        Eop*(C1op*Pop+D12op*Zop)'       B1op         (Aop*Pop+B2op*Zop)*(Eop')+Eop*(Aop*Pop+B2op*Zop)']; 

disp('Parameterize the derivative inequality');


%getting multiplier and kernel for derivative lyapunov function
disp('Parameterizing the negative operators');
[prog, De1op] = sos_posopvar(prog, [nw+nz+nx1, nx2],X,s,theta,dd2,options2);

if override2~=1
    [prog, De2op] = sos_posopvar(prog,[nw+nz+nx1, nx2],X,s,theta, dd3,options3);
    Deop=De1op+De2op;
else
    Deop=De1op;
end
% derivative negativity
% constraints
disp('Setting up the equality constraints');
prog = sosopeq(prog,Deop+Dop); %Dop=-Deop

% choosing a different solver if needed
% option.solver = 'mosek'; 

%solving the sos program
prog = sossolve(prog); 

disp('The an achievable closed-loop H-infty norm of the given system is upper bounded by:')
disp(double(sosgetsol(prog,gamma))); % check the Hinf norm, if the solved successfully