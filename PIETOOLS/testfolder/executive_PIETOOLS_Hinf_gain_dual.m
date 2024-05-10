%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes an alternative `dual' H-infty gain analysis for a 4-PIE System defined
% by the 9 4-PI operator representation
% Eop \dot x=Aop x(t)+B1op w(t)
% z(t)=C1op x(t)+D11 w(t)
%
% Where for now, we assume $z(t)\in \R^{nz}$, $x(t) \in \R^{nx1} \times
% L_2^{nx2}$, w(t)\in \R^{nw}
% The domain is s \in [a,b]
%
% The following must be defined externally:
% Eop,Aop,B1op,C1op,D11,nw,nz,nx1,nx2,n_order1,n_order2,a,b
% s and th must be pvars
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - executive_PIETOOLS_Hinf_gain_dual
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



% setup an SOS program

varlist = [s; theta];
prog = sosprogram(varlist);

% domain of spatial variables

% p-compact form of the interval
%g1 = (X(2)-s)*(s-X(1));

% hinf norm variable which needs to be minimized
pvar gamma;
prog = sosdecvar(prog, gamma); %this sets gamma as decision var
prog = sossetobj(prog, gamma); %this minimizes gamma, comment for feasibility test

% setup the variables for lyapunov function candidate
disp('Parameterizing Positive function');
%nx1=10
%nx2=10


[prog, P1op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd1,options1);

if override1~=1
    [prog, P2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nxb);  

%Assemble the big operator. The LMI is
% Pheq = [ -gamma*I  D'       B*PH
%          D         -gamma*I C
%          H*PB      C*       A*PH+H*PA]

Dop = [-gamma*eye(nz)   D11op              C1op*Pop*(Eop');
        D11op'                -gamma*eye(nw)  B1op';
        Eop*Pop*(C1op')       B1op         Aop*Pop*(Eop')+Eop*Pop*(Aop')]; 

    
disp('Parameterize the derivative inequality');


%getting multiplier and kernel for derivative lyapunov function
disp('Parameterizing the negative operators');


[prog, Deop] = sos_posopvar(prog, [nw+nz+nx1, nx2],X,s,theta,dd2,options2);
[prog, De2op] = sos_posopvar(prog,[nw+nz+nx1, nx2],X,s,theta, dd3,options3);
% derivative negativity
% constraints
disp('Setting up the equality constraints');
prog = sosopeq(prog,Deop+De2op+Dop); %Dop=-Deop-De2op

% choosing a different solver if needed
% option.solver = 'mosek'; 

%solving the sos program
prog = sossolve(prog); 

disp('The H-infty norm of the given system is upper bounded by:')
disp(double(sosgetsol(prog,gamma))); % check the Hinf norm, if the solved successfully
