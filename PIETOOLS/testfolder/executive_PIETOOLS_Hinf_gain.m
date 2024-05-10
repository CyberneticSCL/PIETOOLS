%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes an H-infty gain analysis for a 4-PIE System defined
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

% setup an SOS program

varlist = [s; theta];
prog = sosprogram(varlist);


% hinf norm variable which needs to be minimized
pvar gamma;
prog = sosdecvar(prog, gamma); %this sets gamma as decision var
prog = sossetobj(prog, gamma); %this minimizes gamma, comment for feasibility test

% setup the variables for lyapunov function candidate
disp('Parameterizing Positive Lyapunov Operator');
%nx1=10
%nx2=10

[prog, P1op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd1,options1);

if override1~=1
    [prog, P2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end
%toc

%sospos_RL2RL_MP_7_24_19
% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nxb);  

%Assemble the big operator
% Pheq = [ -gamma*I  D'       B*PH
%          D         -gamma*I C
%          H*PB      C*       A*Ph+H*PA]

disp('Constructing the derivative inequality');

Dop = [-gamma*eye(nw)   D11op'              (B1op')*(Pop*Eop);
        D11op                -gamma*eye(nz)  C1op;
        Eop'*(Pop*B1op)       C1op'             Eop'*(Pop*Aop)+Aop'*(Pop*Eop)]; 


%getting multiplier and kernel for derivative lyapunov function
disp('Parameterizing the negative operators');


% derivative negativity
prog = sosopineq(prog,-Dop,opts);
% choosing a different solver if needed
% option.solver = 'mosek'; 

%solving the sos program
prog = sossolve(prog); 

disp('The H-infty norm of the given system is upper bounded by:')
disp(double(sosgetsol(prog,gamma))); % check the Hinf norm, if the solved successfully
