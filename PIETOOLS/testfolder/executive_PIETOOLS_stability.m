%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a stability analysis for a 4-PIE System defined
% by the 4-PI operator representation
% Eop \dot x=Aop x(t)
%
% Where for now, we assume $x(t) \in \R^{nx1} \times L_2^{nx2}$, 
% The domain is s \in [a,b]
%
% The following must be defined externally:
% Eop,Aop,nx1,nx2,dd1,dd2,dd3,options1,options2,options,a,b
% s and th must be pvars

% setup an SOS program

varlist = [s; theta];
prog = sosprogram(varlist);



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
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  

%Assemble the big operator
% Pheq = [ -gamma*I  D'       B*PH
%          D         -gamma*I C
%          H*PB      C*       A*P*H+H'*PA]

disp('Constructing the derivative inequality');

Dop = [Eop'*(Pop*Aop)+Aop'*(Pop*Eop)]; 
    


%getting multiplier and kernel for derivative lyapunov function
disp('Parameterizing the negative operators');


% [prog, Deop] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd2,options2);
% [prog, De2op] = sos_posopvar(prog, [nx1 ,nx2],X,s,theta,dd3,options3);


% derivative negativity
% constraints
disp('Setting up the equality constraints');
tic
% prog = sosopeq(prog,Deop+De2op+Dop); %Dop=-Deop-De2op
prog = sosopineq(prog,Dop,opts);
toc
% choosing a different solver if needed
% option.solver = 'mosek'; 

%solving the sos program
prog = sossolve(prog); 

