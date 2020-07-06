%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% executive_PIETOOLS_stability_dual.m     PIETOOLS 2020a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a stability analysis for the dual PIE System defined
% by the 4-PI operator representation
% Top^* \dot x(t)=Aop^* x(t)
%
% NOTE: Stability of the dual system is known to be equivalent to stability of
% the primal system.
%
% If any other parts of the PIE are present, these are ignored. Both Top
% and Aop must be properly defined for the script to function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined externally:
%
% Top,Aop - 4-PI operators, typically defined by the conversion script
%
% eppos,eppos2,epneg % stricness terms typically defined by the solver script
%
% sos_opts - options for the SOSSOLVER (e.g. sdp solver), typically defined by the solver script
%
% dd1,dd2,dd3,opts,options1,options2,options - accuracy settings, typically defined by the settings script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=Aop.I;                         % retrieve the domain from Aop
nx1=Aop.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=Aop.dim(2,1);                % retrieve the number of distributed states from Aop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the Lyapunov 
% function candidate
disp('Parameterizing Positive Lyapunov Operator using specified options');

[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, [nx1 ,nx2],X,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the Lyapunov Inequality
%
% i.e. - Assemble the big operator
% Pheq = [ A*P*T'+T*P*A']

disp('Constructing the Negativity Constraint');

Dop = [Top*Pop*(Aop')+Aop*Pop*(Top')+epneg*Top*Top']; 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('Enforcing the Negativity Constraint');

% derivative negativity
% constraints
if sosineq_on
    disp('Using lpi_ineq');
    prog = lpi_ineq(prog,-Dop,opts);
else
    disp('Using an Equality constraint');

    [prog, De1op] = poslpivar(prog, [nx1, nx2],X,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,[nx1, nx2],X, dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
end

disp('Solving the LPI using the specified SDP solver');
%solving the sos program
prog = sossolve(prog,sos_opts); 
% Conclusion:
if exist('prog.solinfo', 'var')
    if norm(prog.solinfo.info.feasratio-1)<=.3 && ~prog.solinfo.info.numerr
        disp('The System of equations was successfully solved.')
    elseif norm(prog.solinfo.info.feasratio-1)<=.3 && prog.solinfo.info.numerr
        disp('The System  of equations was successfully solved. However, Double-check the precision.')
    elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
        disp('The System is  of equations not solved.')
    else
        disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
    end
else
    disp('ODE-PDE converted to PIE. No problem solved because executive file was not selected');
end
