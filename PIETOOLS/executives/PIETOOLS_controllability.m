function [prog, lam,Wc] = PIETOOLS_controllability(PIE,settings,varargin)
if PIE.dim==2
    disp('Requires creating H2 norm script for 2D systems')
    return
end

dd1 = settings.dd1;
dd12 = settings.dd12;
sos_opts = settings.sos_opts;
options1 = settings.options1;
options12 = settings.options12;
override1 = settings.override1;
eppos = settings.eppos;
epneg = settings.epneg;
eppos2 = settings.eppos2;
ddZ = settings.ddZ;
sosineq_on = settings.sosineq_on;
if sosineq_on
    opts = settings.opts;
else
    override2 = settings.override2;
    options2 = settings.options2;
    options3 = settings.options3;
    dd2 = settings.dd2;
    dd3 = settings.dd3;
end

% Dumping relevant 4-PI operators to the workspace
Aop=PIE.A;
Top=PIE.T;
B1op=PIE.B1;
B2op=PIE.B2;    %TB1op = PIE.Tw;
C1op=PIE.C1;
C2op=PIE.C2;
%D11op=PIE.D11;

fprintf('\n --- Searching for the controllability gramian --- \n')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=Aop.I;                         % retrieve the domain from Aop
nx1=Aop.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=Aop.dim(2,1);                % retrieve the number of distributed states from Aop
nw=B1op.dim(1,2);
nu=B2op.dim(1,2);                % retrieve the number of real-valued disturbances
nz=C1op.dim(1,1);                % retrieve the number of real-valued regulated outputs
ny= C2op.dim(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the storage
% function candidate
disp('- Declaring Gramian using specified options...');

dpvar lam;
prog = sosdecvar(prog,lam);
prog = sossetobj(prog,lam);
prog = sosineq(prog,lam);

[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);
if override1~=1
    [prog, P2op] = poslpivar(prog, [nx1 ,nx2],X,dd12,options12);
    Wop=P1op+P2op;
else
    Wop=P1op;
end


Wop.P = Wop.P+eppos*eye(nx1);
Wop.R.R0 = Wop.R.R0+eppos*eye(nx2);
% W2op.R.R0 = W2op.R.R0+eppos*eye(nx2);
% opvar tmp; tmp.R.R1 = Wop.var1; tmp.R.R2 = Wop.var1;
% Wop = Wop+eppos*tmp;
% W2op = W2op+eppos*tmp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Using the controlability gramian

disp('- Constructing the Negativity Constraint...');

% lam = 0;
opvar I; I.dim = PIE.T.dim; I.R.R0 = eye(nx2); I.P = eye(nx1);

Dop =  B1op*B1op'+lam*(PIE.A-I)*(PIE.A-I)'+0.5*lam*(PIE.T+PIE.A+I)*(PIE.T+PIE.A+I)'...
            -PIE.A*Wop*PIE.A'-0.5*(PIE.T-PIE.A+I)*Wop*(PIE.T-PIE.A+I)' - lam*I;
D2op =  B1op*B1op'+(PIE.A-I)*Wop*(PIE.A-I)'+0.5*(PIE.T+PIE.A+I)*Wop*(PIE.T+PIE.A+I)'...
            -lam*PIE.A*PIE.A'-0.5*lam*(PIE.T-PIE.A+I)*(PIE.T-PIE.A+I)' - Wop;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Negativity Constraint...');
if sosineq_on
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq(prog,-Dop,opts);
else
    disp('  - Using an Equality constraint...');
    [prog, De1op] = poslpivar(prog, PIE.T.dim(:,1),X,dd2,options2);
    [prog, De3op] = poslpivar(prog, PIE.T.dim(:,1),X,dd2,options2);
    if override2~=1
        [prog, De2op] = poslpivar(prog,PIE.T.dim(:,1),X, dd3,options3);
        [prog, De4op] = poslpivar(prog,PIE.T.dim(:,1),X, dd3,options3);
        Deop=De1op+De2op;
        Deop2=De3op+De4op;
    else
        Deop=De1op;
        Deop2=De3op;
    end
    prog = lpi_eq(prog,Dop-Deop,'symmetric'); %Dop=-Deop
    prog = lpi_eq(prog,D2op+Deop2,'symmetric'); %Dop=Deop
end


%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = sossolve(prog,sos_opts);
lam = sosgetsol(prog,lam);
Wc = getsol_lpivar(prog,Wop);
end