% This document illiustrates the use of PIETOOLS to verify the 
% stability of a reaction-diffusion equation. We walk through this code in
% the file "PIETOOLS_2021a__How_to_get_started.pdf". 

clc; clear;

% We consider a reaction-diffusion equation
%           \dot(x)(t,s) = \lambda x(t,s) +  \partial_{s}^2 x(t,s),   
% with BCs
%           x(t,a) = 0,  x(t,b) = 0,
% on a domain
%           s \in [a,b] = [0,1]

% We know this PDE to be stable for \lambda < pi^2 \approx 9.8696.
% We wish to test whether this PDE is indeed stable for lambda = 9.
% We may do this in PIETOOLS by taking the following steps:



%%
%%% STEP 1: Specify the Model
% In order for us to study the system with PIETOOLS, we first need to
% describe the model in a standardized format: using the structure "PDE".
% The easiest way to build this structure is with the GUI, as discussed 
% in the manual. However, we can also specify the necessary fields of the
% "PDE" structure directly using the term-based input format, as we will be 
% doing in this code. For a full overview of the fields making up the 
% structure "PDE", please consult the manual.

% The first thing we need to specify is the number of states considered in
% the PDE, separated by order of differentiability. In our case, the PDE
% concerns only a single state variable, which must be twice
% differentiable. We specify this by setting:

PDE.n.n_pde = [0,0,1];

% Next, we specify the PDE dynamics, contained in "PDE.PDE". Since our PDE
% is not coupled to any ODE, nor do any boundary values of the state appear
% in the evolution equation, the only field we need to describe is
% "PDE.PDE.A". This field is a cell, of which each element is a struct
% specifying a particular term in the PDE. To specify 
% "\dot(x)(t,s) = \lambda x(t,s)", with "\lambda = 9"
% we take

lambda = 9;
PDE.PDE.A{1}.Lstate = 2;        PDE.PDE.A{1}.Rstate = 2;
PDE.PDE.A{1}.D = 0;             PDE.PDE.A{1}.I = 0;
PDE.PDE.A{1}.coeff = lambda;

% Here, "Lstate=2" refers to the fact that we are adding a term to the
% evolution equation of our 2nd state component. Note here that we are
% refering to the 2nd state component as the vector of state variables
% differentiable up to order 2, which in this case is the full state x_2=x.
% Since the term we are adding also concerns the 2nd (and only) state
% component, we set "Rstate=2". Then, specifying "D=0" assures we take no
% derivative of this component, and taking "I=0" assures the term involves
% no integral. Finally, by setting "coeff=lambda", we multiply are term
% with "lambda" before adding it to the evolution equation, so that our PDE
% becomes
% \dot{x}_2(t,s) = lambda * x_2(t,s),
% where x_2 = x;
%
% Next, we add the term "\partail_{s}^2 x(t,s)" to our PDE, using

PDE.PDE.A{2}.Lstate = 2;        PDE.PDE.A{2}.Rstate = 2;
PDE.PDE.A{2}.D = 2;             PDE.PDE.A{2}.I = 0;
PDE.PDE.A{2}.coeff = 1;

% Here, we once again set "Lstate=2" and "Rstate=2", as we are still 
% considering the evolution equation of our 2nd state component, and this
% term also involves the 2nd state equation. However, in this case, we want
% to take the 2nd order derivative of our state, which we enforce with
% "D=0". We are still taking no integral though ("I=0"), and in this case
% we are not multiplying the term with any value, so we set "coeff=1".
% Running this line on its own would suggest a PDE
% \dot{x}_2(t,s) = \partial_{s}^2 x_2(t,s),
% but since we have already specified a term PDE.PDE.A{1}, we get the full
% equation
% \dot{x}_2(t,s) = \lambda x_2(t,s) + \partial_{s}^2 x_2(t,s).
% Note that it is important to assign different terms to different elements
% of the cell PDE.PDE.A, but the order of these cells does not matter.

% Having described the evolution equation, it remains to specify the
% boundary conditions. These are described by the field "PDE.PDE.BC",
% where, since we have no ODE connection, nor does an integral of the PDE
% state appear in the BCs, we only need to define "PDE.PDE.BC.Ebb". To
% enforce the first condition, we set

PDE.BC.Ebb{1}.Rstate = 2;
PDE.BC.Ebb{1}.D = 0;        PDE.BC.Ebb{1}.delta = 0;
PDE.BC.Ebb{1}.coeff = [1;0];

% Here, we set "Rstate=2" and "D=0", as the term we are adding to the 
% boundary equation involves the 0th derivative of our 2nd state component. 
% By taking "delta=0", we require this component to be evaluated at the
% lower boundary s = a, so that we add a term "x_2(t,a) = x(t,a)" to the
% boundary equation. Finally, since this term occurs only in the first BC,
% we multiply it with "coeff=[1;0]", so that our BCs become:
% 0 = 1 * x_2(t,a) = x(t,a)
% 0 = 0 * x_2(t,a) = 0
% Note here that the number of rows in "Ebb{i}.coeff" for any cell i must 
% always correspond to the total number of boundary conditions that you 
% enforce, so if you don't want your term to appear in a condition you must
% set the corresponding coefficient to zero.
% Now, to enforce the second BC, we take

PDE.BC.Ebb{2}.Rstate = 2;
PDE.BC.Ebb{2}.D = 0;        PDE.BC.Ebb{2}.delta = 1;
PDE.BC.Ebb{2}.coeff = [0;1];

% In this case, we set "delta=1", to specify a term at the upper boundary
% s=b. In particular, since "Rstate=2" and "D=0", we are adding 
% "x_2(t,b) = x(t,b)" to the boundary equation. Since we only want this
% term to appear in the second BC, we set "coeff=[0;1]", so that our BCs
% become
% 0 = x(t,a) + 0 * x_2(t,b) = x(t,a)
% 0 = 0      + 1 * x_2(t,b) = x(t,b)
% et voila, we have our boundary conditions! As with the cell
% PDE.PDE.A, though, be careful that each term you add to the BCs must be 
% specified with a different element of the cell PDE.BC.Ebb.

% We have now specified our evolution equation and boundary conditions. The
% only thing left to do is specifying the domain of our system,

PDE.dom = [0,1];

% and we are done!



%%
%%% STEP 2: Construct and view the PIE Representation
% Once the model is defined with the struct "PDE" (which, remember, you can
% also do using the GUI), we can convert it to a format we can easily
% analyze, namely the PIE! This conversion can be done using PIETOOLS,
% simply calling:

PIE = convert_PIETOOLS_PDE_terms(PDE);

% This script will take our "PDE" structure, fill in some of the blanks,
% and compute the PI operators describing the equivalent PIE. Note that, in
% the process, you may receive several warnings that you have not defined
% certain fields of the "PDE" structure, and that these are defaulted to
% zero. This is because the conversion requires all the fields of "PDE" to
% be specified, whilst we only specified the nonzero terms (the terms that
% actually appear in the PDE), so the converter fills in all the zeros.
% This also means that the struct "PDE" is adjusted in the process,
% with (for example) elements of the cell "PDE.PDE.A" being assigned
% different indices, so be careful if you want to make changes to your PDE
% afterwards.

% Having converted the PDE to a PIE, let's view the PI operators describing
% this equation. Since we are only considering a "pure" PDE, the PIE will
% be of the form 
% Top * \dot{\hat{x}}(t,s) = Aop * \hat{x}(t,s)
% where Top and Aop are 4-PI operators, and \hat{x} = \partial_s^2 x, 
% noting that "x" is differentiable up to order 2. To view the PI operators
% describing this equation, we call

Top = PIE.T;
Aop = PIE.A;

% Note that, since we have no ODE interconnection, only the 3-PI components
% Top.R and Aop.R of 4-PI operators Top and Aop will be nonzero.



%%
%%% STEP 3: Solve the Stability Problem
% Having converted the PDE to a PIE, we can do what we came for: test
% stability of our PDE. This can be done using the executive file
% "executive_PIETOOLS_stability". However, we can also do this manually, as
% described below:

% To test stability of our PDE, we test feasibility of the PIE
% Pop > 0
% Top'*Pop*Aop + Aop'*Pop*Top < 0
% where Pop is a PI operator. Accordingly, these constraints describe an
% LPI, with as decision variable Pop. To specify and solve this LPI, we
% describe it as an SOS program, which we initialize as follows: 

prog = sosprogram([Top.var1; Top.var2]);

% Here, Top.var1 and Top.var2 denote the independent variables in the
% program, which will be s and theta. To declare our decision variable, the
% positive PI operator Pop, we take:

[prog, Pop] = poslpivar(prog, [0, 1],[0, 1]);

% Here, the first input [0, 1] corresponds to the dimensionality of the 
% problem, with 0 corresponding to the number of state variables in \R (the
% number of ODE state variables), and 1 corresponding to the number of 
% state variables in L_2 (the number of PDE state variables). The second
% input [0,1] corresponds to the domain [a,b] in which the independent
% variables live. Running this line of code, we enhance our program with a
% positive semidefinite decision variable Pop, which we make strictly
% positive definite by taking

Pop = Pop+.0001; 

% Then, it remains to specify the second PI inequality, corresponding to 
% the actual stability test, which we do by taking

options.psatz = 1;
prog = lpi_ineq(prog,-(Top'*Pop*Aop+Aop'*Pop*Top), options);

% Here, providing the optional input psatz = 0 assures we only require the
% operator "-(Top'*Pop*Aop+Aop'*Pop*Top)" to be positive on the domain
% s \in [a,b] = [0,1]. Note the negative sign in specifying this
% inequality, to assure negativity of the operator 
% "Top'*Pop*Aop+Aop'*Pop*Top". Having defined our SOS program, we can solve
% it by calling

prog = sossolve(prog);

% The structure "prog" now contains all information regarding the solved 
% SOS program that we specified, also telling us whether or not it is
% actually feasible. There are several parameters we can look at to see
% whether the problem was successfully solved (and thus feasible), one of
% which is the "feasratio":

if norm(prog.solinfo.info.feasratio-1)<=.1 
    disp('The system of equations was successfully solved :D')
else
    disp('The system of equations was not solved :(')
end

% If everything is working correctly, you should find that the system was
% successfully solved, proving stability of our PDE. If you were to change
% the value of \lambda to 10, however, you will find that the system of
% equations cannot be solved, as the PDE in this case is unstable! 

% Using the basic methodology as outlined in this file, PIETOOLS can be 
% used to solve a wide variety of problems for a wide variety of PDEs (and
% DDEs). For most general problems, such as stability tests or Hinf-gain 
% consider using the PIETOOLS_PDE (or PIETOOLS_DDE) file, which will call
% the appropriate executable once you have specified your problem. For
% implementing your PDE, we suggest using the GUI, or draw inspiration from
% the example files, which describe several systems in the appropriate
% format. If you run into issues, consult the user manual or
% troubleshooting file, and if that doesn't help, contact us (contact info
% is in troubleshooting file). We tried to make this toolbox as easy as
% possible, so let us now if you have any questions or suggestions. 
% Have fun!


