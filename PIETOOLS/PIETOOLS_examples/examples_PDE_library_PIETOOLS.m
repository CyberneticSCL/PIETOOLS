%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIETOOLS.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = examples_PDE_library_PIETOOLS(varargin)
% This library file contains the definition of some common ODE-PDE systems
% drawn from the literature. To use one of the examples, call this function
% with as FIRST input the index of the example you wish to use. Please
% check the table of contents below to see which index corresponds to which
% function. The output of the function then corresponds to a Matlab
% structure describing the desired PDE.
%
% Inputs:
%     (output1,output2) = examples_PDE_library_PIETOOLS(example,option1,option2,option3,...)
%       - example - integer (or possibly double) number (1-27) corresponding to selected example problem    
%       - option1,option2,option3, etc.
%            - option# - 'batch' or 'BATCH' - output# will be the PDE data structure in batch format  
%            - option# - 'terms' or 'TERMS' - output# will be the PDE data structure in terms format
%            - option# - 'gui' or 'GUI' - no output. Executes PIETOOLS PDE GUI
%             and loads saved problem data into GUI.
%            - option# - 'lambda=XXX' overrides the default parameter value
%             lambda with value XXX for the given problem number. The list of allowable parameter 
%             overrides for each example problem is given in the table of
%             contents below. DO NOT ABUSE THIS POWER!
%             NOTE: Changing parameters does not change the GUI file!
%            
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
% 
% Within each type of examples, we have grouped the examples by the type of
% system they represent. Specifically, we have:
%
% 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% 2. Diffusion and Heat-Equation Type Systems
% 3. Beam Type Equations
% 4. Wave Equations
%
% Each example has been defined using both the classical batch input format 
% and the generalized term-based input format. To obtain a PDE struct in
% the classical format, call this function with as second argument the
% value "1" or the character array "batch". If you wish to
% define the PDE in the term-based input format, call this function with as
% second argument the value "2", or the character array "term". You can
% also specify both a value "1" and a value "2" (as second and third 
% argument) if you wish to obtain the PDE in both formats, in which case
% the outputs will be ordered in accordance with your arguments.
%
% In addition, to improve visualization and allow for tinkering, each 
% example has an associated .mat file which can be loaded into the PDE GUI.
% To open the example in the GUI, call the function with an argument valued
% "3" or character array "GUI". Once again, you can load the GUI file in
% addition to producing a batch- or term-based PDE by specifying multiple
% arguments.
%
% The examples are typically called using a line in the main PIETOOLS_PDE.m
% file. Simply save your changes to the library file, uncomment the line
% examples_PDE_library_PIETOOLS.m in PIETOOLS_PDE.m and run PIETOOLS_PDE.m.
% Of course, you can also run the library file directly.
%
% Several examples come with parameters that may be adjusted. Check the
% "Parameters" column in the table of contents below. If you wish to adjust
% one of the parameter values, specify this adjustment using a character 
% array as one of your arguments. For example, add argument "lam = 10" to
% set the parameter "lam" equal to 10. 
% Also, most examples in each problem type include the option to solve the 
% LPI according to their problem type. If you want to invoke a different
% LPI, specify this option after running the function (e.g. set
% "stability = 1", "Hinf_gain=1", etc.).
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
% NOTE: At present, PIETOOLS does not support inputs at the boundary for
% solving the Hinf optimal control problem. Support for this option will be
% included in a future release.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP - 7_1_2020
% Updated to the Nth order format, added 4th order implementation of 
% Timoshenko beam equation. DJ - 05-24-2021
% Adjusted to function format. DJ - 05-31-2021
% Updated zero input case, and fixed some examples. DJ - 12-29-2021
% Update to new terms format and separate files: DJ - 08/25/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process the user inputs
[index,BATCH,TERM,GUI,params] = process_inputs(varargin,nargin);

fprintf([' --- Extracting ODE-PDE example ', num2str(index),' ---\n']);
switch index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Table of Contents
%
%__________________________________________________________________________
%   Problem                                         | Parameters            Notes / References
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
    case 1
% 	PDE: x_{t} = x_{s}                              |   
%   BCs: x(s=0) = 0                                 |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params);
%-----|-------------------------------------------------------------------
    case {2.0,2.1,2.2,2.3,2.4}
%	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1] 
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling
%        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]         | examples 2.1, 2.2, 
%                                                   | 2.3, and 2.4.
	[PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Lamare(index,GUI,params);
%-------------------------------------------------------------------------
    case 3
%   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
%        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
%   The implementation is based on a Linearized Saint–Venant–Exner Model
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diagne(GUI,params);
%--------------------------------------------------------------------------
    case {4.0,4.1,4.2,4.4,4.3,4.5}
% 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}           | Different parameters          Saba 2019 [3] 
%        x2_{t} = sig2*x1 + (1/r2)*x2_{s}           | may be invoked calling  
%   BCs: x1(s=0) = qb*x2(s=0)                       | examples 4.1, 4.2, 
%        x2(s=1) = pb*x1(s=1)                       | 4.3, 4.4, and 4.5.  
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Saba(index,GUI,params);
%--------------------------------------------------------------------------
%        Diffusion and Heat-Equation Type Systems
%--------------------------------------------------------------------------
    case 5
%   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86            (stable for lam < pi^2 = 9.8696)  
%   BCs: x(s=0) = 0,      x(s=1) = 0                |                               Ahmadi 2015 [5] 
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq(GUI,params);
%--------------------------------------------------------------------------
    case 6
%   | PDE: x_{t} = lam*x + x_{ss}                   | lam = 2.466           (unstable for lam > pi^2/4 = 2.467)        
%   | BCs: x(s=0) = 0,      x_{s}(s=1) = 0          |                               Gahlawat 2017 [4]
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat(GUI,params);
%--------------------------------------------------------------------------
    case 7
%   PDE:  x_{t} = a(s)*x_{ss}                       | a = s^3 - s^2 + 2     (unstable for lam > 4.66)            
%                 + b(s)*x_{s}                      | b = 3*s^2 - 2*s               Gahlawat 2017 [4]
%                 + c(s,lam)*x                      | c =-0.5*s^3 + 1.3*s^2 
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            |    - 1.5*s + 0.7 +lam
%                                                   | lam = 4.66                 
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params);
%--------------------------------------------------------------------------
    case 8
% 	ODE:  xo_{t} = A * xo + Bxr * x_{s}(s=a)        | a = 1    b = 1        (stable for lam < pi^2 = 9.8696)
%   PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo      | lam = pi^2-1                  Das 2018 [7] (Example 2)
%   BCs:  x(s=a) = 0,     x(s=b) = 0                | A, Bxr, and Bpv fixed
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params);
%--------------------------------------------------------------------------
    case 9
% 	ODE:  xo_{t} = k*xo                             | k = -1                (stable for k<0)
%   PDE:  x_{t} = x_{ss}                            |
%   BCs:  x(s=0) = 0,     x(s=1) = xo               |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC(GUI,params);
%--------------------------------------------------------------------------
    case 10
%   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2                (stable for k=-2)            
%   PDE: x_{t} = x_{ss}                             |                               Tang 2011 [13]
%   BCs: x(s=0) = -xo,    x(s=1) = k*xo             |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params);
%--------------------------------------------------------------------------
    case 11
% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.93              (stable for R<2.7)
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6]
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi(GUI,params);
%--------------------------------------------------------------------------
    case 12
%   PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = (21+9)            (stable for R<21)
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [0 0 0;                  Ahmadi 2015 [5] Adapted Example B
%                                                   |       s 0 0;
%                                                   |       -s^2 0 0]
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params);
%
%--------------------------------------------------------------------------
%        Beam Type Equations 
%--------------------------------------------------------------------------
    case 13
%   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
%   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |
%   Use states x1 = u_{t}, x2 = u_{ss}.             |
%       =>                                          |
%   PDE: x1_{t} = -c * x2_{ss}                      |
%        x2_{t} = x1_{ss}                           |
%   BCs: x1(s=0) = 0,       x2(s=1) = 0             |
%        x1_{s}(s=0) = 0,   x2_{s}(s=1) = 0         |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq(GUI,params);
%--------------------------------------------------------------------------
%
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
%   
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 14
%   Use states:                                     | k = 1                 Hyperbolic implementation (stable)
%        x1 = w_{t},   x2 = k*aa*g * (w_{s}-phi),   | aa = 1
%        x3 = phi_{t}, x4 = E*II * phi_{s}.         | II = 1
%       =>                                          | g = 1
%   PDE: x1_{t} = (1/r/aa) * x2_{s}                 | E = 1
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3     | r = 1
%        x3_{t} = (1/r/II)*x2 + (1/r/II)*x4_{s}     | 
%        x4_{t} = E*II * x3_{s}                     |
%   BCs: x1(0) = 0,         x2(1) = 0               |
%        x3(0) = 0,         x4(1) = 0               |
%                                                   |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params);
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 15
%   Assume all parameters are 1, and use states:    |                       Hyperbolic/diffusive implementation (unstable)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |
%       =>                                          |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4                           |
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0                          |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params);
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 16
%   We rewrite the system as a single equation:     |                       4th order implementation (stable)
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + w_{tttt}             |
%   Use states:                                     |
%        x1 = w_{ttt},      x2 = w_{t},             |
%        x3 = w_{tt},       x4 = w.                 |
%       =>                                          |
%   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3           |
%                 + gam*x3_{ss}                     |
%        x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2    |
%   BCs: w(s=0) = 0,        w_{s}(s=0) = 0          |
%        w_{ss}(s=1) - w(s=1) = 0                   |
%        w_{sss}(s=1) - w_{s}(s=1) = 0              |
%        w_{t}(s=0) = 0     w_{ts}(s=0) = 0         |
%        w_{tt}(s=0) = 0    w_{tts}(s=0) = 0        |
    if BATCH~=0
        % NOTE: Derivatives of order greater than 2 cannot be represented
        % using the batch format. They can be specified using the GUI.
        disp('Warning: No batch input formulation for this problem exists; reverting to term-based input format')
        BATCH = 0;         TERM = 1;
    end
    PDE_t = PIETOOLS_PDE_Ex_Timoshenko_Beam_3(GUI,params);
%--------------------------------------------------------------------------
%       Wave Equations
%--------------------------------------------------------------------------
    case 17
%   PDE: u_{tt} = u_{ss}                            | k=1                   (stable for k=1)                  
%   BCs: u(s=0) = 0,                                |                               Peet 2019 [8] (Example 8.2)
%        u_{s}(s=1) = -k*u_{t}(s=1)                 |
%   Use states x1 = u_{s}, x2 = u_{t}.              |
%       =>                                          |
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params);
%-----|---------------------------------------------|----------------------
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1                         Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 18
%   Use states x1 = u_{t}, x2 = u.                  | k = 1
%       =>                                          | ad = 1
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params);
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 19
%   Use states x1 = u_{t}, x2 = u, x3 = u_{s}.      | k = 1                
%       =>                                          | ad = 1
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x3(1) = 0                        |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_2(GUI,params);
%
%==========================================================================
%       Hinf-gain analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic/Transport/Balance Type Systems
%--------------------------------------------------------------------------
    case 20
%   PDE: x_{t} = -x_{s} + w(t)                      |                       (gamma = 0.5)
%   BCs: x(s=0) = 0                                 |
%   Out: z(t) = int(x(t,s),s,0,1)                   |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance(GUI,params);
%--------------------------------------------------------------------------
    case 21
%   PDE: phi_{tt} = phi_{ss} + w(t)                 | k = 0.5               (gamma = 2 for k = 0.5)  
%   BCs: phi(s=0) = 0                               |
%        phi_{s}(s=1) = -k*phi_{t}(s=1)             |
%   Out: z(t) = int(phi_{t}(t,s),s,0,1)             |
%   Use states x1 = phi_{s},    x2 = phi_{t}        |
%       =>                                          |
%   PDE: x1_{t} = x2_{s}                            |
%        x2_{t} = x1_{s} + w(t)                     |
%   BCs: x2(0) = 0,     x1(1) + k*x2(1) = 0         |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped(GUI,params);
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 22
%   PDE: x_{t} = x_{ss} + s*w(t)                    |                       (gamma = 0.3333)
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |
%   Out: z(t) = int(x(t,s),s,0,1)                   |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params);
%--------------------------------------------------------------------------
    case 23
%   PDE: x_{t} = A2(s)*x_{ss}                       | A2 = s^3 - s^2 + 2;   (gamma = 15.147 for lamb = 4.6)       
%                + A1(s)*x_{s}                      | A1 = 3*s^2 - 2*s;             Shivakumar 2019 [12] (Example 1)
%                + A0(s)*x + w(t)                   | A0 =-0.5*s^3 +1.3*s^2   (gamma = 23.7-57 using "veryheavy" settings)                
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |     -1.5*s +0.7 +lam
%   Out: z(t) = x(t,1)                              | lam = 4.6;  
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params);
%--------------------------------------------------------------------------
    case 24
%   PDE: xi_{t} = w(t) + lamb*xi        i=1:ne      | lam = (1-1e-2)*pi^2   (gamma = 8.1069 for lamb = (1-1e-2)*pi^2)         
%                 + sum(xk_{ss},k=1,i)              | ne = 1 (state size)           Shivakumar 2019 [12] (Example 3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              |   
%   Out: z(t) = int(x(t,s),s,0,1)                   |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance(GUI,params);
%--------------------------------------------------------------------------
    case 25
%   PDE: x_{t} = Cm*x + (1/R)*x_{ss} + s*w(t)       | R = 2.6               (gamma = 0.8102 for R = 2.6)
%   BCs: x(s=0) = 0,        x(s=1) = 0              | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6]
%   Out: z(t) = int(x1(t,s),s,0,1)                  | 
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance(GUI,params);
%--------------------------------------------------------------------------
    case 26
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + s*w(t)    | R = (21-1e-3)         (gamma =  4.23, for R = 21-1e-3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              | Cm = [0,0,0;                  Shivakumar 2019 [12] (Example 2)
%   Out: z(t) = int(x(t,s),s,0,1)                   |       s,0,0;
%                                                   |       s^2,-s^3,0]
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2(GUI,params);
%     
%==========================================================================
%       Hinf-optimal observer
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
    case 27
%   PDE: x_{t} = x_{s} + w(t)                       | ne = 1 (state size)   (Answer: 1.0012)
%   BCs: x(s=1) = 0                                 |
%   Out: z(t) = int(x(t,s),s,0,1) + w(t)            |
%        y(t) = int(x(t,s),s,0,1)                   |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output(GUI,params);
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 28
%   PDE: x_{t} = x_{ss} + w(t)                      | ne = 1 (state size)   (Answer: 1.0045) 
%   BCs: x(s=0) = 0,        x(s=1) = 0              |
%   Out: z(t) = int(x(t,s),s,0,1) + w(t)            |
%        y(t) = x_{s}(s=1)                          |
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output(GUI,params);
%
%==========================================================================
%       Hinf-optimal controller
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 29
% % PDE                            x_{t} = lam*x + x_{ss} + u(t) + w(t)
% % With BCs                         x(s=0) = 0
% %                                         x(s=1) = 0
% % and regulated output  z = [ int(x(s,t),s,0,1);u]
    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       2D Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Advection Equation
%--------------------------------------------------------------------------
    case 31
%   PDE: x_{t} = c1*x_{s1} + c2*x_{s2}              | c1 = 1; c2 = 1;
%   BCs: x(s1=0) = 0,   x(s2=0) = 0;                |   ne = 1 (state size)
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Transport_Eq(GUI,params);
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 32
%   PDE: x_{t} = lam*x + c1*x_{(2,0)}               | lam = 19;             (stable for lam <= 2*pi^2) 
%                           + c2*x_{(0,2)}          | c1 = 1;   c2 = 1;             Holmes, 1994 [14] Eq. (14)   
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;                |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params);
%--------------------------------------------------------------------------
    case 33
%   PDE: x_{t} = C*(x_{(2,0} + x_{(0,2)})           | C = 1;
%                 - b1*x_{(1,0)} - b2*x_{(0,1)}     | b1 = 0.5; b2 = 2;             Holmes, 1994 [14] Eq. (2)
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;                |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params);
%--------------------------------------------------------------------------
    case 34
%   PDE: x_{t} = -1/(2*lam) * x_{tt}                | lam = 1;
%                 +(C^2/2*lam)*(x_{(2,0}            | C = 1;                        Holmes, 1994 [14] Eq. (3)
%                                 + x_{(0,2)})      | ne = 1 (state size) 
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | 
%        x(s1=1) = 0,   x(s2=1) = 0;                |
%   Use states x1 = x,    x2 = x_{t}                |
%       =>                                          |
%   PDE: x1_{t} = x2                                |
%        x2_{t} = C^2*(x1_{(2,0} + x1_{(0,2)})      |
%                   - 2*lam*x2                      |
%   BCs: x1(s1=0) = 0,   x1(s2=0) = 0,              | 
%        x1(s1=1) = 0,   x1(s2=1) = 0;              |
%        x2(s1=0) = 0,   x2(s2=0) = 0,              | 
%        x2(s1=1) = 0,   x2(s2=1) = 0;              |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params);
%--------------------------------------------------------------------------
    case 35
%   PDE: x_{t} = a*x                                | a = 4.9;              (stable for a <= 0.5*pi^2) 
%                 + b1*x_{(1,0)} + b2*x_{(0,1)}     | b1 = 0;   b2 = 0;             Demetriou, 2019 [17]   
%                  + c1*x_{(2,0)} + c2*x_{(0,2)}    | c1 = 1;   c2 = 1;
%   BCs: x(s1=0) = 0,       x_{s1}(s1=1) = 0,       | ne = 1 (state size)
%        x(s2=0) = 0,       x_{s2}(s2=1) = 0;       |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Parabolic_Eq(GUI,params);
%--------------------------------------------------------------------------
    case 36
%   PDE: x_{tt} = c1*u_{(2,0)} + c2*x_{(0,2)}       | c1 = 1;   c2 = 1;
%   BCs: x(s1=0) = 0;       x(s2=0) = 0;            | ne = 1; (state size)
%        x(s1=1) = 0;       x(s2=1) = 0;            |
%   Use states x1 = x,    x2 = x_{t}                |
%       =>                                          |
%   PDE: x1_{t} = x2;                               | 
%        x2_{t} = c1*x1_{(2,0)} + c2*x1_{(0,2)};    |
%   BCs: x1(s1=0) = 0;      x1(s2=0) = 0;           |
%        x1(s1=1) = 0;      x1(s2=1) = 0;           |
%        x2(s1=0) = 0;      x2(s2=0) = 0;           |
%        x2(s1=1) = 0;      x2(s2=1) = 0;           |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Wave_Eq(GUI,params);
%--------------------------------------------------------------------------
    case 37
%   PDE: p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}    | M = 0.1;
%        v1_{t} = -s2*v1_{s1} -v2 -(1/M^2)*p_{s1}   | lam = 1; nu = 1;              Antonelli, 2021 [17]
%                  +nu*(v1_{s1s1} + v1_{s2s2})      |
%                   +lam*(v1_{s1s1} + v2_{s2s1})    |
%        v2_{t} = -s2*v2_{s1} -(1/M^2)*p_{s2}       |
%                  +nu*(v2_{s1s1} + v2_{s2s2})      |
%                   +lam*(v1_{s1s2} + v2_{s2s2})    |
%   BCs: p(s1=0) = 0;      p(s2=0) = 0;             |
%        v(s1=0) = 0;      v(s2=0) = 0;             |
%        v(s1=1) = 0;      v(s2=1) = 0;             |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_NS_Antonelli(GUI,params);
%--------------------------------------------------------------------------
    case 38
%   ODE: x1_{t} = (A+BK)*x1 + B*x2(s1=0,s2=0)       | A = I; B = I;
%   PDE: x2_{t} = c1*x_{(2,0)} + c2*x_{(0,2)}       | K = -2*I;
%   BCs: x2_{(1,0)}(s1=0) = 0;  x2(s1=1) = 0;       |
%        x2_{(0,1)}(s2=0) = 0;  x2(s2=1) = 0;       |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params);
%--------------------------------------------------------------------------
    case 39
%   PDE: x_{t} = x_{s1s1} + lam * x;    s1 in [0,1] | tau = 1;
%   BCs: x(t,s1=0) = 0;                             | lam = 1;                      Kristic, 2009 [18] 
%        x(t,s1=1) = u(t-tau);                      |
%   Set x1(t)=x(t) and introduce delayed state:     |
%        x2_{t} = x2_{s2};          s2 in [1,1+tau] |
%        u(t) = x2(t,1+D);                          |
%       =>                                          |
%   PDE: x1_{t} = x1_{s1s1} + lam * x1;             |
%        x2_{t} = x2_{s2};                          |
%   BCs: x1(t,s1=0)     = 0;                        |
%        x1(t,s1=1)     = x2(t,s2=1);               |
%        x2(t,s2=1+tay) = u(t);                     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
    [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params);
%--------------------------------------------------------------------------
    case 40
%   ODE: X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0) | tau = 1;
%   PDE: x_{t} = x_{s1s1} + a*x + a2*x(t-tau);      | A = 1;    A1 = 0.5;           Kang, 2017 [19] 
%   BCs: x_{s1}(t,s1=0) = 0;                        | B = -1;
%        x(t,s1=1) = u(t);                          | a = 1;    a2 = 0.5;
%   Set x1(t)=X(t), x2(t,s1)=x(t), and introduce    |
%   transport equation to incorporate delay:        |
%        x3_{t} = x3_{s2}                           |
%        x4_{t} = x4_{s2}           s2 in [1,1+tau] |
%        x3(t,s2=1+tau) = x1(t)                     |
%        x4(t,s1,s2=1+tau) = x2(t,s1)               |
%       =>                                          |
%   ODE:    x1_{t} = A*x1(t) + A1*x3(t,1)           |
%                               + B*x2(t,s1=0);     |
%   PDE:    x2_{t} = x2_{s1s1} + a*x2               |
%                               + a2*x4(t,s1,s2=1); |
%           x3_{t} = x3_{s2};                       |
%           x4_{t} = x4_{s2};                       |
%   BCs:    x2_{s1}(t,s1=0) = 0;                    |
%           x2(t,s1=1) = u(t);                      |
%           x3(t,s2=1+tau) = x1(t);                 |
%           x4(t,s1,s2=1+tau) = x2(t,s1);           |    
%           x4_{s1}(t,s1=0,s2) = 0;
%           x4(t,s1=1,s2) = u(t);
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
    [PDE_t] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params);
%--------------------------------------------------------------------------
    case 41
%   PDE:    x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1)  | tau = 1;
%                       - a1*x(t-tau,s1);           | c = 1;
%   BCs:    x(t,s1=0) = 0;                          | a0 = 1.9;
%           x(t,s1=pi) = 0;                         | a1 = 1;
%   Set x1=x, and x2(t,s1,s2)=x1(t-s2,s1).          |
%   Introduce transport equation:                   | For c=a1=1, a0=1.9,           Caliskan, 2009 [20]  
%           x1(t)  = x2(t,0);                       | stable for tau<1.0347
%           x2_{t} = -x2_{s2},  s2 in [0,tau];      |
%   Then                                            |
%   PDE:    x1_{t}  = c*x1_{s1s1}(t,s1)+a0*x1(t,s1) |
%                       - a1*x2(t,s1,s2=tau);       |
%           x2_{t}  = -x2_{s2};                     |
%   BCs:    x1(t,s1=0)     = 0;                     |
%           x1(t,s1=pi)    = 0;                     |
%           x2(t,s1,s2=0)  = x1(t,s1);              |
%           x2(t,s1=0,s2)  = 0;                     |
%           x2(t,s1=pi,s2) = 0;                     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
    [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params);
%--------------------------------------------------------------------------
    case 42 
%   PDE:    x_{tt}  = x_{s1s1}(t,s1);   s1 in [0,1] | tau = 1;  
%   BCs:    x(t,0) = 0;                             | k = 1;
%           x_{s1}(t,1) = -k*(1-mu)*x_{t}(t,1)      | mu = 0.4;
%                           - k*mu*x_{t}(t-tau,1);  |
% Introduce:                                        | Stable for mu<0.5             Xu, 2006 [21]  
%     x1(t) = - k * (1-mu) * \dot{u}(1,t)           |
%               - k*mu*\dot{u}(1,t-tau);            |
%     x2(t,s1) = x(t,s1);                           |
%     x3(t,s1) = x_{t}(t,s1);                       |
%     x4(t,s2) = x(t-tau*s2,1);                     |
%     x5(t,s2) = x_{t}(t-tau*s2,1);                 |
% Then:                                             |
%   ODE:    x1_{t}(t)    = x2_{s1}(t,1);            |
%   PDE:    x2_{t}(t,s1) = x3(t,s1);                |    s1 in (0,1)
%           x3_{t}(t,s1) = x2_{s1s1}(t,s1);         |       
%           x4_{t}(t,s2) = -(1/tau) x4_{s2}(t,s2);  |    s2 in (0,1)
%           x5_{t}(t,s2) = -(1/tau) x5_{s2}(t,s2);  |
%   BCs:    x2(t,0) = 0;                            |                         
%           x3(t,0) = 0;                            |
%           x4(t,0) = x2(t,1);                      |
%           x5(t,0) = x3(t,1);                      |
%           x1(t) = - k * (1-mu) * x2(t,1)          |
%                           - k * mu * x4(t,1);     |
%           x2_{s1}(t,1) = -k*(1-mu)*x3_{t}(t,1)    |
%                               - k*mu*x5_{t}(t,1); |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
    [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(GUI,params);
% %---------------------------------------------------------------------% %
% 
%==========================================================================
%       Additional Examples (Undocumented)
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
% 101.| PDE: x_{t} = x_{ss} + s*w(t)                | N = 8,    
%     | BCs: x(s=0) = 0,       x_{s}(s=1)=0         | x0 =@(s) s^2-2*s
%     | Out: z(t) = int(x(t,s),s,0,1)               | w =@(t) sin(t+eps)./(t+eps)
%     |                                             | t_int = [0 2*pi]
%--------------------------------------------------------------------------
% 102.| PDE: xi_{t} = lamb*xi + w(t)        i=1:ne  | ne = 1,    
%     |               + sum(xk_{ss},k=1,i) + w(t)   | lam = (1-1e-2)*pi^2  
%     | BCs: x(s=0) = 0,       x_{s}(s=1)=0         | N = 8
%     | Out: z(t) = int(sum(xi(t,s),i=1,ne),s,0,1)  | x0 =@(s) 0
%     |                                             | w =@(t) sin(t+eps)./(t+eps)
%     |                                             | t_int = [0 2*pi]
%__________________________________________________________________________
%==========================================================================







case index==100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional examples (undocumented) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % %
% % % % %------------------------------------------------------------------

p = 1; q = 0; qc = 3;

A = [ 0,   -1/4, -1/5,  1/5,  1/6;
      1/2,  1,   -4,    9/2,  7/2;
     -9/4, -1/2, -14,   23,   16;
     -1/5, -1/2, -11/4, 1/10, 5/4;
     -4/3, -4/3, -9,    9,    5/2];

B = [-7/2, -3/2, -1/10, 1/2, 1]';
C = [1/10, -1/3, -4, 7/8, 7/8];
    
 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if TERM~=0
 disp('No terms-based input format available for this system, using batch format instead.')
 BATCH = 1;
end
if BATCH~=0
 %%% Batch input format
 PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nx = 5;   PDE_b.nu = 0;
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
 PDE_b.dom = [0,1];

 PDE_b.A2 = p;   PDE_b.A0 = qc-q;
 PDE_b.B22 = 0;    
 
 PDE_b.A = A; 

 % (3)
 PDE_b.E0 = [B, zeros(5,3)];
 PDE_b.B = [0 0 1 0;
            0 1 0 0];
 PDE_b.Bx = [zeros(1,5);C];

%  % (4)
%  PDE_b.E0 = [B, zeros(5,3)];
%  PDE_b.B = [0 0 1 0;
%             0 0 0 1];
%  PDE_b.Bx = [zeros(1,5);C];

end
if GUI~=0
 disp('No GUI representation available for this system.')
end



    case 101
% % % % %------------------------------------------------------------------
% % % % % 
% % % % %------------------------------------------------------------------

p = 1; q = 0; qc = 3;

A = [-1/4, -1/6, 2, 1, 1/12;
-3/2, -3/2, 5, 5, 1/6;
3/2, -4, -15/2, -5, -1/3;
-13/2, 22, 22, -14, -1/2;
1/7, -1/2, -1/2, 1/5, -5/2];

B = [-5/4, 2/3, 1/6, -1/6, 0]';
C = [-2/5, -5/4, 3/2, 1/3, 1/40];
    
 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if TERM~=0
 disp('No terms-based input format available for this system, using batch format instead.')
 BATCH = 1;
end
if BATCH~=0
 %%% Batch input format
 PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nx = 5;   PDE_b.nu = 0;
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
 PDE_b.dom = [0,1];

 PDE_b.A2 = p;   PDE_b.A0 = qc-q;
 PDE_b.B22 = 0;    
 
 PDE_b.A = A; 
 
% %  % (18)
%  PDE_b.E0 = [zeros(5,2),B, zeros(5,1)];
%  PDE_b.B = [1 0 0 0;
%             0 1 0 0];
%  PDE_b.Bx = [zeros(1,5);C];

% (19)
 PDE_b.E0 = [zeros(5,2),B, zeros(5,1)];
 PDE_b.B = [1 0 0 0;
            0 0 0 1];
 PDE_b.Bx = [zeros(1,5);C];

end
if GUI~=0
 disp('No GUI representation available for this system.')
end

%__________________________________________________________________________
%==========================================================================
end

% % % Processing of the outputs.
% Check if the number of desired outputs is reasonable
if nargout>0 && BATCH==0 && TERM==0
    error('No PDE structure has been produced. Use "Get PDE Object" in the GUI or specify ''batch'' or ''terms'' to obtain corresponding structure')
elseif nargout==2 && max(BATCH,TERM)==1
    error('To obtain PDE struct in both batch and term-based format, submit both ''batch'' and ''terms'' as arguments');
elseif nargout>2
    error('At most 2 outputs can be produced: A PDE struct in batch format, and a PDE struct in term-based format');
end

% Specify the outputs of the function
if BATCH~=0
    varargout{BATCH} = PDE_b;
end
if TERM~=0
    PDE_t = initialize_PIETOOLS_PDE(PDE_t,true);
    display_PDE(PDE_t);
    varargout{TERM} = PDE_t;
end

% Check if the user wants to run the executive
userinp = input('\n Would you like to run the executive associated to this problem? (y/n) \n ---> ','s');
if strcmpi(userinp,'y') || strcmpi(userinp,'yes')
    PIE = convert_PIETOOLS_PDE(varargout{1});
    assignin('base','PIE',PIE);
    evalin('base','PIETOOLS_auto_execute');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % [1] - 
% @article{lamare2016optimisation,
%   title={An optimisation approach for stability analysis and controller synthesis of linear hyperbolic systems},
%   author={Lamare, Pierre-Olivier and Girard, Antoine and Prieur, Christophe},
%   journal={ESAIM: Control, Optimisation and Calculus of Variations},
%   volume={22},
%   number={4},
%   pages={1236--1263},
%   year={2016},
%   publisher={EDP Sciences}
% }
%
% % [2] - 
% @article{diagne2012lyapunov,
%   title={Lyapunov exponential stability of 1-D linear hyperbolic systems of balance laws},
%   author={Diagne, Ababacar and Bastin, Georges and Coron, Jean-Michel},
%   journal={Automatica},
%   volume={48},
%   number={1},
%   pages={109--114},
%   year={2012},
%   publisher={Elsevier}
% }
%
% % [3] - 
% @article{saba2019stability,
%   title={Stability Analysis for a Class of Linear 2x2 Hyperbolic PDEs Using a Backstepping Transform},
%   author={Saba, David Bou and Argomedo, Federico Bribiesca and Auriol, Jean and Di Loreto, Michael and Di Meglio, Florent},
%   journal={IEEE Transactions on Automatic Control},
%   year={2019},
%   publisher={IEEE}
% }
%
% % [4] - 
% @article{gahlawat2016convex,
%   title={A convex sum-of-squares approach to analysis, state feedback and output feedback control of parabolic PDEs},
%   author={Gahlawat, Aditya and Peet, Matthew M},
%   journal={IEEE Transactions on Automatic Control},
%   volume={62},
%   number={4},
%   pages={1636--1651},
%   year={2016},
%   publisher={IEEE}
% }
% 
% % [5] - 
% @article{valmorbida2015stability,
%   title={Stability analysis for a class of partial differential equations via semidefinite programming},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   journal={IEEE Transactions on Automatic Control},
%   volume={61},
%   number={6},
%   pages={1649--1654},
%   year={2015},
%   publisher={IEEE}
% }
% 
% % [6] - 
% @inproceedings{valmorbida2014semi,
%   title={Semi-definite programming and functional inequalities for distributed parameter systems},
%   author={Valmorbida, Giorgio and Ahmadi, Mohamadreza and Papachristodoulou, Antonis},
%   booktitle={53rd IEEE conference on decision and control},
%   pages={4304--4309},
%   year={2014},
%   organization={IEEE}
% }
%
% % [7] -
% @article{das2018representation,
% author = {Das, Amritam and Shivakumar, Sachin and Weiland, Siep and Peet, Matthew},
% year = {2018},
% pages = {},
% title = {Representation and Stability Analysis of PDE-ODE Coupled Systems}
% }
%
% % [8] - 
% @article{peet2019discussion,
%   title={Discussion paper: A new mathematical framework for representation and analysis of coupled pdes},
%   author={Peet, Matthew M and Shivakumar, Sachin and Das, Amritam and Weiland, Seip},
%   journal={IFAC-PapersOnLine},
%   volume={52},
%   number={2},
%   pages={132--137},
%   year={2019},
%   publisher={Elsevier}
% }
% 
% % [9] - 
% @article{datko1986example,
%   title={An example on the effect of time delays in boundary feedback stabilization of wave equations},
%   author={Datko, Richard and Lagnese, John and Polis, MP},
%   journal={SIAM Journal on Control and Optimization},
%   volume={24},
%   number={1},
%   pages={152--156},
%   year={1986},
%   publisher={SIAM}
% }
%
% % [12] -
% @inproceedings{shivakumar2019computing,
%   title={Computing input-ouput properties of coupled linear pde systems},
%   author={Shivakumar, Sachin and Peet, Matthew M},
%   booktitle={2019 American Control Conference (ACC)},
%   pages={606--613},
%   year={2019},
%   organization={IEEE}
% }
%
% % [13] -
% @article{tang2011state,
% title = {State and output feedback boundary control for a coupled PDE–ODE system},
% journal = {Systems & Control Letters},
% volume = {60},
% number = {8},
% pages = {540-545},
% year = {2011},
% issn = {0167-6911},
% doi = {https://doi.org/10.1016/j.sysconle.2011.04.011},
% url = {https://www.sciencedirect.com/science/article/pii/S0167691111000922},
% author = {Shuxia Tang and Chengkang Xie},
% keywords = {Coupled system, PDEs, Boundary control, Output feedback},
% abstract = {This note is devoted to stabilizing a coupled PDE–ODE system with interaction at the interface. First, a state feedback boundary controller is designed, and the system is transformed into an exponentially stable PDE–ODE cascade with an invertible integral transformation, where PDE backstepping is employed. Moreover, the solution to the resulting closed-loop system is derived explicitly. Second, an observer is proposed, which is proved to exhibit good performance in estimating the original coupled system, and then an output feedback boundary controller is obtained. For both the state and output feedback boundary controllers, exponential stability analyses in the sense of the corresponding norms for the resulting closed-loop systems are provided. The boundary controller and observer for a scalar coupled PDE–ODE system as well as the solutions to the closed-loop systems are given explicitly.}
% }
%
% % [14] - 
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.
%
% % [15] -
% @inproceedings{meyer2015stability,
%   title={Stability analysis of parabolic linear PDEs with two spatial dimensions using Lyapunov method and SOS},
%   author={Meyer, Evgeny and Peet, Matthew M},
%   booktitle={2015 54th IEEE Conference on Decision and Control (CDC)},
%   pages={1884--1890},
%   year={2015},
%   organization={IEEE}
% }
%
% % [16] - 
% @inproceedings{demetriou2019feedback,
% title={Feedback kernel approximations and sensor selection for controlled 2D parabolic PDEs using computational geometry methods},
% author={Demetriou, Michael A and Hu, Weiwei},
% booktitle={2019 IEEE 58th Conference on Decision and Control (CDC)},
% pages={2144--2150},
% year={2019},
% organization={IEEE}
% }
%
% % [17] -
% @article{antonelli2021linear,
%  title={Linear stability analysis of the homogeneous Couette flow in a 2D isentropic compressible fluid},
%  author={Antonelli, Paolo and Dolce, Michele and Marcati, Pierangelo},
%  journal={Annals of PDE},
%  volume={7},
%  number={2},
%  pages={1--53},
%  year={2021},
%  publisher={Springer}
%}
%
% % [18] -
% @article{krstic2009control,
%   title={Control of an unstable reaction--diffusion PDE with long input delay},
%   author={Krstic, Miroslav},
%   journal={Systems \& Control Letters},
%   volume={58},
%   number={10-11},
%   pages={773--782},
%   year={2009},
%   publisher={Elsevier}
% }
%
% % [19] -
% @article{kang2017boundary,
%   title={Boundary control of delayed ODE--heat cascade under actuator saturation},
%   author={Kang, Wen and Fridman, Emilia},
%   journal={Automatica},
%   volume={83},
%   pages={252--261},
%   year={2017},
%   publisher={Elsevier}
% }
%
% % [20] -
% @article{CALISKAN2009,
%  title = {Stability Analysis of the Heat Equation with Time-Delayed Feedback},
%  journal = {IFAC Proceedings Volumes},
%  volume = {42},
%  number = {6},
%  pages = {220-224},
%  year = {2009},
%  note = {6th IFAC Symposium on Robust Control Design},
%  issn = {1474-6670},
%  doi = {https://doi.org/10.3182/20090616-3-IL-2002.00038},
%  url = {https://www.sciencedirect.com/science/article/pii/S1474667015404057},
%  author = {Sina Yamaç çalişkan and Hitay özbay},
% }
%
% % [21] -
% @article{XU_2006,
%      author = {Xu, Gen Qi and Yung, Siu Pang and Li, Leong Kwan},
%      title = {Stabilization of wave systems with input delay in the boundary control},
%      journal = {ESAIM: Control, Optimisation and Calculus of Variations},
%      pages = {770--785},
%      publisher = {EDP-Sciences},
%      volume = {12},
%      number = {4},
%      year = {2006},
%      doi = {10.1051/cocv:2006021},
%      zbl = {1105.35016},
%      mrnumber = {2266817},
%      language = {en},
%      url = {http://www.numdam.org/articles/10.1051/cocv:2006021/}
% }


function [index,BATCH,TERM,GUI,params] = process_inputs(varargin0,nargin1)
% Subroutine to process the user inputs for the 
% examples_PDE_library_PIETOOLS function.

n_examples = 42;

BATCH = 0;      % if nonzero, batch-based PDE is assigned as output number BATCH of this function
TERM = 0;       % if nonzero, term-based PDE is assigned as output number TERM of this function
GUI = 0;        % 1 to load GUI, 0 not to

% Suppress (some) warnings in PDE initialization
evalin('base','silent_initialize_pde = true;');
% Each example comes with its own executives, so initially disable all
evalin('base','stability = 0;');
evalin('base','stability_dual = 0;');
evalin('base','Hinf_gain = 0;')
evalin('base','Hinf_gain_dual = 0;')
evalin('base','Hinf_estimator = 0;')
evalin('base','Hinf_control = 0;')

params = {};

% Collect the inputs
if nargin1==0 %<-- If there is no input, pause the script and let the user specify an example
    userinp = input('\n Select an example (1 through 40) to convert \n ---> ','s');
    varargin0 = split(userinp,[" ",","]);
    index = str2double(varargin0{1});
    if ~isnan(index) && (index>=0 && index<=n_examples)
        varargin0{1} = index;
    else
        userinp = input('\n No existing example specified... Please input an integer value 1 through 29 to extract the example \n ---> ','s');
        index = str2double(varargin0{1});
        if isnan(index)
            error('Please specify the desired example as first argument when calling ''examples_PDE_library_PIETOOLS''')
        else
            varargin0 = [str2double(userinp);varargin];
        end
    end
end
nargin0 = length(varargin0);

if nargin0==1 %<-- If there is one input, this input MUST correspond to the desired PDE
    if isdouble(varargin0{1}) && varargin0{1}<=n_examples
        index = varargin0{1};
        TERM = 1;
    elseif contains(varargin0{1},'batch','IgnoreCase',true)
        index = randi(n_examples,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        BATCH = 1;
    elseif contains(varargin0{1},'term','IgnoreCase',true)
        index = randi(n_examples,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        TERM = 1;
    elseif contains(varargin0{1},'gui','IgnoreCase',true)
        index = randi(n_examples,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        GUI = 1;
    else
        error('First argument must correspond to an example listed in this file')
    end
elseif nargin0>=2
    if isdouble(varargin0{1}) && varargin0{1}<=n_examples
        index = varargin0{1};
    else
        error('First argument must correspond to an example listed in this file')
    end
    pindx = 1;
    for j=2:nargin0
        if contains(varargin0{j},'batch','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==1))
            if BATCH==0
                BATCH = TERM+1;
            end
        elseif contains(varargin0{j},'term','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==2))
            if TERM==0
                TERM = BATCH+1;
            end
        elseif contains(varargin0{j},'gui','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==3))
            GUI = 1;
        elseif ischar(varargin0{j})
            try eval(varargin0{j});      %<-- In this case we assume the input defines certain parameters
                if contains(varargin0{j},';')
                    params{pindx} = varargin0{j};
                else
                    params{pindx} = [varargin0{j},';'];
                end
                pindx = pindx+1;
            catch
                disp(['Warning: Argument ',num2str(j),' is not understood, and is therefore ignored']);
            end
        else
            disp(['Warning: Argument ',num2str(j),' is not understood, and is therefore ignored']);
        end
    end
end
if BATCH==0 && TERM==0 && GUI==0
    disp('Warning: The input format has not been (appropriately) defined, reverting to terms input format...');
    TERM=1;
end

end