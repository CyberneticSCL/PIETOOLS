%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIETOOLS.m     PIETOOLS 2021b
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Table of Contents
%
%__________________________________________________________________________
%     | Problem                                     | Parameters            Notes / References
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
% 1.  | PDE: x_{t} = x_{s}                          |   
%     | BCs: x(s=0) = 0                             |
%-----|-------------------------------------------------------------------
% 2.  | PDE: x_{t} = Fm*x - Lm*x_{s}                | Different parameters          Lamare 2016 [1] 
%     | BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]     | may be invoked calling
%     |      [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]     | examples 2.1, 2.2, 
%     |                                             | 2.3, and 2.4.
%-------------------------------------------------------------------------
% 3.  | PDE: x_{t} = Mm*x - Lm*x_{s}                | Mm, Lm, K00, K01,             Diagne 2012 [2] 
%     | BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]     | K10 and K11 may be
%     |      [x_-(s=1)] = [K10, K11] [x_-(s=0)]     | specified
%-----|--------------------------------------------------------------------
% 4.  | PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}       | Different parameters          Saba 2019 [3] 
%     |      x2_{t} = sig2*x1 + (1/r2)*x2_{s}       | may be invoked calling  
%     | BCs: x1(s=0) = qb*x2(s=0)                   | examples 4.1, 4.2, 
%     |      x2(s=1) = pb*x1(s=1)                   | 4.3, 4.4, and 4.5.  
%
%--------------------------------------------------------------------------
%        Diffusion and Heat-Equation Type Systems
%--------------------------------------------------------------------------
% 5.  | PDE: x_{t} = lam*x + x_{ss}                 | lam = 9.86            (stable for lam < pi^2 = 9.8696)  
%     | BCs: x(s=0) = 0,      x(s=1) = 0            |                               Ahmadi 2015 [5] 
%--------------------------------------------------------------------------
% 6.  | PDE: x_{t} = lam*x + x_{ss}                 | lam = 2.466           (unstable for lam > pi^2/4 = 2.467)        
%     | BCs: x(s=0) = 0,      x_{s}(s=1) = 0        |                               Gahlawat 2017 [4]
%--------------------------------------------------------------------------
% 7.  | PDE:  x_{t} = a(s)*x_{ss}                   | a = s^3 - s^2 + 2     (unstable for lam > 4.66)            
%     |               + b(s)*x_{s}                  | b = 3*s^2 - 2*s               Gahlawat 2017 [4]
%     |               + c(s,lam)*x                  | c =-0.5*s^3 + 1.3*s^2 
%     | BCs:  x(s=0) = 0,     x_{s}(s=1) = 0        |    - 1.5*s + 0.7 +lam
%     |                                             | lam = 4.66                     
%--------------------------------------------------------------------------
% 8.  | ODE:  xo_{t} = A * xo + Bxr * x_{s}(s=a)    | a = 1    b = 1        (stable for lam < pi^2 = 9.8696)
%     | PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo  | lam = pi^2-1                  Das 2018 [7] (Example 2)
%     | BCs:  x(s=a) = 0,     x(s=b) = 0            | A, Bxr, and Bpv fixed
%--------------------------------------------------------------------------
% 9.  | ODE:  xo_{t} = k*xo                         | k = -1                (stable for k<0)
%     | PDE:  x_{t} = x_{ss}                        |
%     | BCs:  x(s=0) = 0,     x(s=1) = xo           |
%--------------------------------------------------------------------------
% 10. | ODE: xo_{t} = xo + x_{s}(s=0)               | k = -2                (stable for k=-2)            
%     | PDE: x_{t} = x_{ss}                         |                               Tang 2011 [13]
%     | BCs: x(s=0) = -xo,    x(s=1) = k*xo         |
%--------------------------------------------------------------------------
% 11. | PDE:  x_{t} = Cm*x + (1/R)*x_{ss}           | R = 2.93              (stable for R<2.7)
%     | BCs:  x(s=0) = 0,     x_{s}(s=1) = 0        | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6]
%--------------------------------------------------------------------------
% 12. | PDE:  x_{t} = Cm*x + (1/R)*x_{ss}           | R = (21+9)            (stable for R<21)
%     | BCs:  x(s=0) = 0,     x_{s}(s=1) = 0        | Cm = [0 0 0;                  Ahmadi 2015 [5] Adapted Example B
%     |                                             |       s 0 0;
%     |                                             |       -s^2 0 0]
%
%--------------------------------------------------------------------------
%        Beam Type Equations 
%--------------------------------------------------------------------------
% 13. | PDE: u_{tt} = -c*u_{ssss}                   | c = 0.1
%     | BCs: u(s=0) = 0,        u_{ss}(s=1) = 0     |                               Peet 2019 [8] (Example 8.1.0.1)
%     |      u_{s}(s=0) = 0,    u_{sss}(s=1) = 0    |
%     | Use states x1 = u_{t}, x2 = u_{ss}.         |
%     | =>                                          |
%     | PDE: x1_{t} = -c * x2_{ss}                  |
%     | BCs: u(s=0) = 0,        u_{ss}(s=1) = 0     |
%     |      u_{s}(s=0) = 0                         |
%     |                                             |
%--------------------------------------------------------------------------
%     | PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                           Peet 2019 [8] (Example 8.1.0.2)
%     |      r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%     | BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%     |      w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
%     |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% 14. | Use states:                                 | k = 1                 Hyperbolic implementation (stable)
%     |     x1 = w_{t},   x2 = k*aa*g * (w_{s}-phi),| aa = 1
%     |     x3 = phi_{t}, x4 = E*II * phi_{s}.      | II = 1
%     | =>                                          | g = 1
%     | PDE: x1_{t} = (1/r/aa) * x2_{s}             | E = 1
%     |      x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3 | r = 1
%     |      x3_{t} = (1/r/II)*x2 + (1/r/II)*x4_{s} | 
%     |      x4_{t} = E*II * x3_{s}                 |
%     | BCs: x1(0) = 0,         x2(1) = 0           |
%     |      x3(0) = 0,         x4(1) = 0           |
%     |                                             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% 15. | Assume all parameters are 1, and use states:|                       Hyperbolic/diffusive implementation (unstable)
%     |      x1 = w_{t},    x2 = w_{s},             |
%     |      x3 = phi_{t},  x4 = phi.               |
%     | =>                                          |
%     | PDE: x1_{t} = x2_{s} - x4_{s}               |
%     |      x2_{t} = x1_{s}                        |
%     |      x3_{t} = x2 - x4                       |
%     |      x4_{t} = x3                            |
%     | BCs: x4(0) = 0,     x4_{s}(1) = 0,          |
%     |      x3(0) = 0,     x1(0) = 0,              |
%     |      x2(1) - x4(1) = 0                      |
%     |                                             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% 16. | We rewrite the system as a single equation: |                       4th order implementation (stable)
%     | PDE: 0 = alp*w_{ssss} + bet*w_{tt}          |
%     |           - gam*w_{ttss} + w_{tttt}         |
%     | Use states:                                 |
%     |      x1 = w_{ttt},      x2 = w_{t},         |
%     |      x3 = w_{tt},       x4 = w.             |
%     | =>                                          |
%     | PDE: x1_{t} = -alp*x4_{ssss} - bet*x3       |
%     |               + gam*x3_{ss}                 |
%     |      x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2|
%     | BCs: w(s=0) = 0,        w_{s}(s=0) = 0      |
%     |      w_{ss}(s=1) - w(s=1) = 0               |
%     |      w_{sss}(s=1) - w_{s}(s=1) = 0          |
%     |      w_{t}(s=0) = 0     w_{ts}(s=0) = 0     |
%     |      w_{tt}(s=0) = 0    w_{tts}(s=0) = 0    |
%     
%--------------------------------------------------------------------------
%       Wave Equations
%--------------------------------------------------------------------------
% 17. | PDE: u_{tt} = u_{ss}                        | k=1                   (stable for k=1)                  
%     | BCs: u(s=0) = 0,                            |                               Peet 2019 [8] (Example 8.2)
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
%     | Use states x1 = u_{s}, x2 = u_{t}.          |
%     | =>                                          |
%     | PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}     |
%     | BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0 |
%     |                                             |
%-----|---------------------------------------------|----------------------
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1                         Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% 18. | Use states x1 = u_{t}, x2 = u.              | k = 1
%     | =>                                          | ad = 1
%     | PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}  |
%     | BCs: x1(0) = 0,     x2(0) = 0,              |
%     |      k*x1(1) + x2_{s}(1) = 0                |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% 19. | Use states x1 = u_{t}, x2 = u, x3 = u_{s}.  | k = 1                
%     | =>                                          | ad = 1
%     | PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}   |
%     | BCs: x1(0) = 0,     x2(0) = 0,              |
%     |      k*x1(1) + x3(1) = 0                    |
%
%==========================================================================
%       Hinf-gain analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic/Transport/Balance Type Systems
%--------------------------------------------------------------------------
% 20. | PDE: x_{t} = -x_{s} + w(t)                  |                       (gamma = 0.5)
%     | BCs: x(s=0) = 0                                |
%     | Out: z(t) = int(x(t,s),s,0,1)               |
%--------------------------------------------------------------------------
% 21. | PDE: phi_{tt} = phi_{ss} + w(t)             | k = 0.5               (gamma = 2 for k = 0.5)  
%     | BCs: phi(s=0) = 0                           |
%     |      phi_{s}(s=1) = -k*phi_{t}(s=1)         |
%     | Out: z(t) = int(phi_{t}(t,s),s,0,1)         |
%     | Use states x1 = phi_{s},    x2 = phi_{t}    |
%     | =>                                          |
%     | PDE: x1_{t} = x2_{s}                        |
%     |      x2_{t} = x1_{s} + w(t)                 |
%     | BCs: x2(0) = 0,     x1(1) + k*x2(1) = 0     |
%     |                                             |
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
% 22. | PDE: x_{t} = x_{ss} + s*w(t)                |                       (gamma = 0.3333)
%     | BCs: x(s=0) = 0,        x_{s}(s=1) = 0      |
%     | Out: z(t) = int(x(t,s),s,0,1)               |
%--------------------------------------------------------------------------
% 23. | PDE: x_{t} = A2(s)*x_{ss}                   | A2 = s^3 - s^2 + 2;   (gamma = 15.147 for lamb = 4.6)       
%     |              + A1(s)*x_{s}                  | A1 = 3*s^2 - 2*s;             Shivakumar 2019 [12] (Example 1)
%     |              + A0(s)*x + w(t)               | A0 =-0.5*s^3 +1.3*s^2
%     | BCs: x(s=0) = 0,        x_{s}(s=1) = 0      |     -1.5*s +0.7 +lam
%     | Out: z(t) = x(t,1)                          | lam = 4.6;  
%--------------------------------------------------------------------------
% 24. | PDE: xi_{t} = w(t) + lamb*xi        i=1:ne  | lam = (1-1e-2)*pi^2   (gamma = 8.1069 for lamb = (1-1e-2)*pi^2)         
%     |               + sum(xk_{ss},k=1,i)          | ne = 1 (state size)           Shivakumar 2019 [12] (Example 3)
%     | BCs: x(s=0) = 0,        x(s=1) = 0          |   
%     | Out: z(t) = int(x(t,s),s,0,1)               |
%--------------------------------------------------------------------------
% 25. | PDE: x_{t} = Cm*x + (1/R)*x_{ss} + s*w(t)   | R = 2.6               (gamma = 0.8102 for R = 2.6)
%     | BCs: x(s=0) = 0,        x(s=1) = 0          | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6]
%     | Out: z(t) = int(x1(t,s),s,0,1)              | 
%--------------------------------------------------------------------------
% 26. | PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + s*w(t)| R = (21-1e-3)         (gamma =  4.23, for R = 21-1e-3)
%     | BCs: x(s=0) = 0,        x(s=1) = 0          | Cm = [0,0,0;                  Shivakumar 2019 [12] (Example 2)
%     | Out: z(t) = int(x(t,s),s,0,1)               |       s,0,0;
%     |                                             |       s^2,-s^3,0]
%     
%==========================================================================
%       Hinf-optimal observer
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
% 27. | PDE: x_{t} = x_{s} + w(t)                   | ne = 1 (state size)   (Answer: 1.0012)
%     | BCs: x(s=1) = 0                             |
%     | Out: z(t) = int(x(t,s),s,0,1) + w(t)        |
%     |      y(t) = int(x(t,s),s,0,1)               |
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
% 28. | PDE: x_{t} = x_{ss} + w(t)                  | ne = 1 (state size)   (Answer: 1.0045) 
%     | BCs: x(s=0) = 0,        x(s=1) = 0          |
%     | Out: z(t) = int(x(t,s),s,0,1) + w(t)        |
%     |      y(t) = x_{s}(s=1)                      |
%
%==========================================================================
%       Hinf-optimal controller
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
% 29. | PDE: x_{t} = lam*x + x_{ss} + u(t)          | lam = 10;
%     | BCs: x(s=0) = 0,        x(s=1)=0            |
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
% 31. | PDE: x_{t} = c1*x_{s1} + c2*x_{s2}          | c1 = 1; c2 = 1;
%     | BCs: x(s1=0) = 0,   x(s2=0) = 0;            | ne = 1 (state size)
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
% 32. | PDE: x_{t} = lam*x + c1*x_{(2,0)}           | lam = 19;             (stable for lam <= 2*pi^2) 
%     |                         + c2*x_{(0,2)}      | c1 = 1;   c2 = 1;             Holmes, 1994 [14] Eq. (14)   
%     | BCs: x(s1=0) = 0,   x(s2=0) = 0,            | ne = 1 (state size)
%     |      x(s1=1) = 0,   x(s2=1) = 0;            |
%--------------------------------------------------------------------------
% 33. | PDE: x_{t} = C*(x_{(2,0} + x_{(0,2)})       | C = 1;
%     |               - b1*x_{(1,0)} - b2*x_{(0,1)} | b1 = 0.5; b2 = 2;             Holmes, 1994 [14] Eq. (2)
%     | BCs: x(s1=0) = 0,   x(s2=0) = 0,            | ne = 1 (state size)
%     |      x(s1=1) = 0,   x(s2=1) = 0;            |
%--------------------------------------------------------------------------
% 34. | PDE: x_{t} = -1/(2*lam) * x_{tt}            | lam = 1;
%     |               +(C^2/2*lam)*(x_{(2,0}        | C = 1;                        Holmes, 1994 [14] Eq. (3)
%     |                               + x_{(0,2)})  | ne = 1 (state size)
%     | BCs: x(s1=0) = 0,   x(s2=0) = 0,            | 
%     |      x(s1=1) = 0,   x(s2=1) = 0;            |
%     | Use states x1 = x,    x2 = x_{t}            |
%     | =>                                          |
%     | PDE: x1_{t} = x2                            |
%     |      x2_{t} = -2*lam*x2                     |
%     |              + C^2*(x1_{(2,0} + x1_{(0,2)}) |
%--------------------------------------------------------------------------
% 35. | PDE: x_{t} = a*x                            | a = 4.9;              (stable for a <= 0.5*pi^2) 
%     |               + b1*x_{(1,0)} + b2*x_{(0,1)} | b1 = 0;   b2 = 0;    
%     |                + c1*x_{(2,0)} + c2*x_{(0,2)}| c1 = 1;   c2 = 1;
%     | BCs: x(s1=0) = 0,   x_{s1}(s1=1) = 0,       | ne = 1 (state size)
%     |      x(s2=0) = 0,   x_{s2}(s2=1) = 0;       |
%--------------------------------------------------------------------------
% 36. | PDE: x_{tt} = c1*u_{(2,0)} + c2*x_{(0,2)}   | c1 = 1;   c2 = 1;
%     | BCs: x(s1=0) = 0;       x(s2=0) = 0;        | ne = 1; (state size)
%     |      x(s1=1) = 0;       x(s2=1) = 0;        |
%--------------------------------------------------------------------------
% 37. | PDE: p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}| M = 0.1;
%     |      v1_{t} = -s2*v1_{s1} -v2               | lam = 1; nu = 1;              Antonelli, 2021 [17]
%     |               -(1/M^2)*p_{s1}               |
%     |                 +nu*(v1_{s1s1} + v1_{s2s2}) |
%     |                 +lam*(v1_{s1s1} + v2_{s2s1})|
%     |       v2_{t} = -s2*v2_{s1} -(1/M^2)*p_{s2}  |
%     |                 +nu*(v2_{s1s1} + v2_{s2s2}) |
%     |                 +lam*(v1_{s1s2} + v2_{s2s2})|
%     | BCs: p(s1=0) = 0;      p(s2=0) = 0;         |
%     |      v(s1=0) = 0;      v(s2=0) = 0;         |
%     |      v(s1=1) = 0;      v(s2=1) = 0;         |
%--------------------------------------------------------------------------
% 38. | ODE: x1_{t} = (A+BK)*x1 + B*x2(s1=0,s2=0)   | A = I; B = I;
%     | PDE: x2_{t} = c1*x_{(2,0)} + c2*x_{(0,2)}   | K = -2*I;
%     | BCs: x2_{(1,0)}(s1=0) = 0;  x2(s1=1) = 0;   |
%     |      x2_{(0,1)}(s2=0) = 0;  x2(s2=1) = 0;   |
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The currently defined PDE struct will be cleared to avoid conflict
pvar s theta s1 s2 theta1 theta2

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);


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
if nargin==0 %<-- If there is no input, pause the script and let the user specify an example
    userinp = input('\n Select an example (1 through 29) to convert \n ---> ','s');
    varargin0 = split(userinp,[" ",","]);
    index = str2double(varargin0{1});
    if ~isnan(index) && (index>=0 && index<=29)
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
    nargin0 = length(varargin0);
else
    nargin0 = nargin;
    varargin0 = varargin;
end

if nargin0==1 %<-- If there is one input, this input MUST correspond to the desired PDE
    if isdouble(varargin0{1}) && varargin0{1}<=30
        index = varargin0{1};
        BATCH = 1;
    elseif contains(varargin0{1},'batch','IgnoreCase',true)
        index = randi(29,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        BATCH = 1;
    elseif contains(varargin0{1},'term','IgnoreCase',true)
        index = randi(29,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        TERM = 1;
    elseif contains(varargin0{1},'gui','IgnoreCase',true)
        index = randi(29,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        GUI = 1;
    else
        error('First argument must correspond to an example listed in this file')
    end
elseif nargin0>=2
    if isdouble(varargin0{1}) && varargin0{1}<=40
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
    disp('Warning: The input format has not been (appropriately) defined, reverting to batch input format...');
    BATCH=1;
end
npars = length(params);

fprintf([' --- Extracting ODE-PDE example ', num2str(index),' ---\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if index==0
    % % Example 2 from amritams paper
lambda = pi^2-1; %stable for lambda<pi^2
a=0; b= 1; PDE.dom = [a b];
PDE.nx = 4; PDE.n0=0; PDE.n1=0; PDE.n2 =2;
PDE.A = [-1.2142, 1.9649, 0.2232, 0.5616;
            -1.8042, -0.7260, -0.3479, 5.4355;
            -0.2898, 0.7381, -1.7606, 0.8294;
            -0.9417, -5.3399, -1.0704, -0.7590]; 
PDE.E0 = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0]*[zeros(PDE.n2) zeros(PDE.n2) eye(PDE.n2) zeros(PDE.n2)];
PDE.A2 = eye(PDE.n2); PDE.A0 = lambda*eye(PDE.n2); 
PDE.B = [eye(PDE.n2), zeros(PDE.n2), zeros(PDE.n2),zeros(PDE.n2); 
           zeros(PDE.n2), eye(PDE.n2), zeros(PDE.n2), zeros(PDE.n2)]; 
PDE.E = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
%PDE.Bx = zeros(4,4);
PDE_b=PDE;
elseif index==1
%%
% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% % % % %------------------------------------------------------------------
% % Transport equation x_{t} = x_{s}
% % with BC            x(s=0) = 0

 %%% Executive Function:
 evalin('base','stability = 1;');
 evalin('base','stability_dual = 1;');


if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];
 
 PDE_b.A1= 1; 
 
 on = eye(PDE_b.n1);   ze = zeros(PDE_b.n1);
 PDE_b.B = [ze on];

end
if TERM~=0
 %%% Term-based input format
 % Single 1D state component
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];

 % PDE: x_{t} = x_{s}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = 1;

 % BCs: 0 = x(0)
 PDE_t.BC{1}.term{1}.loc = 0;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Transport_Eq.mat'));
 logval = app.loadData(data);
 if logval
    disp("Failed to load data object. Incorrect structure");
 end
end


elseif index>=2 && index<3
% % % %====================================================================
%%%% Examples from Lamare [1]:
% % PDE          x_{t} = Fm*x - Lm*x_{s}
% % with BC [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]
% %         [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
% %
%%%% There are several examples from [1] included here. Add a decimal to
%%%% your example input to specify a particular set of parameters.

if index==2 || index==2.1
% % Example 5.1
 Gm1=[.2];
 Gm2=[-.3];
 Gm3=[.6];
 Gm4=[.1];
 Lm=[-3 0;0 1];
 Fm=[.2 -.3; .6 .1];
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Lamare_Ex_5_1.mat');
elseif index==2.2
% % Example 5.2
 Gm1=[.1];
 Gm2=[-.8];
 Gm3=[.6];
 Gm4=[-.4];
 Lm=[-1 0;0 1];
 Fm=[-.3 .1; .1 -.3];
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Lamare_Ex_5_2.mat');
elseif index==2.3
% % Example 5.3
 Gm1=[-.2202];
 Gm2=[1.3955];
 Gm3=[-.0596];
 Gm4=[.2090];
 Lm=[-1 0;0 2];
 Fm=[-.1 .1; .5 -.8];
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Lamare_Ex_5_3.mat');
elseif index==2.4
% % Example 5.4
 Gm1=[.5];
 Gm2=[-.4];
 Gm3=[.2];
 Gm4=[.8];
 Lm=[-2 0;0 1];
 Fm=[-.6565 -.3743; -.113 -.6485];
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Lamare_Ex_5_4.mat');
else
 error(['Example ',num2str(indx),'does not exist, please choose from 2.1, 2.2, 2.3 or 2.4']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Examples from [1]

 %%% Executive function
 evalin('base','stability = 1;')
%%% evalin('base','stability_dual = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2=0;
 PDE_b.dom = [0,1];
 PDE_b.A0 = Fm; 
 PDE_b.A1 = -Lm;
 ny1 = size(Gm1,1);   ny2 = size(Gm3,1);
 on0 = eye(ny1);   on1 = eye(ny2);   zer12 = zeros(ny1,ny2);
 PDE_b.B = [-Gm1 zer12 on0 -Gm2;
          -Gm3 on1 zer12' -Gm4];

end
if TERM~=0
 %%% Term-based input format
 % A single state component.
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = Fm*x
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = Fm;
 
 % PDE: x_{t} = ... -Lm*x_{s}
 PDE_t.x{1}.term{2}.x = 1;
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = -Lm;
 
 ny1 = size(Gm1,1);    ny2 = size(Gm3,1);
 on0 = eye(ny1);       on1 = eye(ny2);     zer12 = zeros(ny1,ny2);
 
 % BCs: 0 = [Gm1,0;-Gm3,1] * x(0)       
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [-Gm1, zer12; -Gm3, on1];
 
 % BCs: 0 = ... + [1,-Gm2;0,-Gm4] * x(1)
 PDE_t.BC{1}.term{2}.x = 1;
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [on0, -Gm2; zer12', -Gm4];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(dir);
 logval = app.loadData(data);
 if logval
    disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==3
% % % %====================================================================
% % Example from Diagne [2]
% % PDE          x_{t} = Mm*x - Lm*x_{s}
% % with BC [x_+(s=0)] = [K00, K01] [x_+(s=1)]
% %         [x_-(s=1)] = [K10, K11] [x_-(s=0)]
%
% The implementation is based on a Linearized Saint–Venant–Exner Model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NOTE: The values here have been chosen arbitrarily (feel free to adjust)

 Hs = 1;       Vs = 1;     Bs = 0.1;
 g = 10;       Cf = 1;

 Lm = diag(sort(eig([Vs,Hs,0;10,Vs,10;0,1*Vs^2,0])));
 lm1 = Lm(1,1);    lm2 = Lm(2,2);      lm3 = Lm(3,3);

 al1 = Cf*(3*Vs - 2*lm1) * Vs/Hs * lm1/((lm1-lm2)*(lm1-lm3));
 al2 = Cf*(3*Vs - 2*lm2) * Vs/Hs * lm2/((lm2-lm3)*(lm2-lm1));
 al3 = Cf*(3*Vs - 2*lm3) * Vs/Hs * lm3/((lm3-lm1)*(lm3-lm2));

 Mm = [al1*ones(3,1) , al2*ones(3,1) , al3*ones(3,1)];

 a21 = (lm1 - lm2) * (1+ ((lm1-Vs) * (lm2-Vs))/(g*Hs));
 a32 = (lm2 - lm3) * (1+ ((lm2-Vs) * (lm3-Vs))/(g*Hs));
 a13 = (lm3 - lm1) * (1+ ((lm3-Vs) * (lm1-Vs))/(g*Hs));
 c21 = (lm3/g) * (lm1 - lm2);
 c32 = (lm1/g) * (lm2 - lm3);
 c13 = (lm2/g) * (lm3 - lm1);

 k1 = 1;       k2 = 2;
 pi2 = (a21 - c21*k1)/(a32 - c32*k1);
 pi3 = (a13 - c13*k1)/(a32 - c32*k1);
 chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
 chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));

 K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
 K00 = K(1,1);   K01 = K(1,2:3);
 K10 = K(2:3,1); K11 = K(2:3,2:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%% Executive function
 evalin('base','stability = 1;')
%%% evalin('base','stability_dual = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

ne = size(Lm,1);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = ne;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = Mm; 
 PDE_b.A1 = -Lm;

 ny1 = size(K00,1);   ny2 = size(K10,1);
 on0 = eye(ny1);   on1 = eye(ny2);   zer12 = zeros(ny1,ny2);
 PDE_b.B = [ on0 -K01 -K00 zer12;
          zer12' -K11 -K10  on1;];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = Mm * x
 PDE_t.x{1}.term{1}.C = Mm;

 % PDE: x_{t} = ... -Lm * x_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = -Lm;
 
 ny1 = size(K00,1);    ny2 = size(K10,1);
 on0 = eye(ny1);       on1 = eye(ny2);     zer12 = zeros(ny1,ny2);

 % BCs: 0 = [1,-K01;0,-K11] * x(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [on0, -K01; zer12', -K11];      

 % BCs: 0 = ... + [-K00,0;-K10,1] * x(1)
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [-K00, zer12; -K10, on1];   

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Diagne_Ex.mat'));
 logval = app.loadData(data);
 if logval
    disp("Failed to load data object. Incorrect structure");
 end
end



elseif index>=4 && index<5
% % % %====================================================================
% % Example from Saba [3]
% % PDE        x1_{t} = sig1*x2 - (1/r1)*x1_{s} 
% %            x2_{t} = sig2*x1 + (1/r2)*x2_{s}
% % with BCs   x1(s=0) = qb*x2(s=0)
% %            x2(s=1) = pb*x1(s=1)
%
%%%% There are several examples from [3] included here. Add a decimal to
%%%% your example input to specify a particular set of parameters.
if index==4 || index==4.1
 r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5;  % Stable with stripped settings
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_1.mat');
 
elseif index==4.2
 r1=.5;   r2=1.1;   sig1=1;     sig2=-.1;     qb=1.2;   pb=0;   % Stable with stripped settings
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_2.mat');
 
elseif index==4.3
 r1=.5;   r2=1.1;   sig1=1;     sig2=1;       qb=1.2;   pb=-.7;
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_3.mat');
 
elseif index==4.4
 r1=.5;   r2=1.1;   sig1=1;     sig2=.663;    qb=1.2;   pb=0;   % max sig2=.663 
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_4.mat');
 
elseif index==4.5
 r1=.5;   r2=1.1;   sig1=1;     sig2=2.048;   qb=1.2;   pb=-.4; % max sig2=1.049 
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_5.mat');
 
else
 error(['Example ',num2str(indx),'does not exist, please choose from 4.1, 4.2, 4.3, 4.4 or 4.5']);
end

 %%% Executive function
 evalin('base','stability = 1;')
%%% evalin('base','stability_dual = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 sig1; sig2 0]; 
 PDE_b.A1 = [-1/r1 0; 0 1/r2];
 PDE_b.B = [1 -qb, 0, 0; 0, 0, -pb, 1];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = [0,sig1;sig2,0] * x
 PDE_t.x{1}.term{1}.C = [0, sig1; sig2, 0];

 % PDE: x_{t} = ... + [-1/r1,0;0,1/r2] * x_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = [-1/r1, 0; 0, 1/r2];
 
 % BCs: 0 = [1,-qb;0,0] * x(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [1 -qb; 0 0];

 % BCs: 0 = ... + [0,0;-pb,1] * x(1)
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [0 0; -pb 1];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(dir);
 logval = app.loadData(data);
 if logval
    disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==5
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusion and Heat-Equation Type Systems
% % % % %------------------------------------------------------------------
% % Scalable Diffusion Equation on [0,1] adapted from Ahmadi 2015 [5]
% % PDE        x_{t} = lam*x + x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = 0
%
% % Stable for lam < pi^2 = 9.8696 (this line should remain commented)

 ne = 1;   lam = 0; %9.86
 
 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end
 
if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
 PDE_b.dom = [0,1];

 PDE_b.A0 = lam*eye(ne);
 PDE_b.A2 = 1*eye(ne);

 PDE_b.B=[eye(ne) zeros(ne) zeros(ne) zeros(ne);
        zeros(ne) eye(ne) zeros(ne)  zeros(ne)];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

 % PDE: x_{t} = lam*x
 PDE_t.x{1}.term{1}.C = lam * eye(ne);

 % PDE: x_{t} = ... + x_{ss}
 PDE_t.x{1}.term{2}.D = 2;
 PDE_t.x{1}.term{2}.C = eye(ne);

 % BCs: 0 = x(s=0)
 PDE_t.BC{1}.term{1}.loc = 0;
 
 % BCs: 0 = x(s=1)
 PDE_t.BC{2}.term{1}.loc = 1;              

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Reaction_Diffusion_Eq_1.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==6
% % % %====================================================================
% % Scalable Diffusion Equation Problem 1 from Gahlawat_2017 [4]
% % PDE        x_{t} = lam*x + x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
%
% % Unstable for lam > pi^2/4 = 2.467 (this line should remain commented)

 ne = 1;   lam = 2.466;

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2=ne;
 PDE_b.dom = [0,1];

 PDE_b.A0 = lam*eye(ne);
 PDE_b.A2 = 1*eye(ne);

 PDE_b.B = [eye(ne) zeros(ne) zeros(ne) zeros(ne);
          zeros(ne)  zeros(ne) zeros(ne) eye(ne)];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = lam*x
 PDE_t.x{1}.term{1}.C = lam*eye(ne);

 % PDE: x_{t} = ... + x_{ss}
 PDE_t.x{1}.term{2}.D = 2;

 % BCs: 0 = x(s=0)
 PDE_t.BC{1}.term{1}.loc = 0;      

 % BCs: 0 = x_{s}(s=1)
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;       

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Reaction_Diffusion_Gahlawat.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==7
% % % %====================================================================
% % Diffusion Equation Problem 2 from Gahlawat 2017 [4]
% % PDE        x_{t} = c(s)*x + b(s)*x_{s} + a(s)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
%
% % and with   a(s) = s^3 - s^2 + 2
% %            b(s) = 3*s^2 - 2*s
% %            c(s) = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam
%
% % Unstable for lam > 4.66 (this line should remain commented)

 ne = 1; lam = 4.66;
 afun = s^3 - s^2 + 2;
 bfun = 3*s^2 - 2*s;
 cfun = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
 PDE_b.dom = [0,1];

 PDE_b.A2 = afun*eye(ne);
 PDE_b.A1 = bfun*eye(ne);
 PDE_b.A0 = cfun*eye(ne);
 
 PDE_b.B = [eye(ne), zeros(ne), zeros(ne), zeros(ne);
            zeros(ne), zeros(ne), zeros(ne), eye(ne)];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

 % PDE: x_{t} = cfun * x
 PDE_t.x{1}.term{1}.C = cfun*eye(ne);

 % PDE: x_{t} = ... + bfun * x_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = bfun*eye(ne);
 
 % PDE: x_{t} = ... + afun * x_{ss}
 PDE_t.x{1}.term{3}.D = 2;
 PDE_t.x{1}.term{3}.C = afun;
 
 % BCs: 0 = x(0)
 PDE_t.BC{1}.term{1}.loc = 0;  

 % BCs: 0 = x_{s}(1)
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Parabolic_Eq_Gahlawat.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==8
% % % %====================================================================
% % Heat Equation coupled with ODE at the boundary from Das [7] (Example 2)
% % ODE        xo_{t} = A * xo + Bxr * x_{s}(s=a)
% % PDE        x_{t} = lam * x + x_{ss} + Bpv * xo
% % with BCs   x(s=a) = 0
% %            x(s=b) = 0
%
% % Stable for lam<pi^2
    
 lam = pi^2-1;  a = 0;  b = 1;

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;   PDE_b.nx = 4;
 PDE_b.dom = [a,b];
 
 % ODE
 PDE_b.A = [-1.2142,  1.9649,  0.2232,  0.5616;
            -1.8042, -0.7260, -0.3479,  5.4355;
            -0.2898,  0.7381, -1.7606,  0.8294;
            -0.9417, -5.3399, -1.0704, -0.7590]; 
 PDE_b.E0 = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0] * [zeros(PDE_b.n2) zeros(PDE_b.n2) eye(PDE_b.n2) zeros(PDE_b.n2)];

 % PDE
 PDE_b.A0 = lam*eye(PDE_b.n2);   PDE_b.A2 = eye(PDE_b.n2);
 PDE_b.E = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
 
 % BCs
 PDE_b.B = [eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2),zeros(PDE_b.n2); 
          zeros(PDE_b.n2), eye(PDE_b.n2), zeros(PDE_b.n2), zeros(PDE_b.n2)]; 
end
if TERM~=0
 %%% Term-based input format
 % Initialize ODE state component.
 PDE_t.x{1}.vars = [];
 % Initialize PDE state component.
 PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [a,b];

 % ODE: xo_{t} = A * xo;
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = [-1.2142,  1.9649,  0.2232,  0.5616;
                         -1.8042, -0.7260, -0.3479,  5.4355;
                         -0.2898,  0.7381, -1.7606,  0.8294;
                         -0.9417, -5.3399, -1.0704, -0.7590]; 
            
 % ODE: xo_{t} = ... + Bxr * x_{s}(s=0)
 PDE_t.x{1}.term{2}.x = 2;
 PDE_t.x{1}.term{2}.loc = 0;
 PDE_t.x{1}.term{2}.C = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
 
 % PDE: x_{t} = lam * x
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.C = lam*eye(2);
 
 % PDE: x_{t} = ... + x_{ss}
 PDE_t.x{2}.term{2}.x = 2;
 PDE_t.x{2}.term{2}.D = 2;
 
 % PDE: x_{t} = ... + Bpv * v
 PDE_t.x{2}.term{3}.x = 1;
 PDE_t.x{2}.term{3}.C = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
 
 % BCs: 0 = x(0)   
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = a;
 PDE_t.BC{1}.term{1}.C = eye(2);
 
 % BCs: 0 = x(1)
 PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{2}.term{1}.loc = b;
 PDE_t.BC{2}.term{1}.C = eye(2);
 
end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Heat_Eq_with_ODE_Das.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==9
% % % %====================================================================
% % Heat Equation coupled with ODE at the boundary
% % ODE        xo_{t} = k * xo
% % PDE        x_{t} = x_{ss} 
% % with BCs   x(s=0) = 0
% %            x(s=1) = xo
%
% % Stable for k<0

 k = -1;

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;   PDE_b.nx = 1;
 PDE_b.dom = [0,1];

 PDE_b.A2 = 1;  PDE_b.A = k;

 PDE_b.B = [1 0 0 0; 0 1 0 0]; PDE_b.Bx = [0;1];

end
if TERM~=0
 %%% Term-based input format
 % Initialize ODE state component.
 PDE_t.x{1}.vars = [];
 % Initialize 1D PDE state component.
 PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];

 % ODE: xo_{t} = k * xo;
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = k;

 % PDE: x_{t} = x_{ss}
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.D = 2;

 % BCs: 0 = x(0)
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = 0;

 % BCs: 0 = x(1) - xo
 PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{2}.term{1}.loc = 1;
 
 PDE_t.BC{2}.term{2}.x = 1;
 PDE_t.BC{2}.term{2}.C = -1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Heat_Eq_with_ODE_in_BC.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==10
% % % %====================================================================
% % % Coupled ODE-PDE system from [13]
% % ODE             xo_{t} = xo + x_{s}(s=0)
% % PDE             x_{t} = x_{ss}
% % With BCs        x(s=0) = -xo
% %                 x(s=1) = k * xo


 k = -2;

 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.nx = 1;   PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nu = 0;
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
 PDE_b.dom = [0,1];

 PDE_b.A = 1;   PDE_b.E0 = [0 0 1 0];
 PDE_b.A2 = 1;
 PDE_b.Bx = [-1; k];
 
 PDE_b.B = [1 0 0 0;
            0 1 0 0];

end
if TERM~=0
 %%% Term-based input format
 % Initialize ODE state component.
 PDE_t.x{1}.vars = [];
 % Initialize 1D PDE state component.
 PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];

 % ODE: xo_{t} = xo
 PDE_t.x{1}.term{1}.x = 1;
 
 % ODE: xo_{t} = ... + x_{s}(s=0)
 PDE_t.x{1}.term{2}.x = 2;
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.loc = 0;

 % PDE: x_{t} = x_{ss}
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.D = 2;

 % BC 1: 0 = x(0)      
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = 0;
 
 % BC 1: 0 = ... + xo   
 PDE_t.BC{1}.term{2}.x = 1;
 
 % BC 2: 0 = x(1) 
 PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{2}.term{1}.loc = 1;
 
 % BC 2: 0 = ... -k * xo
 PDE_t.BC{2}.term{2}.x = 1;
 PDE_t.BC{2}.term{2}.C = -k;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Heat_Eq_with_ODE_Tang.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==11
% % % %====================================================================
% % Example Diffusion-reaction from Ahmadi [6]
% % PDE        x_{t} = Cm*x + (1/R)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
%
 Cm = [1, 1.5; 5, 0.2];
 % Stable for R<2.7 (this line should remain commented)
 R = 2.93; % [2.93 2.94]

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 0;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 ne = size(Cm,2);
 on = eye(ne);      ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
 PDE_b.dom = [0,1];

 PDE_b.A0 = Cm;
 PDE_b.A2 = (1/R)*on;

 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = Cm * x
 PDE_t.x{1}.term{1}.C = Cm;

 % PDE: x_{t} = ... + (1/R) * x_{ss}
 PDE_t.x{1}.term{2}.D = 2;
 PDE_t.x{1}.term{2}.C = (1/R)*eye(ne);

 % BCs: 0 = x(0)
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BCs: 0 = x_{s}(1)
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Reaction_Diffusion_Ahmadi.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==12
% % % %====================================================================
% % Adapted Example B (Diffusion-Reaction equation) from [5]
% % PDE        x_{t} = Cm*x + (1/R)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
%
 Cm = [0 0 0; s 0 0; s -s^2 0];
 % Stable for R<21 (this line should remain commented)
 R = (21+9);

 %%% Executive function
 evalin('base','stability = 1;');
 evalin('base','stability_dual = 0;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 ne = size(Cm,2);  
 on = eye(ne);      ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
 PDE_b.dom = [0,1];

 PDE_b.A0 = Cm;  
 PDE_b.A2 = (1/R)*on;

 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = Cm * x
 PDE_t.x{1}.term{1}.C = Cm;

 % PDE: x_{t} = ... + (1/R) * x_{ss}
 PDE_t.x{1}.term{2}.D = 2;
 PDE_t.x{1}.term{2}.C = (1/R)*eye(ne);

 % BCs: 0 = x(0)
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BCs: 0 = x_{s}(1)
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;   

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Reaction_Diffusion_Ahmadi_2.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==13
% % % % %------------------------------------------------------------------
% % % % % 3. Beam Type Equations
% % % % %------------------------------------------------------------------
% % Euler-Bernoulli beam equation [8] (Example 8.1.0.1)
% % PDE        u_{tt} = -c*u_{ssss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=0) = 0 
% %            u_{ss}(s=1) = 0 
% %            u_{sss}(s=1) = 0 
%
% % We use states x1 = u_{t}, x2 = u_{ss}, so that our system becomes
% % PDE        x1_{t} = -c * x2_{ss}
% %            x2_{t} = x1_{ss}
% % with BCs   x1(0) = 0, x2(1) = 0, x1_{s}(0) = 0, x2_{s}(1) = 0
%
 c = 0.1;  %.01;   %c=EI/mu

 %%% Executive function
 evalin('base','stability = 1;')
 evalin('base','stability_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;
 PDE_b.dom = [0,1];

 PDE_b.A2 = [0 -c; 1 0];

 PDE_b.B = [ 1 0 0 0 0 0 0 0;
             0 0 0 1 0 0 0 0;
             0 0 0 0 1 0 0 0;
             0 0 0 0 0 0 0 1];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = [0,-c;1,0] * x_{ss}
 PDE_t.x{1}.term{1}.D = 2;
 PDE_t.x{1}.term{1}.C = [0, -c; 1, 0];

 % BCs: 0 = x1(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [1,0];

 % BCs: 0 = x2(1)
 PDE_t.BC{2}.term{1}.loc = 1;
 PDE_t.BC{2}.term{1}.C = [0,1];      

 % BCs: 0 = x1_{s}(0)
 PDE_t.BC{3}.term{1}.D = 1;
 PDE_t.BC{3}.term{1}.loc = 0;
 PDE_t.BC{3}.term{1}.C = [1,0];

 % BCs: 0 = x2_{s}(1)
 PDE_t.BC{4}.term{1}.D = 1;
 PDE_t.BC{4}.term{1}.loc = 1;
 PDE_t.BC{4}.term{1}.C = [0,1];
 
%  % Note that BCs can also be coupled together:
%  
%  PDE_t.BC{1}.term{1}.loc = 0;
%  PDE_t.BC{1}.term{1}.C = [1,0;0,0];
%  
%  PDE_t.BC{1}.term{2}.loc = 1;
%  PDE_t.BC{1}.term{1}.C = [0,0;0,1];
%  
%  PDE_t.BC{2}.term{1}.D = 1;
%  PDE_t.BC{2}.term{1}.loc = 0;
%  PDE_t.BC{2}.term{1}.C = [1,0;0,0];
%  
%  PDE_t.BC{2}.term{2}.D = 1;
%  PDE_t.BC{2}.term{2}.loc = 1;
%  PDE_t.BC{2}.term{1}.C = [0,0;0,1];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Euler_Bernoulli_Beam_Eq.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==14
% % % %====================================================================
% % Timoschenko beam equation (hyperbolic) [8] (Example 8.1.0.2)
% % PDE        r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})
% %            r*II * phi_{tt} = E*II * phi_{ss} + k*aa*g * (w_{s} - phi)
% % with BCs   phi(s=0) = 0
% %            w(s=0) = 0 
% %            phi_{s}(s=1) = 0 
% %            w_{s}(s=1) - phi(s=1) = 0 
%
% % We use states x1 = w_{t}, x2 = k*aa*g * (w_{s}-phi), x3 = phi_{t}, x4 = E*II * phi_{s}.
% % Then, our system becomes:
% % PDE        x1_{t} = (1/r/aa) * x2_{s}
% %            x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3
% %            x3_{t} = (1/r/II) * x2 + (1/r/II) * x4_{s}
% %            x4_{t} = E*II * x3_{s}
% % with BCs   x1(0) = 0, x3(0) = 0, x4(1) = 0, x2(1) = 0

 k = 1;    aa = 1;    II = 1;    g = 1;    E = 1;   r = 1;  %5-30kPa
 

 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 4;   PDE_b.n2 = 0; 
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 0 0 0; 
             0 0 -k*aa*g 0;
             0 1/r/II 0 0;
             0 0 0 0];
 PDE_b.A1 = [0 1/r/aa 0 0; 
            k*aa*g 0 0 0;
            0 0 0 1/r/II;
            0 0 E*II 0];

 PDE_b.B = [ 1 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0;
            0 0 0 0 0 0 0 1;
            0 0 0 0 0 1 0 0];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];
 
 % PDE: x2_{t} = -k*aa*g * x3,   x3_{t} = (1/r/II) * x2
 PDE_t.x{1}.term{1}.C = [0, 0,      0,       0; 
                         0, 0,      -k*aa*g, 0;
                         0, 1/r/II, 0,       0;
                         0, 0,      0,       0];

 % PDE: x1_{t} = (1/r/aa) * x2_{s},        x2_{t} = ... + k*aa*g * x1_{s}
 %      x3_{t} = ... + (1/r/II) * x4_{s}   x4_{t} = E*II * x3_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = [0,      1/r/aa, 0,    0; 
                         k*aa*g, 0,      0,    0;
                         0,      0,      0,    1/r/II;
                         0,      0,      E*II, 0];

 % BCs: 0 = x1(0),   0 = x3(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [1,0,0,0;0,0,1,0;0,0,0,0;0,0,0,0];

 % BCs: 0 = x4(1),   0 = x2(1)
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [0,0,0,0;0,0,0,0;0,0,0,1;0,1,0,0];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Timoshenko_Beam_1.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==15
% % % %--------------------------------------------------------------------
% % Timoschenko Beam equation Example (hyperbolic/diffusive) - unstable [8]
%
% % We now use states x1 = w_{t}, x2 = w_{s}, x3 = phi_{t}, x4 = phi.
% % We also assume all coefficients are 1. Then, our system becomes:
% % PDE        x1_{t} = x2_{s} - x4_{s}
% %            x2_{t} = x1_{s}
% %            x3_{t} = x2 - x4
% %            x4_{t} = x3
% % with BCs   x4(0) = 0, x4_{s}(1) = 0, x3(0) = 0, x1(0) = 0, x2(1)-x4(1) = 0

 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 % Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 3;   PDE_b.n2=1;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 0 0 0;
             0 0 0 0;
             0 1 0 -1;
             0 0 1 0];
 PDE_b.A1 = [0 1 0 -1;
             1 0 0 0;
             0 0 0 0;
             0 0 0 0];
 PDE_b.A2 = [0;0;1;0];

 PDE_b.B = [1 0 0 0 0 0 0 0 0 0;
            0 0 1 0 0 0 0 0 0 0;
            0 0 0 0 0 0 1 0 0 0;
            0 0 0 0 0 0 0 0 0 1;
            0 0 0 0 1 0 0 -1 0 0];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.vars = s;    PDE_t.dom = [0,1];
 
 % PDE 1: x1_{t} = x2_{s} - x4_{s}
 PDE_t.x{1}.term{1}.x = [2;4];
 PDE_t.x{1}.term{1}.D = 1;
 PDE_t.x{1}.term{1}.C = [1,-1];
 
 % PDE 2: x2_{t} = x1_{s}
 PDE_t.x{2}.term{1}.x = 1;
 PDE_t.x{2}.term{1}.D = 1;
 
 % PDE 3: x3_{t} = x2 - x4
 PDE_t.x{3}.term{1}.x = [2;4];
 PDE_t.x{3}.term{1}.C = [1,-1];
 
 % PDE 4: x4_{t} = x3
 PDE_t.x{4}.term{1}.x = 3;
 
 % BC 1: x4(0) = 0
 PDE_t.BC{1}.term{1}.x = 4;
 PDE_t.BC{1}.term{1}.loc = 0;
 
 % BC 2: x4_{s}(1) = 0 
 PDE_t.BC{2}.term{1}.x = 4;
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;
 
 % BC 3: x3(0) = 0 
 PDE_t.BC{3}.term{1}.x = 3;
 PDE_t.BC{3}.term{1}.loc = 0;
 
 % BC 4: x1(0) = 0
 PDE_t.BC{4}.term{1}.x = 1;
 PDE_t.BC{4}.term{1}.loc = 0;
 
 % BC 5: x2(1)-x4(1) = 0
 PDE_t.BC{5}.term{1}.x = [2;4];
 PDE_t.BC{5}.term{1}.loc = 1;
 PDE_t.BC{5}.term{1}.C = [1,-1];
 
% % % Alternatively, group together x1 x2 and x3
%  PDE_t.vars = s;    PDE_t.dom = [0,1];
%  
%  % PDE: x3_{t} = x2
%  PDE_t.x{1}.term{1}.x = 1;
%  PDE_t.x{1}.term{1}.C = [0 0 0; 0 0 0; 0 1 0];
% 
%  % PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}
%  PDE_t.x{1}.term{2}.x = 1;
%  PDE_t.x{1}.term{2}.D = 1;
%  PDE_t.x{1}.term{2}.C = [0 1 0; 1 0 0; 0 0 0];
%  
%  % PDE: x3_{t} = ... - x4
%  PDE_t.x{1}.term{3}.x = 2;
%  PDE_t.x{1}.term{3}.C = [0; 0; -1];
%  
%  % PDE: x1_{t} = ... - x4_{s}
%  PDE_t.x{1}.term{4}.x = 2;
%  PDE_t.x{1}.term{4}.D = 1;
%  PDE_t.x{1}.term{4}.C = [-1; 0; 0];
%  
%  % PDE: x3_{t} = ... + x4_{ss}
%  PDE_t.x{1}.term{5}.x = 2;
%  PDE_t.x{1}.term{5}.D = 2;
%  PDE_t.x{1}.term{5}.C = [0; 0; 1];
% 
%  % PDE: x4_{t} = x3
%  PDE_t.x{2}.term{1}.x = 1;
%  PDE_t.x{2}.term{1}.C = [0 0 1];
% 
%  % BC 1: 0 = x1(0) 
%  %       0 = x3(0)
%  PDE_t.BC{1}.term{1}.x = 1;
%  PDE_t.BC{1}.term{1}.loc = 0;
%  PDE_t.BC{1}.term{1}.C = [1,0,0;0,0,1];
% 
%  % BC 2: 0 = x4(0)
%  PDE_t.BC{2}.term{1}.x = 2;
%  PDE_t.BC{2}.term{1}.loc = 0;
% 
%  % BC 3: 0 = x4_{s}(1)
%  PDE_t.BC{3}.term{1}.x = 2;
%  PDE_t.BC{3}.term{1}.D = 1;
%  PDE_t.BC{3}.term{1}.loc = 1;
% 
%  % BC 4: 0 = x2(1)
%  PDE_t.BC{4}.term{1}.x = 1;
%  PDE_t.BC{4}.term{1}.loc = 1;
%  PDE_t.BC{4}.term{1}.C = [0,1,0];
% 
%  % BC 4: 0 = ... - x4(1)
%  PDE_t.BC{4}.term{2}.x = 2;
%  PDE_t.BC{4}.term{2}.loc = 1;
%  PDE_t.BC{4}.term{2}.C = -1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Timoshenko_Beam_2.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==16
% % % %--------------------------------------------------------------------
% % Timoschenko Beam equation 4th order implementation - stable
%
% % Assuming sufficiently smooth solutins, we may rewrite the system as a single equation
% % PDE        0 = alp*w_{ssss} + bet*w_{tt} - gam*w_{ttss} + w_{tttt}
%
% % We use states x1 = w_{ttt}, x2 = w_{t}, x3 = w_{tt}, x4 = w,
% % and consider BCs
% % w(s=0) = w_{s}(s=0) = 0
% % w_{ss}(s=1) - w(s=1) = w_{sss}(s=1) - w_{s}(s=1) = 0
% % w_{t}(s=0) = 0      w_{ts}(s=0) = 0
% % w_{tt}(s=0) = 0     w_{tts}(s=0) = 0
 
alp = 1;     bet = 1;     gam = 1;

% %%% Executive function
 evalin('base','stability = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
% %%% NOTE: Derivatives of order greater than 2 cannot be represented using 
% %%% the batch input format. They can be specified using the GUI.
 disp('Warning: No batch input formulation for this problem exists; reverting to term-based input format')
 BATCH = 0;
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.vars = s;    PDE_t.dom = [0,1];
 
 % PDE: x1 = -bet * x3 + gam * x3_{ss}
 PDE_t.x{1}.term{1}.x = 3;
 PDE_t.x{1}.term{1}.D = [0;2];
 PDE_t.x{1}.term{1}.C = [-bet, gam];

 % PDE: x1 = ... - alp * x4_{ssss}
 PDE_t.x{1}.term{2}.x = 4;
 PDE_t.x{1}.term{2}.D = 4;
 PDE_t.x{1}.term{2}.C = -alp;

 % PDE: x2 = x3
 PDE_t.x{2}.term{1}.x = 3;
 
 % PDE: x3 = x1
 PDE_t.x{3}.term{1}.x = 1;

 % PDE: x4 = x2
 PDE_t.x{4}.term{1}.x = 2;

 % BC 1: 0 = x4(0)
 PDE_t.BC{1}.term{1}.x = 4;
 PDE_t.BC{1}.term{1}.loc = 0;

 % BC 2: 0 = x4(1)
 PDE_t.BC{2}.term{1}.x = 4;
 PDE_t.BC{2}.term{1}.loc = 1;

 % BC 3: 0 = x4_{ss}(1) - x4(1)
 PDE_t.BC{3}.term{1}.x = 4;
 PDE_t.BC{3}.term{1}.D = [2;0];
 PDE_t.BC{3}.term{1}.loc = 1;
 PDE_t.BC{3}.term{1}.C = [1,-1];
      
 % BC 4: 0 = x4_{sss}(1) - x4_{s}(1)
 PDE_t.BC{4}.term{1}.x = 4;
 PDE_t.BC{4}.term{1}.D = [3;1];
 PDE_t.BC{4}.term{1}.loc = 1;
 PDE_t.BC{4}.term{1}.C = [1,-1]; 

 % BC 5: 0 = x2(0),             % BC 7: 0 = x3(0)      
 PDE_t.BC{5}.term{1}.x = 2;     PDE_t.BC{7}.term{1}.x = 3;
 PDE_t.BC{5}.term{1}.loc = 0;   PDE_t.BC{7}.term{1}.loc = 0;
 
 % BC 6: 0 = x2_{s}(0),         % BC 8: 0 = x3_{s}(0)
 PDE_t.BC{6}.term{1}.x = 2;     PDE_t.BC{8}.term{1}.x = 3;
 PDE_t.BC{6}.term{1}.D = 1;     PDE_t.BC{8}.term{1}.D = 1;
 PDE_t.BC{6}.term{1}.loc = 0;   PDE_t.BC{8}.term{1}.loc = 0;                 

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Timoshenko_Beam_3.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==17
% % % % %------------------------------------------------------------------
% % % % % 4. Wave Equations
% % % % %------------------------------------------------------------------
% % % %--------------------------------------------------------------------
% % Boundary-Damped Wave equation (Hyperbolic) [8] (Example 8.2)
% % PDE        u_{tt} = u_{ss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
%
% % We use states x1 = u_{s}, x2 = u_{t}.
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.

 k = 1;

 %%% Executive function
 evalin('base','stability = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end
 
if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 0; 0 0];
 PDE_b.A1 = [0 1; 1 0];

 PDE_b.B = [0 1 0 0;
          0 0 1 k];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.vars = s;    PDE_t.dom = [0,1];
 
 % PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}
 PDE_t.x{1}.term{1}.x = 2;
 PDE_t.x{1}.term{1}.D = 1;
 
 PDE_t.x{2}.term{1}.x = 1;
 PDE_t.x{2}.term{1}.D = 1;

 % BC1: 0 = x2(0)
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = 0;       

 % BC2: 0 = x1(1) + k*x2(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;
 PDE_t.BC{2}.term{2}.x = 2;
 PDE_t.BC{2}.term{2}.loc = 1;
 PDE_t.BC{2}.term{2}.C = k;  

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Wave_Eq_Boundary_Damped.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==18
% % % %====================================================================
% % Test 7.5c - Datko Boundary-Damped Wave equation [9]
% % PDE        u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} 
% % with BCs   u(s=0) = 0
% %            u_{s}(s=1) = -k*u_{t}(s=1)
%
% % We use states x1 = u_{t}, x2 = u.
% % Then x1(0) = 0, x2(0) = 0,  k*x1(1) + x2_{s}(1) = 0.

 k = 1;    ad = 1;

 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 1;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [-2*ad, -ad^2; 1, 0];
 PDE_b.A2 = [1; 0];

 PDE_b.B = [0 0 1 0 0 0;
          1 0 0 0 0 0;
          0 k 0 0 0 1];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.vars = s;        PDE_t.dom = [0,1];
 
 % PDE: x1_{t} = -2*ad * x1 - ad^2 * x2
 PDE_t.x{1}.term{1}.x = [1;2];
 PDE_t.x{1}.term{1}.C = [-2*ad, -ad^2];
 
 % PDE: x1_{t} = ... + x2_{ss}
 PDE_t.PDE.A{4}.coeff = 1; PDE_t.PDE.A{4}.Lstate = 1; PDE_t.PDE.A{4}.Rstate = 2; PDE_t.PDE.A{4}.D = 2; 
 PDE_t.x{1}.term{2}.x = 2;
 PDE_t.x{1}.term{2}.D = 2;

 % PDE: x2_{t} = x1
 PDE_t.x{2}.term{1}.x = 1; 

 % BC 1: 0 = x1(0)
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;

 % BC 2: 0 = x2(0)      
 PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{2}.term{1}.loc = 0;

 % BC 3: 0 = k * x1(1) + x2_{s}(1)
 PDE_t.BC{3}.term{1}.x = [1;2];
 PDE_t.BC{3}.term{1}.loc = 1;
 PDE_t.BC{3}.term{1}.D = [0;1];
 PDE_t.BC{3}.term{1}.C = [k, 1];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Wave_Eq_Datko_Boundary_Damped_1.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==19
% % % %--------------------------------------------------------------------
% % Test 7.5d - Datko Boundary-Damped Wave equation (Hyperbolic) [9]
%
% % Now we use states x1 = u_{t}, x2 = u, x3 = u_{s}.
% % Then x1(0) = 0, x2(0) = 0,  k*x1(1) + x3(1) = 0.

 k = 1;    ad = 1;

 %%% Executive function
 evalin('base','stability = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 3;   PDE_b.n2=0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [-2*ad -ad^2 0; 1 0 0;0 0 0];
 PDE_b.A1 = [0 0 1; 0 0 0;1 0 0];

 PDE_b.B = [1 0 0 0 0 0;
          0 1 0 0 0 0;
          0 0 0 k 0 1];

end
if TERM~=0
 %%% Term-based input format
 PDE_t.vars = s;    PDE_t.dom = [0,1];
 
 % PDE: x1_{t} = -2*ad * x1 - ad^2 * x2
 PDE_t.x{1}.term{1}.x = [1; 2];
 PDE_t.x{1}.term{1}.C = [-2*ad, -ad^2];
 
 % PDE: x1_{t} = ... + x3_{s}
 PDE_t.x{1}.term{2}.x = 3;
 PDE_t.x{1}.term{2}.D = 1;
 
 % PDE: x2_{t} = x1
 PDE_t.x{2}.term{1}.x = 1;

 % PDE: x3_{t} = x1_{s}
 PDE_t.x{3}.term{1}.x = 1;
 PDE_t.x{3}.term{1}.D = 1;

 % BC 1: 0 = x1(0)              BC 2: 0 = x2(0) 
 PDE_t.BC{1}.term{1}.x = 1;     PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = 0;   PDE_t.BC{2}.term{1}.loc = 0;

 % BC 3: 0 = k * x1(1) + x3(1)
 PDE_t.BC{3}.term{1}.x = [1; 3];
 PDE_t.BC{3}.term{1}.loc = 1;
 PDE_t.BC{3}.term{1}.C = [k,1];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Stability_Examples/Wave_Eq_Datko_Boundary_Damped_2.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf gain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic/Transport/Balance Type Systems
% % % % %------------------------------------------------------------------
% % Example pure transport equation 1D: 
% % PDE         x_{t} = -x_{s} + w(t)
% % With BC     x(s=0) = 0
% % And output  z(t) = int(x(t,s),s,0,1)
%
% % gamma = 0.5 (this line should remain commented)

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 evalin('base','Hinf_gain_dual = 1;')
 
 ne = 1;
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 on = ones(ne);     ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = ne;   PDE_b.n2 = 0;   PDE_b.nw = ne;   PDE_b.nz = ne;
 PDE_b.dom = [0,1];

 PDE_b.A1 = -1*on; 
 PDE_b.B21 = on;   PDE_b.B21 = on;
 PDE_b.Ca1 = 1;

 on = eye(PDE_b.n1); ze = zeros(PDE_b.n1);
 PDE_b.B = [on ze];

end
if TERM~=0
 %%% Term-based input format
 % Initialize a PDE state.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];
 
 % PDE: x_{t} = -x_{s}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = 1;
 PDE_t.x{1}.term{1}.C = -eye(ne);
 
 % PDE: x_{t} = ... + w
 PDE_t.x{1}.term{2}.w = 1;
 
 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BCs: 0 = x(0)        
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Transport_Eq_with_Disturbance.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==21
% % % %====================================================================
% % Example tip damped wave equation:
% % PDE         phi_{tt} = phi_{ss} + w(t)
% % With BCs    phi(s=0) = 0
% %             phi_{s}(s=1) = -k*phi_{t}(s=1)
% % And output  z(t) = int(phi_{t}(t,s),s,0,1)
%
% % We use states x1 = phi_{s}, x2 = phi_{t}
% % Then x2(0) = 0,  x1(1) + k*x2(1) = 0.
%
% % gamma = 2 for k = 0.5 (this line should remain commented)
 k = 0.5;

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 evalin('base','Hinf_gain_dual = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;   PDE_b.nw = 1;   PDE_b.nz = 1;
 PDE_b.dom = [0,1];
 np = PDE_b.n0 + PDE_b.n1 + PDE_b.n2;

 PDE_b.A1= [0 1; 1 0];
 PDE_b.B21 = [0; 1];
 PDE_b.Ca1 = [0 1];

 PDE_b.B = [0 1 0 0; 0 0 1 k];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state components.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 PDE_t.x{2}.vars = s;   PDE_t.x{2}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];
 
 % PDE: x1_{t} = x2_{s}
 PDE_t.x{1}.term{1}.x = 2;
 PDE_t.x{1}.term{1}.D = 1;

 % PDE: x2_{t} = x1_{s}
 PDE_t.x{2}.term{1}.x = 1;
 PDE_t.x{2}.term{1}.D = 1;
 
 % PDE: x2_{t} = ... + w
 PDE_t.x{2}.term{2}.w = 1;
 
 % Output: z = int_{0}^{1} x2(s) ds
 PDE_t.z{1}.term{1}.x = 2;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{2}.dom;

 % BC 1: 0 = x2(0)         
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.loc = 0;        

 % BC 2: 0 = x1(1) + k*x2(1)
 PDE_t.BC{2}.term{1}.x = [1; 2];
 PDE_t.BC{2}.term{1}.loc = 1;
 PDE_t.BC{2}.term{1}.C = [1, k];

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Wave_Eq_Tip_Damped.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==22
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% % Heat equation with distributed disturbance:
% % PDE         x_{t} = x_{ss} + s*w(t)
% % With BCs    x(s=0) = 0
% %             x_{s}(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
%
% % gamma = 0.3333  (this line should remain commented)

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 
 ne = 1;
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end
on = eye(ne);     ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = ne;   PDE_b.nz = ne;
 PDE_b.dom = [0,1];

 PDE_b.A2 = on;    PDE_b.B21 = s*on;
 PDE_b.Ca1 = on;

 PDE_b.B = [on ze ze ze;
            ze ze ze on];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];
 
 % PDE: x_{t} = x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = 2;
 
 % PDE: x_{t} = ... + s * w,
 PDE_t.x{1}.term{2}.w = 1;
 PDE_t.x{1}.term{2}.C = s*eye(ne);

 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BC 2: 0 = x_{s}(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Heat_Eq_with_Distributed_Disturbance.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==23
% % % %====================================================================
% %  Parabolic PDE example from [12] (Example 1)
% % PDE         x_{t} = A0(s)*x + A1(s)*x_{s} + A2(s)*x_{ss} + w(t)
% % With BCs    x(s=0) = 0
% %             x_{s}(s=1) = 0
% % And Output  z(t) = x(t,1)
%
% % gamma = 15.147 for lamb = 4.6 (this line should remain commented)
 lam = 4.6; 
 A0 = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;
 A1 = 3*s^2 - 2*s; 
 A2 = s^3 - s^2 + 2;

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
%%% evalin('base','Hinf_gain_dual = 1;')

if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;   PDE_b.nw = 1;   PDE_b.nz = 1;
 PDE_b.dom = [0,1];

 PDE_b.A0 = A0;  PDE_b.A1 = A1;   PDE_b.A2 = A2;
 PDE_b.B21 = s;
 PDE_b.C10 = [0,1,0,0];

 PDE_b.B = [1 0 0 0;
          0 0 0 1];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];
 
 % PDE: x_{t} = A0 * x + A1 * x_{s} + A2 * x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = [0; 1; 2];
 PDE_t.x{1}.term{1}.C = [A0, A1, A2];
 
 % PDE: x_{t} = ... + w;
 PDE_t.x{1}.term{2}.w = 1;
 
 % Output: z = x(s=1)
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.loc = 1;
 
 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;
 
 % BC 2: 0 = x_{s}(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.D = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Parabolic_Eq_with_Disturbance.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==24
% % % %====================================================================
% % Diffusion-Reaction equation from [12] (Example 3)
% % PDE         xi_{t} = lamb*xi + sum(xk_{ss},k=1,i) + w(t)
% % With BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)
%
% % gamma = 8.1069 for lamb = (1-1e-2)*pi^2 (this line should remain commented)
 lam = (1-1e-2)*pi^2; 
 ne = 1;

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

on = eye(ne);     ze = zeros(ne);     lt = tril(ones(ne));
 
if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz=1;
 PDE_b.dom = [0,1];
 PDE_b.A0 = lam*on;   PDE_b.A2 = lt;
 PDE_b.B21 = ones(ne,1);
 PDE_b.Ca1 = ones(1,ne);

 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = [];

 % PDE: x_{t} = lam * xi
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = lam;
 
 % PDE: x_{t} = ... + sum_{k=1}^{i}(xk_{ss})
 PDE_t.x{1}.term{2}.x = 1;
 PDE_t.x{1}.term{2}.D = 2;
 PDE_t.x{1}.term{2}.C = lt; % lower-triangular matrix
 
 % PDE: x_{t} = ... + w
 PDE_t.x{1}.term{3}.w = 1;
 
  % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BC 2: 0 = x(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Reaction_Diffusion_Eq_with_Disturbance.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==25
% % % %====================================================================
% % Diffusion-reation PDE from [6]
% % PDE         x_{t} = Cm*x + (1/R)*x_{ss} + [s;s]*w(t)
% % with BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x1(t,s),s,0,1)

 Cm = [1, 1.5; 5, 0.2];
 % gamma = 0.8102 for R = 2.6 (this line should remain commented)
 R = 2.6;


 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 ne = size(Cm,2);
 on = eye(ne);     ze = zeros(ne); 

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz = 1;
 PDE_b.dom = [0,1];

 PDE_b.A0 = Cm;
 PDE_b.A2 = (1/R)*eye(ne);
 PDE_b.B21 = s*ones(ne,1);
 PDE_b.Ca1 = [1 zeros(1,ne-1)];

 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = []; 
 
 % PDE: x_{t} = Cm * x + (1/R) * x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = [0; 2];
 PDE_t.x{1}.term{1}.C = [Cm, (1/R)*eye(ne)];
 
 % PDE: x_{t} = ... + [s;s] * w
 PDE_t.x{1}.term{2}.w = 1;
 PDE_t.x{1}.term{2}.C = s*ones(ne,1);
 
 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.C = [1,zeros(1,ne-1)];
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BC 2: 0 = x(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;
 
end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Reaction_Diffusion_Eq_with_Distributed_Disturbance.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==26
% % % %====================================================================
% % Diffusion-reaction PDE from [12]
% % PDE         x_{t} = Cm(s)*x + (1/R)*x_{ss} + [s;s;s]*w(t)
% % with BCs    x(s=0) = 0
% %             x(s=1) = 0
% % And Output  z(t) = int(x(t,s),s,0,1)

 Cm = [0,0,0; s,0,0; s^2,-s^3,0];
 % gamma =  4.23, for R = 21-1e-3 (this line should remain commented)
 R = (21-1e-3); 
 

 %%% Executive function
 evalin('base','Hinf_gain = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

ne = size(Cm,2);
on = eye(ne);     ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;   PDE_b.nw = 1;   PDE_b.nz = ne;
 PDE_b.dom = [0,1];

 PDE_b.A0 = Cm;
 PDE_b.A2= (1/R)*on;
 PDE_b.B21 = s*ones(ne,1);
 PDE_b.Ca1 = eye(PDE_b.n2);

 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.w{1}.vars = [];  PDE_t.z{1}.vars = []; 
 
 % PDE: x_{t} = Cm * x + (1/R) * x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = [0; 2];
 PDE_t.x{1}.term{1}.C = [Cm, (1/R)*eye(ne)];
 
 % PDE: x_{t} = ... + [s;s;s] * w
 PDE_t.x{1}.term{2}.w = 1;
 PDE_t.x{1}.term{2}.C = s*ones(ne,1);
 
 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BC 2: 0 = x(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Gain_Examples/Reaction_Diffusion_Eq_with_Distributed_Disturbance_2.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal observer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % % 1. Hyperbolic Transport, Balance Laws, Conservation Equations
% % % % %------------------------------------------------------------------
% % % Example 1D transport:
% % PDE             x_{t} = x_{s} + w(t)
% % with BC         x(s=1) = 0
% % And Outputs     z(t) = int(x(t,s),s,0,1) + w(t)
% %                 y(t) = int(x(t,s),s,0,1)
%
% % Answer: 1.0012  (this line should remain commented)

 %%% Executive function
 evalin('base','Hinf_estimator = 1;')
 ne = 1;
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 on = eye(ne);     ze = zeros(ne); 

if BATCH~=0
 %%% Batch input format
 PDE_b.nw = ne;   PDE_b.ny = ne;   PDE_b.nz = ne;   PDE_b.nx = 0;
 PDE_b.n0 = 0;    PDE_b.n1 = ne;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];
 
 PDE_b.A1 = on;
 PDE_b.Ca1 = on;   PDE_b.Ca2 = on;
 PDE_b.B21 = on;   PDE_b.D11 = on; 
   
 PDE_b.B = [ze on];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.z{1}.vars = [];  PDE_t.w{1}.vars = [];
 PDE_t.y{1}.vars = [];
 
 % PDE: x_{t} = x_{s}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = 1;
 
 % PDE: x_{t} = ... + w;
 PDE_t.x{1}.term{2}.w = 1;
 
 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;
 
 % Output: z = ... + w
 PDE_t.z{1}.term{2}.w  =1;
 
 % Output: y = int_{0}^{1} x(s) ds
 PDE_t.y{1}.term{1}.x = 1;
 PDE_t.y{1}.term{1}.I{1} = PDE_t.x{1}.dom;

 % BC 1: 0 = x(1)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 1;
 
end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Optimal_Observer_Examples/Transport_Eq_with_Observed_Output.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==28
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% % % Example Heat Equation
% % PDE             x_{t} = x_{ss} + w(t)
% % With BCs        x(s=0) = 0
% %                 x(s=1) = 0
% % And Outputs     z(t) = int(x(t,s),s,0,1) + w(t)
% %                 y(t) = x_{s}(s=1)
%
% % Answer: 1.0045  (this line should remain commented)

 %%% Executive function
 evalin('base','Hinf_estimator = 1;')
 ne = 1;
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 on = eye(ne);     ze = zeros(ne);

if BATCH~=0
 %%% Batch input format
 PDE_b.nw = ne;   PDE_b.ny = ne;   PDE_b.nz = ne;   PDE_b.nx = 0;
 PDE_b.n0 = 0;    PDE_b.n1 = 0;    PDE_b.n2 = ne;
 PDE_b.dom = [0,1];
 
 PDE_b.A2 = on; PDE_b.Ca1 = on;
 PDE_b.C20 = [ze ze ze on];  PDE_b.B21 = on;   D11 = on; 
 
 PDE_b.B = [on ze ze ze;
            ze on ze ze];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize finite-dimensional inputs and outputs.
 PDE_t.z{1}.vars = [];  PDE_t.w{1}.vars = [];
 PDE_t.y{1}.vars = [];
 
 % PDE: x_{t} = x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = 2;
 PDE_t.x{1}.term{1}.C = eye(ne);
 
 % PDE: x_{t} = ... + w
 PDE_t.x{1}.term{2}.w = 1;
 
 % Output: z = int_{0}^{1} x(s) ds
 PDE_t.z{1}.term{1}.x = 1;
 PDE_t.z{1}.term{1}.I{1} = PDE_t.x{1}.dom;
 
 % Output: z = ... + w
 PDE_t.z{1}.term{2}.w  =1;
 
 % Output: y = x_{s}(s=1)
 PDE_t.y{1}.term{1}.x = 1;
 PDE_t.y{1}.term{1}.D = 1;
 PDE_t.y{1}.term{1}.loc = 1;

 % BC 1: 0 = x(0)     
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;         

 % BC 2: 0 = x(1)
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Optimal_Observer_Examples/Heat_Eq_with_Observed_Output.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hinf optimal controller %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % % 2. Diffusive/Heat Equation Type Systems
% % % % %------------------------------------------------------------------
% % % Stabilizing controller for Heat Equation with z=0 and w=0:
% % PDE             x_{t} = lam*x + x_{ss} + u(t)
% % With BCs        x(s=0) = 0
% %                 x(s=1) = 0

 lam = 10;

 %%% Executive function
 evalin('base','Hinf_control = 1;')
 
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 %%% Batch input format
 PDE_b.nw = 0;   PDE_b.ny = 0;   PDE_b.nz = 0;   PDE_b.nx = 0;   PDE_b.nu = 1;
 PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 1;
 PDE_b.dom = [0,1];

 PDE_b.A0 = lam;   PDE_b.A2 = 1;
 PDE_b.B22 = 1;    %PDE_b.D12 = 1;
 
 PDE_b.B = [1 0 0 0;
            0 1 0 0];

end
if TERM~=0
 %%% Term-based input format
 % Initialize 1D PDE state component.
 PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];
 % Initialize a finite-dimensional input.
 PDE_t.u{1}.vars = [];
 
 % PDE: x_{t} = lam * x + x_{ss}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.D = [0; 2];
 PDE_t.x{1}.term{1}.C = [lam, 1];
 
 % PDE: x_{t} = ... + w
 PDE_t.x{1}.term{2}.u = 1;
 
 % BC 1: 0 = x(0)
 PDE_t.BC{1}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = 0;      

 % BC 2: 0 = x(1)        
 PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = 1;

end
if GUI
 %%% Associated GUI save file
 app = PIETOOLS_PDE_GUI; 
 load(fullfile(root,'Examples_PDE_GUI/Hinf_Optimal_Controller_Examples/Reaction_Diffusion_Eq_with_Controlled_Input.mat'));
 logval = app.loadData(data);
 if logval
     disp("Failed to load data object. Incorrect structure");
 end
end



elseif index==31
%%%####################################################################%%%%
% % # # # # # # # # # # # #     2D EXAMPLES     # # # # # # # # # # # # % %
%%%####################################################################%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % % % %==================================================================
% % % 2D transport equation
% % % PDE         x_{t}  = c1*x_{s1} + c2*x_{s2}
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;

%%% Executive Function:
 evalin('base','stability = 1;');

 c1 = 1; c2 = 1; 
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 
 % PDE: x_{t} = [c1, c2] * [x_{s1}; x_{s2}]
 PDE_t.x{1}.term{1}.D = [1,0; 0,1];
 PDE_t.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];
 
 % BC1: 0 = x(s1,0)
 PDE_t.BC{1}.term{1}.loc = [s1,0];
 % BC2: 0 = x(0,s2)
 PDE_t.BC{2}.term{1}.loc = [0,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end
 
 

elseif index==32
% % % % %==================================================================
% % % 2D reaction-diffusion equation (KISS Model)
% % % PDE         x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2s2}
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% % %             x(s1=1) = 0;   x(s2=1) = 0;
%
% % Stable when lam <= 2*pi^2 (Holmes, 1994 [14]).
    
%%% Executive Function:
 evalin('base','stability = 1;');

 c1 = 1; c2 = 2; lam = 19;
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 
 % PDE: x_{t} = [lam, c1, c2] * [x; x_{s1s1}; x_{s2s2}]
 PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
 PDE_t.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];
 
 % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
 PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
 PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==33
% % % % %==================================================================
% % % 2D advection-diffusion equation ((Holmes, 1994 [14])
% % % PDE         x_{t}   = C*(x_{s1s1} + x_{s2s2}) - b1*x_{s1} - b2*x_{s2}
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% % %             x(s1=1) = 0;   x(s2=1) = 0;

%%% Executive Function:
 evalin('base','stability = 1;');

 C = 1;  b1 = 0.5;  b2 = 2;
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 
 % PDE: x_{t} = [C, C] * [x_{s1s1}; x_{s2s2}];
 PDE_t.x{1}.term{1}.D = [2,0; 0,2];
 PDE_t.x{1}.term{1}.C = [C*eye(ne), C*eye(ne)];
 
 % PDE: x_{t} = ... + [b1, b2] * [x_{s1}; x_{s2}];
 PDE_t.x{1}.term{2}.D = [1,0; 0,1];
 PDE_t.x{1}.term{2}.C = [b1*eye(ne), b2*eye(ne)];
 
 % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
 PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
 PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==34
% % % % %==================================================================
% % % Telegraph equation (Holmes, 1994 [14])
% % % PDE         x_{t}   = -1/(2*lam) * x_{tt} + (C^2/2*lam)*(x_{s1s1} + x_{s2s2})
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% % %             x(s1=1) = 0;   x(s2=1) = 0;
% % %
% % % Use states x1 = x, x2 = x_{t}, then:
% % %
% % % PDE         x1_{t}   = x2
% % %             x2_{t}   = -2*lam * x2 + C^2*(x1_{s1s1} + x1_{s2s2})
% % % With BCs    x1(s1=0) = 0;   x1(s2=0) = 0;
% % %             x1(s1=1) = 0;   x1(s2=1) = 0;

%%% Executive Function:
 evalin('base','stability = 1;');

 C = 1;  lam = 1;
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
 
 % PDE: x1_{t} = x2;
 PDE_t.x{1}.term{1}.x = 2;
 
 % PDE: x2_{t} = -2*lam*x2
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.C = -2*lam;
 
 % PDE: x2_{t} = ... + [C^2, C^2] * [x1_{s1s1}; x1_{s2s2}];
 PDE_t.x{2}.term{2}.x = 1;
 PDE_t.x{2}.term{2}.D = [2,0; 0,2];
 PDE_t.x{2}.term{2}.C = [C^2 * eye(ne), C^2 * eye(ne)];
 
 % BC1: 0 = x1(s1,0)                    % BC3: 0 = x1(s1,1)
 PDE_t.BC{1}.term{1}.x = 1;             PDE_t.BC{3}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC2: 0 = x1(0,s2)                    % BC4: 0 = x1(1,s2)
 PDE_t.BC{2}.term{1}.x = 1;             PDE_t.BC{4}.term{1}.x = 1;
 PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==35
% % % 2D parabolic equation
% % % PDE         x_{t}  = a*x + b1*x_{s1} + b2*x_{s2} + c1*x_{s1s1} + c2*x_{s2s2}
% % % With BCs    x(s1=0) = 0;        x(s2=0) = 0;
% % %             x_{s1}(s1=1) = 0;   x_{s2}(s2=1) = 0;
%
% For c1=c2=1, b1=b2=0, will be stable when a <= 0.5*pi^2
% (Demetriou, 2019 [17]).
    
%%% Executive Function:
 evalin('base','stability = 1;');

 a = 4.9;  b1 = 0;  b2 = 0;  c1 = 1;  c2 = 1;
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 
 % PDE: x_{t} = a * x
 PDE_t.x{1}.term{1}.C = a*eye(ne);
 
 % PDE: x_{t} = ... + [b1, b2] * [x_{(1,0)}; x_{(0,1)}]
 PDE_t.x{1}.term{2}.D = [1,0; 0,1];
 PDE_t.x{1}.term{2}.C = [b1*eye(ne), b2*eye(ne)];
 
 % PDE: x_{t} = ... + [c1, c2] * [x_{(2,0)}; x_{(0,2)}]
 PDE_t.x{1}.term{3}.D = [2,0; 0,2];
 PDE_t.x{1}.term{3}.C = [c1*eye(ne), c2*eye(ne)]; 
 
 % BC1: 0 = x(s1,0)                     
 PDE_t.BC{1}.term{1}.loc = [s1,0];
 % BC2: 0 = x(0,s2)
 PDE_t.BC{2}.term{1}.loc = [0,s2];
 % BC3: 0 = x_{(0,1)}(s1,1)
 PDE_t.BC{3}.term{1}.D = [0,1];
 PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC4: 0 = x_{(1,0)}(1,s2)
 PDE_t.BC{4}.term{1}.D = [1,0];
 PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end
 
 

elseif index==36
% % % % %==================================================================
% % % 2D wave equation
% % % PDE         x_{tt} = c1*u_{(2,0)} + c2*x_{(0,2)};
% % % With BCs    x(s1=0) = 0;       x(s2=0) = 0;
% % %             x(s1=1) = 0;       x(s2=1) = 0;

%%% Executive Function:
 evalin('base','stability = 1;');

 c1 = 1;  c2 = 1;
 ne = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 
 % PDE: x_{tt} = [c1, c2] * [x_{s1s1}; x_{s2s2}]
 PDE_t.x{1}.tdiff = 2;
 PDE_t.x{1}.term{1}.D = [2,0; 0,2];
 PDE_t.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];
 
 % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
 PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
 PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==37
% % % % %==================================================================
% % % Isentropic compressible Navier-Stokes, linearized around pE=1, and
% % % vE=[vE1; vE2] = [s2;0] (Antonelli, 2021 [17])
% % % PDE         p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}
% % %             v1_{t} = -s2*v1_{s1} - v2 - (1/M^2)*p_{s1} + nu*(v1_{s1s1} + v1_{s2s2}) + lam*(v1_{s1s1} + v2_{s2s1})
% % %             v2_{t} = -s2*v2_{s1}      - (1/M^2)*p_{s2} + nu*(v2_{s1s1} + v2_{s2s2}) + lam*(v1_{s1s2} + v2_{s2s2})
% % % With BCs    p(s1=0) = 0;          p(s2=0) = 0;
% % %             v(s1=0) = 0;          v(s2=0) = 0;
% % %             v(s1=1) = 0;          v(s2=1) = 0;

%%% Executive Function:
 evalin('base','stability = 1;');

 M = 0.1;  lam = 1;  nu = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
 PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
 PDE_t.x{3}.vars = [s1;s2];   PDE_t.x{3}.dom = [0,1;0,1];
 
 % PDE: x1_{t} = -s2*x1_{s1}
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = -s2;
 
 % PDE: x1_{t} = ... - x2
 PDE_t.x{1}.term{2}.x = 2;
 PDE_t.x{1}.term{2}.C = -1;
 
 % PDE: x1_{t} = ... -(1/M^2)*p_{s1}
 PDE_t.x{1}.term{3}.x = 3;
 PDE_t.x{1}.term{3}.C = -(1/M^2);
 
 % PDE: x1_{t} = ... + nu*(v1_{s1s1} + v1_{s2s2})
 PDE_t.x{1}.term{4}.x = [1; 1];
 PDE_t.x{1}.term{4}.D = [2,0; 0,2];
 PDE_t.x{1}.term{4}.C = [nu, nu];
 
 % PDE: x1_{t} = ... + lam*(v1_{s1s2} + v2_{s2s2})
 PDE_t.x{1}.term{5}.x = [1; 2];
 PDE_t.x{1}.term{5}.D = [1,1; 0,2];
 PDE_t.x{1}.term{5}.C = [lam, lam];
 
 
 % PDE: x2_{t} = -s2*x2_{s1}
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.C = -s2;
 
 % PDE: x2_{t} = ... -(1/M^2)*p_{s1}
 PDE_t.x{2}.term{2}.x = 3;
 PDE_t.x{2}.term{2}.C = -(1/M^2);
 
 % PDE: x2_{t} = ... + nu*(v2_{s1s1} + v2_{s2s2})
 PDE_t.x{2}.term{3}.x = [2; 2];
 PDE_t.x{2}.term{3}.D = [2,0; 0,2];
 PDE_t.x{2}.term{3}.C = [nu, nu];
 
 % PDE: x2_{t} = ... + lam*(v1_{s1s1} + v2_{s2s1})
 PDE_t.x{2}.term{4}.x = [1; 2];
 PDE_t.x{2}.term{4}.D = [2,0; 1,1];
 PDE_t.x{2}.term{4}.C = [lam, lam];
 
 
 % PDE: x3_{t} = -s2*x3
 PDE_t.x{3}.term{1}.x = 3;
 PDE_t.x{3}.term{1}.C = -s2;
 
 % PDE: x3_{t} = ... - x1_{s1} - x2_{s2}
 PDE_t.x{3}.term{2}.x = [1;2];
 PDE_t.x{3}.term{2}.D = [1,0; 0,1];
 PDE_t.x{3}.term{2}.C = [-1, -1];
 
 
 % BC1: 0 = x1(0,s2);                   % BC2: 0 = x1(1,s2);
 PDE_t.BC{1}.term{1}.x = 1;             PDE_t.BC{2}.term{1}.x = 1;
 PDE_t.BC{1}.term{1}.loc = [0,s2];      PDE_t.BC{2}.term{1}.loc = [1,s2];
 % BC3: 0 = x1(s1,0);                   % BC4: 0 = x1(s1,1);
 PDE_t.BC{3}.term{1}.x = 1;             PDE_t.BC{4}.term{1}.x = 1;
 PDE_t.BC{3}.term{1}.loc = [s1,0];      PDE_t.BC{4}.term{1}.loc = [s1,1];
 
 % BC5: 0 = x2(0,s2);                   % BC6: 0 = x2(1,s2);
 PDE_t.BC{5}.term{1}.x = 2;             PDE_t.BC{6}.term{1}.x = 2;
 PDE_t.BC{5}.term{1}.loc = [0,s2];      PDE_t.BC{6}.term{1}.loc = [1,s2];
 % BC7: 0 = x2(s1,0);                   % BC8: 0 = x2(s1,1);
 PDE_t.BC{7}.term{1}.x = 2;             PDE_t.BC{8}.term{1}.x = 2;
 PDE_t.BC{7}.term{1}.loc = [s1,0];      PDE_t.BC{8}.term{1}.loc = [s1,1];
 
 % BC5: 0 = x3(0,s2);                   % BC6: 0 = x3(s1,0);
 PDE_t.BC{9}.term{1}.x = 3;             PDE_t.BC{10}.term{1}.x = 3;
 PDE_t.BC{9}.term{1}.loc = [0,s2];      PDE_t.BC{10}.term{1}.loc = [s1,0];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==38
% % % % %==================================================================
% % % 2D heat equation coupled to ODE
% % % ODE         x1_{t}  = (A+BK)*x1 + B*x2(s1=0,s2=0)
% % % PDE         x2_{t}  = c1*x_{(2,0)} + c2*x_{(0,2)}
% % % With BCs    x2_{(1,0)}(s1=0) = 0;     x2(s1=1) = 0;
% % %             x2_{(0,1)}(s2=0) = 0;     x2(s2=1) = 0;    

%%% Executive Function:
 evalin('base','stability = 1;');

 c1 = 1;  c2 = 1;   k = -2;
 nx = 1;    ne = 1;
 
 A = eye(nx);    B = eye(nx,ne);  K = k*eye(ne,nx);
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

nx = size(A,1);     ne = size(K,1);

if BATCH~=0
 disp('No batch input format available for this system, using terms-based format instead.')
 TERM = 1;
end
if TERM~=0
 %%% Term-based input format
 % Initialize an ODE state component.
 PDE_t.x{1}.vars = [];
 % Initialize a 2D PDE state component.
 PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
 
 % ODE: x1_{t} = (A+B*K)*x1
 PDE_t.x{1}.term{1}.x = 1;
 PDE_t.x{1}.term{1}.C = A+B*K;
 
 % ODE: x1_{t} = ... + B*x2(s1=0,s2=0)
 PDE_t.x{1}.term{2}.x =2;
 PDE_t.x{1}.term{2}.loc = [0,0];
 PDE_t.x{1}.term{2}.C = B;
 
 % PDE: x2_{t} = [c1, c2] * [x2_{s1s1}; x2_{s2s2}]
 PDE_t.x{2}.term{1}.x = 2;
 PDE_t.x{2}.term{1}.D = [2,0; 0,2];
 PDE_t.x{2}.term{1}.C = [c1*eye(ne), c2*eye(ne)];
 
 % BC1: 0 = x2_{s2}(s1,0)              
 PDE_t.BC{1}.term{1}.x = 2;
 PDE_t.BC{1}.term{1}.D = [0,1];
 PDE_t.BC{1}.term{1}.loc = [s1,0];      
 % BC2: 0 = x2_{s1}(0,s2)
 PDE_t.BC{2}.term{1}.x = 2;
 PDE_t.BC{2}.term{1}.D = [1,0];
 PDE_t.BC{2}.term{1}.loc = [0,s2];
 % BC3: 0 = x2(s1,1)
 PDE_t.BC{3}.term{1}.x = 2;
 PDE_t.BC{3}.term{1}.loc = [s1,1];
 % BC4: 0 = x2(1,s2)
 PDE_t.BC{4}.term{1}.x = 2;
 PDE_t.BC{4}.term{1}.loc = [1,s2];
 
end
if GUI~=0
 disp('No GUI representation available for this system.')
end



elseif index==39





elseif index==100
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



elseif index==101
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


end

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
    varargout{TERM} = PDE_t;
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