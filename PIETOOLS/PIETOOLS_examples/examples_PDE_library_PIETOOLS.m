 function varargout = examples_PDE_library_PIETOOLS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIETOOLS.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARARGOUT = EXAMPLES_PDE_LIBRARY_PIETOOLS(VARARGIN) returns the PDE of
% the example whose index is the first input, in the formats the other
% inputs request: a 'pde_struct' object in the term-based format (default;
% initialized and displayed) and/or a 'struct' in the batch format. It
% first resets the executive flags and forms in the base workspace; the
% example file then declares the system, sets the flags and the form of
% its test there (stability_form, Hinf_gain_form, H2_norm_form) and
% returns its recommended settings st. Finally the function asks whether
% to run the executive: on 'y' it converts the PDE to a PIE and, in the
% base workspace, assigns PIE and settings = st (clears settings when the
% example returns no st) and runs the executives the flags name. With no
% input, it asks for an example index.
%
% This library file contains the definition of some common ODE-PDE systems
% drawn from the literature. To use one of the examples, call this function
% with as FIRST input the index of the example you wish to use. Please
% check the table of contents below to see which index corresponds to which
% function. The output of the function then corresponds to a Matlab
% structure describing the desired PDE.
%
% Inputs:
%     (output1,output2) = examples_PDE_library_PIETOOLS(example,option1,option2,option3,...)
%       - example - integer (or possibly double) number (1-72; the variants 2.1-2.4, 4.1-4.5, 13.1, 13.2; 16 stops with an error) corresponding to selected example problem
%       - option1,option2,option3, etc.
%            - option# - 'batch' or 'BATCH' (or 1) - output# will be the PDE data structure in batch format  
%            - option# - 'terms' or 'TERMS' (or 2) - output# will be the PDE data structure in terms format
%            - option# - 'gui' or 'GUI' (or 3) - no output. Executes PIETOOLS PDE GUI
%             and loads saved problem data into GUI (examples 1-15, 17-28
%             and 51; the others print that no GUI representation exists).
%            - option# - 'lambda=XXX' overrides the default parameter value
%             lambda with value XXX for the given problem number. The list of allowable parameter 
%             overrides for each example problem is given in the table of
%             contents below. DO NOT ABUSE THIS POWER!
%             NOTE: Changing parameters does not change the GUI file!
%            - A format option given without an index selects a random
%             example.
% Outputs:
%     output1, output2 - the PDE in the formats requested, in their order
%             (at most two; none with 'gui' alone).
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
% Examples 1-15 and 17-29 have been defined using both the classical batch
% input format and the generalized term-based input format; examples 30-72
% only in the term-based format (a batch request returns the term-based
% PDE, with a message). To obtain a PDE struct in
% the classical format, call this function with as second argument the
% value "1" or the character array "batch". If you wish to
% define the PDE in the term-based input format, call this function with as
% second argument the value "2", or the character array "term". You can
% also specify both a value "1" and a value "2" (as second and third 
% argument) if you wish to obtain the PDE in both formats, in which case
% the outputs will be ordered in accordance with your arguments.
%
% In addition, to improve visualization and allow for tinkering, examples
% 1-15, 17-28 and 51 have an associated .mat file which can be loaded into
% the PDE GUI.
% To open the example in the GUI, call the function with an argument valued
% "3" or character array "GUI". Once again, you can load the GUI file in
% addition to producing a batch- or term-based PDE by specifying multiple
% arguments.
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
% PIETOOLS - examples_PDE_library_PIETOOLS
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% MP, 07/01/2020: Initial coding;
% DJ, 05/24/2021: Update to the Nth order format, and add 4th order
%                   implementation of Timoshenko beam equation;
% DJ, 05/31/2021: Adjust to function format;
% DJ, 12/29/2021: Update zero input case, and fix some examples;
% DJ, 08/25/2022: Update to new terms format and separate files;
% DB, MMP, DJ, 11/21/2026: Changes since PIETOOLS 2025: examples 53-55 (3-D)
%   and 56-72 (mixed domains) added, n_examples = 72; variants 13.1, 13.2 of
%   the Euler-Bernoulli beam; case 16 (the 4th order Timoshenko beam of the
%   entry of 05/24/2021) stops with an error: its BCs give an unstable
%   system; each case collects st, the recommended settings of its example
%   file; on 'y' settings = st is assigned (cleared without st) and the
%   executive runs through PIETOOLS_auto_execute_sop; the forms of the tests
%   are reset with the flags; the 2-D message names st; table of contents
%   (test, form, stability notion of each case), help and license revised.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process the user inputs
[index,BATCH,TERM,GUI,params] = process_inputs(varargin,nargin);



% Build the PDE of the example (and its st) in the case of its index
fprintf([' --- Extracting ODE-PDE example ', num2str(index),' ---\n']);
switch index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Table of Contents
%
% Each case calls its example file, which declares the system, sets the
% executive flags of its test (Test:) and the form of the test in the base
% workspace, and returns st, its recommended settings. Forms: 'P' a
% storage <T x, P T x> with P >= eppos I (for stability with epneg > 0:
% exponential stability); 'Q' the generalized storage T'Q; a dual test
% appends 'd'. H2 norm: the forms 'o', 'oco' and their duals 'c', 'cco'.
%__________________________________________________________________________
%   Problem                                         | Parameters            Notes / References
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
    case 1
% 	PDE: x_{t} = v*x_{s}                            | v=-1;   (Exponentially PDE stable ) 
%   BCs: x(s=0) = 0                                 |                (Finite Energy PDE stable ) 
%   Test: stability and its dual, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params);              % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params);            % DB, MMP, DJ, 11/21/2026
%-----|-------------------------------------------------------------------
    case {2.0,2.1,2.2,2.3,2.4}
%	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1] 
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling 
%        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]           | examples 2.1, 2.2, 
%                                                                            | 2.3, and 2.4.    
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable )   
%   Test: stability, form 'P'
%	[PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Lamare(index,GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
	[PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Lamare(index,GUI,params);               % DB, MMP, DJ, 11/21/2026
%-------------------------------------------------------------------------
    case 3
%   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
%        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
%   The implementation is based on a Linearized Saint–Venant–Exner Model
%   The dissipativity condition of Diagne 2012 [2] fails for default
%   parameters  k1 = -.5;       k2 = -4.5;                                      
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable )
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diagne(GUI,params);                    % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Diagne(GUI,params);                  % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case {4.0,4.1,4.2,4.4,4.3,4.5}
% 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}           | Different parameters          Saba 2019 [3] 
%        x2_{t} = sig2*x1 + (1/r2)*x2_{s}           | may be invoked calling  
%   BCs: x1(s=0) = qb*x2(s=0)                       | examples 4.1, 4.2, 
%        x2(s=1) = pb*x1(s=1)                       | 4.3, 4.4, and 4.5.  
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Saba(index,GUI,params);                % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Saba(index,GUI,params);              % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%        Diffusion and Heat-Equation Type Systems
%--------------------------------------------------------------------------
    case 5
%   PDE: x_{t} = lam*x + x_{ss}                     | lam = 9.86   (unstable for lam > pi^2,
%                                                                                          exponentially stable and FE PDE stable for lam < pi^2 = 9.8696)  
%   BCs: x(s=0) = 0,      x(s=1) = 0                |                               Ahmadi 2015 [5] 
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy FEPIE2PDE stable ) 
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq(GUI,params);     % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq(GUI,params);   % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 6
%   | PDE: x_{t} = lam*x + x_{ss}                   | lam = 2.466           (unstable for lam > pi^2/4 \approx 2.467,  
%                                                                                                     exponentially stable and FE PDE stable for lam < pi^2/4)
%   | BCs: x(s=0) = 0,      x_{s}(s=1) = 0          |                               Gahlawat 2017 [4]

%   Test: stability and its dual, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Gahlawat(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 7
%   PDE:  x_{t} = a(s)*x_{ss}                    | a = s^3 - s^2 + 2     (unstable for lam > 4.66)                                                                                                      
%                 + b(s)*x_{s}                          | b = 3*s^2 - 2*s                 exponentially stable and FE PDE stable for lam < 4.66 
%                 + c(s,lam)*x                          | c =-0.5*s^3 + 1.3*s^2            Gahlawat 2017 [4]
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0     |    - 1.5*s + 0.7 +lam
%                                                              | lam = 4.65                 
%   Test: stability and its dual, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params);     % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params);   % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 8
% 	ODE:  xo_{t} = A * xo - Bxr * x_{s}(s=a)        | a = 0    b = 1        (stable for lam < 9.7946; -Bxr: the sign of [7] reversed)
%   PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo      | lam = pi^2-1                  Das 2018 [7] (Example 2)
%   BCs:  x(s=a) = 0,     x(s=b) = 0                | A, Bxr, and Bpv fixed
%   Test: stability and its dual, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params);      % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das(GUI,params);    % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 9
% 	ODE:  xo_{t} = k*xo                             | k = -1                (stable for k<0)
%   PDE:  x_{t} = x_{ss}                            |
%   BCs:  x(s=0) = 0,     x(s=1) = xo               |
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 
%   Test: stability and its dual, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC(GUI,params);    % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_in_BC(GUI,params);  % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 10
%   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2                (stable for k=-2)            
%   PDE: x_{t} = x_{ss}                             |                               Tang 2011 [13]
%   BCs: x(s=0) = -xo,   x(s=1) = k*xo              |                           (PDE stable )
%   Test: stability, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params);     % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Tang(GUI,params);   % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 11
% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.9               (stable for R<2.7 (not tight))
%   BCs:  x(s=0) = 0,     x(s=1) = 0                | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6] Example D
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 
%   Test: stability and its dual, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 12
%   PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 30            (stable for R<21 (not tight))
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [0 0 0;                  Ahmadi 2015 [5] Adapted Example B
%                                                   |       s 0 0;
%                                                   |       s -s^2 0]
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Ahmadi_2(GUI,params); % DB, MMP, DJ, 11/21/2026
%
%--------------------------------------------------------------------------
%        Beam Type Equations 
%--------------------------------------------------------------------------
    case {13.0,13.1,13.2}
%   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
%   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |
%   13.1 use PDE states x1 = u_{t}, x2 = u_{ss}.             | (PDE stable)
%       =>                                          |
%   PDE: x1_{t} = -c * x2_{ss}                      |
%        x2_{t} = x1_{ss}                           |
%   BCs: x1(s=0) = 0,       x2(s=1) = 0             |
%        x1_{s}(s=0) = 0,   x2_{s}(s=1) = 0         |
% 13.2  use PDE states x1= u and x_2=u_{t} |                (PIE to PDE stable)

%   Test: stability and its dual, form 'P' (13, 13.1) or 'Q' (13.2)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq(index,GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                             Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi)-c
%        phi_{t}  
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
%   
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 14
%   Use states:                                     | k = 1                 Hyperbolic implementation | (PDE stable without damping c=0)
%        x1 = w_{t},   x2 = k*aa*g * (w_{s}-phi),   | aa = 1 
%        x3 = phi_{t},  x4 = E*II * phi_{s}.        | II = 1;  c = 0
%       =>                                          | g = 1
%   PDE: x1_{t} = (1/r/aa) * x2_{s}                 | E = 1   (PDE stable and Finite-energy PDE stable with damping c >0)
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3     | r = 1
%        x3_{t} = (1/r/II)*x2 + (1/r/II)*x4_{s}     |
%                 - c*x3                            |
%        x4_{t} = E*II * x3_{s}                     |
%   BCs: x1(0) = 0,         x2(1) = 0               |
%        x3(0) = 0,         x4(1) = 0               |
%                                                   |
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params);         % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Timoshenko_Beam_1(GUI,params);       % DB, MMP, DJ, 11/21/2026
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 15 
%   Assume all parameters are 1, and use states:    | Hyperbolic/diffusive implementation (c = 0: stable, not exponentially stable)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |                  (PDE stable with damping c >0)
%       =>                                          |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4 + x4_{ss} - c x3              |                   % DB, MMP, DJ, 11/21/2026
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0                          |
%   Test: dual stability, form 'Qd'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params);         % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Timoshenko_Beam_2(GUI,params);       % DB, MMP, DJ, 11/21/2026
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 16
%   We rewrite the system as a single equation:     |                       4th order implementation 
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
%    if BATCH~=0                                                            % DB, MMP, DJ, 11/21/2026 (was)
%        disp('Warning: No batch input formulation for this problem exists; reverting to term-based input format') % DB, MMP, DJ, 11/21/2026 (was)
%        BATCH = 0;         TERM = 1;                                       % DB, MMP, DJ, 11/21/2026 (was)
%    end                                                                    % DB, MMP, DJ, 11/21/2026 (was)
%    PDE_t = PIETOOLS_PDE_Ex_Timoshenko_Beam_3(GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
    % Not a stability example: with the BCs above, those of a published     % DB, MMP, DJ, 11/21/2026
    % simulation (Shivakumar et al., TAC 2025, Sec. VII-B), the system is   % DB, MMP, DJ, 11/21/2026
    % unstable (eigenvalues with real parts 1.52, 2.27, 2.75, ...); the     % DB, MMP, DJ, 11/21/2026
    % case stops with an error.                                             % DB, MMP, DJ, 11/21/2026
    error('Example 16 was removed from the stability examples on 11/19/2026.'); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Wave Equations
%--------------------------------------------------------------------------
    case 17
%   PDE: u_{tt} = u_{ss}                            | k=1    
%   BCs: u(s=0) = 0,                                |               Peet 2019 [8] (Example 8.2)
%        u_{s}(s=1) = -k*u_{t}(s=1)                 |
%   Use states x1 = u_{s}, x2 = u_{t}.              |     (PDE exponentially stable for k>0)                 
%       =>                                          |                      (Finite Energy PDE stable  for k>0) 
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params);   % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped(GUI,params); % DB, MMP, DJ, 11/21/2026
%-----|---------------------------------------------|----------------------
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1     Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 18
%   Use states x1 = u_{t}, x2 = u.                  | k = 1         (PDE exponentially dual stable )                 
%       =>                                          | ad = 1                    (Finite Energy PIE to PDE stable )   
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |
%   Test: dual stability, form 'Pd'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_1(GUI,params); % DB, MMP, DJ, 11/21/2026
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
    case 19
%   Use states x1 = u_{t}, x2 = u, x3 = u_{s}.      | k = 1  ( PDE exponentially stable)
%       =>                                          | ad = 1                    
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x3_{s}       |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x3(1) = 0                        |
%   Test: stability, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_2(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Datko_Boundary_Damped_2(GUI,params); % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-gain analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic/Transport/Balance Type Systems
%--------------------------------------------------------------------------
    case 20
%   PDE: x_{t} = -x_{s} + w(t)                      |                       (L2 gain gamma = 0.5)
%   BCs: x(s=0) = 0                                 |
%   Out: z(t) = int(x(t,s),s,0,1)                   |
%   Test: H-infinity gain and its dual, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 21
%   PDE: phi_{tt} = phi_{ss} + w(t)                 | k = 0.5               (L2 gain 0.8419 for k = 0.5)
%   BCs: phi(s=0) = 0                               |
%        phi_{s}(s=1) = -k*phi_{t}(s=1)             |
%   Out: z(t) = int(phi_{t}(t,s),s,0,1)             |
%   Use states x1 = phi_{s},    x2 = phi_{t}        |
%       =>                                          |
%   PDE: x1_{t} = x2_{s}                            |
%        x2_{t} = x1_{s} + w(t)                     |
%   BCs: x2(0) = 0,     x1(1) + k*x2(1) = 0         |
%   Test: H-infinity gain and its dual, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped(GUI,params);        % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped(GUI,params);      % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 22
%   PDE: x_{t} = x_{ss} + s*w(t)                    |
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |
%   Out: z(t) = int(x(t,s),s,0,1)                   |
%   Test: H-infinity gain, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 23
%   PDE: x_{t} = A2(s)*x_{ss}                       | A2 = s^3 - s^2 + 2;
%                + A1(s)*x_{s}                      | A1 = 3*s^2 - 2*s;
%                + A0(s)*x + w(t)                   | A0 =-0.5*s^3 +1.3*s^2
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |     -1.5*s +0.7 +lam          Shivakumar 2019 [12] (Example 1)
%   Out: z(t) = x(t,1)                              | lam = 4.6;               (gamma = 23.706 for lam = 4.6)
%   Test: H-infinity gain, form 'P'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 24
%   PDE: xi_{t} = w(t) + lamb*xi        i=1:ne      | lam = (1-1e-2)*pi^2   (gamma = 8.2141 for lamb = (1-1e-2)*pi^2, ne = 1)         
%                 + sum(xk_{ss},k=1,i)              | ne = 1 (state size)           Shivakumar 2019 [12] (Example 3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              |   
%   Out: z(t) = int(x(t,s),s,0,1)                   |
%   Test: H-infinity gain, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 25
%   PDE: x_{t} = Cm*x + (1/R)*x_{ss} + s*w(t)       | R = 2.6               (gamma = 0.8102 for R = 2.6)
%   BCs: x(s=0) = 0,        x(s=1) = 0              | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6]
%   Out: z(t) = int(x1(t,s),s,0,1)                  | 
%   Test: H-infinity gain, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 26
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + s*w(t)    | R = (21-1e-3)         (gamma =  2.1457, for R = 21-1e-3)
%   BCs: x(s=0) = 0,        x(s=1) = 0              | Cm = [0,0,0;                  Shivakumar 2019 [12] (Example 2)
%   Out: z(t) = int(x(t,s),s,0,1)                   |       s,0,0;
%                                                   |       s^2,-s^3,0]
%   Test: H-infinity gain, form 'Q'
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2(GUI,params); % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-optimal observer
%==========================================================================
%--------------------------------------------------------------------------
%       Hyperbolic Transport, Balance Laws, Conservation Equations
%--------------------------------------------------------------------------
    case 27
%   PDE: x_{t} = x_{s} + w(t)                       |                       (gain: 1.0004)
%   BCs: x(s=1) = 0                                 |
%   Out: z(t) = int(x(t,s),s,0,1) + w(t)            |
%        y(t) = int(x(t,s),s,0,1)                   |
%   Test: H-infinity optimal estimator
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 28
%   PDE: x_{t} = x_{ss} + w(t)                      | ne = 1 (state size)   (gain: 1.0005) 
%   BCs: x(s=0) = 0,        x(s=1) = 0              |
%   Out: z(t) = int(x(t,s),s,0,1) + w(t)            |
%        y(t) = x_{s}(s=1)                          |
%   Test: H-infinity optimal estimator
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output(GUI,params); % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-optimal controller
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 29
% % PDE: x_{t} = lam*x + x_{ss} + (s-s^2)u(t) + (s-s^2)w    | lam = 10;     (gain: 1.0000)  
% % BCs: x(s=0) = 0,        x(s=1) = 0
% % Out: z1(t) = int(x(s,t),s,0,1)
% %      z2(t) = u(t)
%   Test: H-infinity optimal controller
%    [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 30
% % ODE: x_{t}(t) = u(t);                                   | lam = 5;
% % PDE: v_{t}(t,s) = lam*v(t,s) +v_{ss}(t,s) +s*(2-s)*w(t) |                       adapted from Shivakumar 2022 [11] (Example 22)
% % BCs: v(t,s=0) = 0                                       |
% %      v_{s}(t,s=1) = x(t)                                |
% % Out: z1(t) = x(t)                                       |
% %      z2(t) = int_{0}^{1}v(t,s)ds                        |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity optimal controller
%    PDE_t = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%        Beam Type Equations 
%--------------------------------------------------------------------------
    case 31
% % PDE: v1_{t}(t,s) = -c*v2_{ss}(t,s) +s^2*w(t) +u(t,s);   | c = 0.1;
% %      v2_{t}(t,s) = v1_{ss}(t,s)                         |                       adapted from Shivakumar 2022 [11] (Example 21)  
% % BCs: v1(t,s=0) = v1_{s}(t,s=0) = 0;                     |
% %      v2(t,s=1) = v2_{s}(t,s=1) = 0;                     |
% % Out: z(t) = int_{0}^{1}(1-s)*v2(t,s)ds;                 |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity optimal controller
%    PDE_t = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Wave Equations
%--------------------------------------------------------------------------
% % ODE: xo_{t}(t) = u(t);                                  | c = 1;
% % PDE: x1_{t}(t,s) = x2(t,s);                             |                       adapted from Shivakumar 2022 [11] (Example 23) 
% %      x2_{t}(t,s) = c*x1_{ss}(t,s) + s*(2-s)*w;          |
% % BCs: x1(t,s=0) = 0;                                     |
% %      x2(t,s=0) = 0;                                     |
% %      x1_{s}(t,s=1) = xo(t);                             |
% % Out: z(t) = [int(x1(t,s),s,0,1); xo(t)]                 |
    case 32
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity optimal controller
%    PDE_t = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control(GUI,params);      % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control(GUI,params);  % DB, MMP, DJ, 11/21/2026
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
    case 33
%   PDE: x_{t} = c1*x_{s1} + c2*x_{s2}              | c1 = -1; c2 = -1;
%   BCs: x(s1=0) = 0,   x(s2=0) = 0;                |   ne = 1 (state size)
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Transport_Eq(GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Transport_Eq(GUI,params);               % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 34
%   PDE: x_{t} = lam*x + c1*x_{(2,0)}               | lam = 19;             (stable for lam <= 2*pi^2) 
%                           + c2*x_{(0,2)}          | c1 = 1;   c2 = 1;             Holmes, 1994 [14] Eq. (14)   
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;                |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params);                   % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params);                 % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 35
%   PDE: x_{t} = C*(x_{(2,0} + x_{(0,2)})           | C = 1;
%                 - b1*x_{(1,0)} - b2*x_{(0,1)}     | b1 = 0.5; b2 = 2;             Holmes, 1994 [14] Eq. (2)
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,                | ne = 1 (state size)
%        x(s1=1) = 0,   x(s2=1) = 0;                |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params);        % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params);      % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 36
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
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params);               % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 37
%   PDE: x_{t} = a*x                                | a = 4.0;              (stable for a <= 0.5*pi^2) 
%                 + b1*x_{(1,0)} + b2*x_{(0,1)}     | b1 = 0;   b2 = 0;             Demetriou, 2019 [16]
%                  + c1*x_{(2,0)} + c2*x_{(0,2)}    | c1 = 1;   c2 = 1;
%   BCs: x(s1=0) = 0,       x_{s1}(s1=1) = 0,       | ne = 1 (state size)
%        x(s2=0) = 0,       x_{s2}(s2=1) = 0;       |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Parabolic_Eq(GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Parabolic_Eq(GUI,params);               % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 38
%   PDE: x_{tt} = c1*x_{(2,0)} + c2*x_{(0,2)}       | c1 = 1;   c2 = 1;
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
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Wave_Eq(GUI,params);                      % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Wave_Eq(GUI,params);                    % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 39
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
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_NS_Antonelli(GUI,params);                 % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_NS_Antonelli(GUI,params);               % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 40
%   ODE: x1_{t} = (A+BK)*x1 + B*x2(s1=0,s2=0)       | A = I; B = I;
%   PDE: x2_{t} = c1*x2_{(2,0)} + c2*x2_{(0,2)}     | K = -2*I;
%   BCs: x2_{(1,0)}(s1=0) = 0;  x2(s1=1) = 0;       |
%        x2_{(0,1)}(s2=0) = 0;  x2(s2=1) = 0;       |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params);             % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params);           % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 41
%   PDE: x_{t} = x_{s1s1} + lam * x;    s1 in [0,1] | tau = 1;
%   BCs: x(t,s1=0) = 0;                             | lam = 1;                      Krstic, 2009 [18]
%        x(t,s1=1) = u(t-tau);                      |
%   Set x1(t)=x(t) and introduce delayed state:     |
%        x2_{t} = x2_{s2};          s2 in [1,1+tau] |
%        x2(t,1+tau) = u(t);                        |
%       =>                                          |
%   PDE: x1_{t} = x1_{s1s1} + lam * x1;             |
%        x2_{t} = x2_{s2};                          |
%   BCs: x1(t,s1=0)     = 0;                        |
%        x1(t,s1=1)     = x2(t,s2=1);               |
%        x2(t,s2=1+tau) = u(t);                     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 42
%   ODE: X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0) | tau = 1;
%   PDE: x_{t} = x_{s1s1} + a*x + a2*x(t-tau);      | A = -1;    A1 = 0.5;          Kang, 2017 [19] 
%   BCs: x_{s1}(t,s1=0) = 0;                        | B = -1;
%        x(t,s1=1) = 0;                             | a = 1;    a2 = 0.5;
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
%           x2(t,s1=1) = 0;                         |
%           x3(t,s2=1+tau) = x1(t);                 |
%           x4(t,s1,s2=1+tau) = x2(t,s1);           |    
%           x4_{s1}(t,s1=0,s2) = 0;                 |
%           x4(t,s1=1,s2) = 0;                      |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params);     % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params);   % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 43
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
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params);        % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params);      % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 44
%   PDE:    x_{tt}  = x_{ss}(t,s);   s in [0,1]     | tau = 1;  
%   BCs:    x(t,0) = 0;                             | k = 1;
%           x_{s}(t,1) = -k*(1-mu)*x_{t}(t,1)       | mu = 0.4;
%                           - k*mu*x_{t}(t-tau,1);  |
% Introduce:                                        | Stable for mu<0.5             Xu, 2006 [21]  
%     x1(t)    = - k*(1-mu)*x(1,t)                  |
%                   - k*mu*x(1,t-tau);              |
%     x2(t,s) = x(t,s);                             |
%     x3(t,s) = x_{t}(t,s);                         |
%     x4(t,s) = x(t-tau*s,1);                       |
%     x5(t,s) = x_{t}(t-tau*s,1);                   |
% Then:                                             |
%   ODE:    x1_{t}(t)    = x2_{s}(t,1);             |
%   PDE:    x2_{t}(t,1) = x3(t,s);                  |    s in (0,1)
%           x3_{t}(t,1) = x2_{ss}(t,s);             |       
%           x4_{t}(t,2) = -(1/tau) x4_{s}(t,2);     |
%           x5_{t}(t,2) = -(1/tau) x5_{s}(t,2);     |
%   BCs:    x2(t,0) = 0;                            |                         
%           x3(t,0) = 0;                            |
%           x4(t,0) = x2(t,1);                      |
%           x5(t,0) = x3(t,1);                      |
%           x1(t) = - k * (1-mu) * x2(t,1)          |
%                           - k * mu * x4(t,1);     |
%           x2_{s}(t,1) = -k*(1-mu)*x3(t,1)         |
%                               - k*mu*x5(t,1);     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end    
%   Test: stability, form 'Q'
%    [PDE_t] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(GUI,params);        % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(GUI,params);      % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-gain analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 45
%   PDE: x_{t} = r*x + nu*x_{(2,0)}                 | r = 5;  nu = 1;           (stable for r <= nu*pi^2*(1/(b-a)^2 +(1/4)/(d-c)^2)) 
%                           + nu*x_{(0,2)} + w      | a = 0;   b = 1;           (Exact L2 gain gam = 0.11101 for default params) 
%   OUT:     z = int(int(x,s1,a,b),s2,c,d)          | c = 0;   d = 1;
%   BCs: x(s1=a) = 0,   x(s2=c) = 0,                | ne = 1 (state size)
%        x(s1=b) = 0,   x_{s2}(s2=d) = 0;           |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity gain, form 'P'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN(GUI,params);       % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN(GUI,params);     % DB, MMP, DJ, 11/21/2026
%
%--------------------------------------------------------------------------
    case 46
%   PDE: x_{t} = r*x + nu*x_{(2,0)}                 | r = 15;  nu = 1;          (stable for r <= nu*pi^2*(1/(b-a)^2 +1/(d-c)^2)) 
%                           + nu*x_{(0,2)} + w      | a = 0;   b = 1;           (Exact L2 gain gam = 0.17110 for default params) 
%   OUT:     z = int(int(x,s1,a,b),s2,c,d)          | c = 0;   d = 1;
%   BCs: x(s1=a) = 0,   x(s2=c) = 0,                | ne = 1 (state size)
%        x(s1=b) = 0,   x(s2=d) = 0;                |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity gain, form 'P'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD(GUI,params);       % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD(GUI,params);     % DB, MMP, DJ, 11/21/2026
%
%--------------------------------------------------------------------------
    case 47
%   PDE: x_{t} = r*x + nu*x_{(2,0)}                 | r = 1;  nu = 1;           (stable for r <= nu*pi^2*((1/4)/(b-a)^2 +(1/4)/(d-c)^2)) 
%                           + nu*x_{(0,2)} + w      | a = 0;   b = 1;           (Exact L2 gain gam = 0.20668 for default params) 
%   OUT:     z = int(int(x,s1,a,b),s2,c,d)          | c = 0;   d = 1;
%   BCs: x(s1=a) = 0,        x(s2=c) = 0,           | ne = 1 (state size)
%        x_{s1}(s1=b) = 0,   x_{s2}(s2=d) = 0;      |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity gain, form 'P'
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN(GUI,params);       % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN(GUI,params);     % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-optimal observer
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 48
%   PDE: x_{t} = r*x + nu*x_{(2,0)}                 | r = 8;  nu = 1;           (stable for r <= nu*pi^2*((1/4)/(b-a)^2 +(1/4)/(d-c)^2)) 
%                      + nu*x_{(0,2)}               | a = 0;   b = 1;           
%                        + Cw(s1,s2)*w;             | c = 0;   d = 1;
%   OUT:    y1 = x(s2=d) + w1;                      | ne = 1 (state size)
%           y2 = x(s1=b) + w2;                      | Cw = (s1^2-b)*(s2^2-d); 
%           z = int(int(x,s1,a,b),s2,c,d)           |
%   BCs: x(s1=a) = 0,        x(s2=c) = 0,           |
%        x_{s1}(s1=b) = 0,   x_{s2}(s2=d) = 0;      |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity optimal estimator
%    [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve(GUI,params); % DB, MMP, DJ, 11/21/2026
% 
%==========================================================================
%       H2 gain
%==========================================================================
%--------------------------------------------------------------------------
%       Transport Equation Type Systems
%--------------------------------------------------------------------------
    case 49
%  PDE :        x_{t} = x_{s} + (s-s^2)*w(t)        |                           (H2 norm 0.1016, by numerical integration)
%  With BC:     x(s=1) = 0                          |
%  And output:  z(t) = int(x(t,s),s,0,1)            |
%                                                   |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H2 norm, form 'o', and its dual, form 'c'
%    [PDE_t] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2(GUI,params);  % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2(GUI,params); % DB, MMP, DJ, 11/21/2026
    %%
    case 50
% Original 1D PDE: x1_{t}(t,s1) = x1_{s1s1}(t,s1) - 20 x1(t,s1) - 4 x2(t,s1,-1) - 0.1x3(t,s1,-1) + s1(s1-\pi)w(t);
% 2D Transport Equation resulting from delay 1: x2_{t}(t,s1,s2) = 4 x2_{s2}(t,s1,s2);
%  2D Transport Equation resulting from delay 2: x3_{t}(t,s1,s2) = 3.3333 x3_{s2}(t,s1,s2);
%  (s1 in [0,pi], s2 in [-1,0]; x2(t,s1,-1) = x1(t-0.25,s1), x3(t,s1,-1) = x1(t-0.3,s1))
%  output z(t) = int_0^pi x1(t,s1) ds1;
%  BC1 0 = x1(t,0);
%  BC2 0 = x1(t,3.1416);
% BC3 0 = x2(t,s1,0) - x1(t,s1);
% BC4  0 = x3(t,s1,0) - x1(t,s1);
if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H2 norm, form 'oco', and its dual, form 'cco'
%    [PDE_t] =  PIETOOLS_PDDE_Ex(GUI,params);                               % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDDE_Ex(GUI,params);                              % DB, MMP, DJ, 11/21/2026
% 
%==========================================================================
%       H2 estimator
%==========================================================================
%--------------------------------------------------------------------------
%       Unstable Reaction-Diffusion Type Systems
%--------------------------------------------------------------------------
    case 51
%  PDE:         x_{t} = 3x+(s^2+0.2)x_{ss}          |
%                                   -(s^2)/2w(t)    |
%  With BC:    x(s=0) = 0, x_{s}(s=1)=0             |
%  And outputs:  z(t) = int(x(t,s),s,0,1)           |
%                y(t) = x(s=1) + w                  |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H2 optimal estimator
%    [PDE_t] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion(GUI,params);      % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion(GUI,params);    % DB, MMP, DJ, 11/21/2026
% 
%==========================================================================
%       H2 controller
%==========================================================================
%--------------------------------------------------------------------------
%       Unstable Reaction-Diffusion Type Systems
%--------------------------------------------------------------------------
    case 52
%  PDE :        x_{t} = 3x+(s^2+0.2)x_{ss} +s*u(t) +s^2w(t)
%  With BC     x(s=0) = 0, x_{s}(s=1)=0
%  And outputs  z(t) = [int(x(t,s),s,0,1)
%                                            u]
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H2 optimal controller
%    [PDE_t] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control(GUI,params); % DB, MMP, DJ, 11/21/2026 (was)
    [PDE_t,st] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control(GUI,params); % DB, MMP, DJ, 11/21/2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       3D Examples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 53
%   PDE: x_{t} = lam*x + C*(x_{(2,0,0)}             | lam = 1;  C = 1;          (stable iff lam < 3*pi^2*C/2 = 14.8044 for C = 1)
%                         + x_{(0,2,0)}             |
%                         + x_{(0,0,2)})            | ne = 1 (state size)
%   BCs: x(s1=0) = 0,        x(s1=1) = 0,           |
%        x(s2=0) = 0,        x_{s2}(s2=1) = 0,      |
%        x_{s3}(s3=0) = 0,   x(s3=1) = 0;           |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq(GUI,params);      % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 54
%   PDE: x_{t} = lam*x + C*(x_{(2,0,0)}             | lam = 1;  C = 1;          (stable iff lam < 3*pi^2*C + |b|^2/(4*C) = 30.4838
%                         + x_{(0,2,0)}             | b1 = 0.5; b2 = 1;                     for the default parameters)
%                         + x_{(0,0,2)})            | b3 = 1.5;
%                 + b1*x_{(1,0,0)} + b2*x_{(0,1,0)} | ne = 1 (state size)
%                 + b3*x_{(0,0,1)}                  | 
%   BCs: x(s1=0) = 0,   x(s2=0) = 0,   x(s3=0) = 0, |
%        x(s1=1) = 0,   x(s2=1) = 0,   x(s3=1) = 0; |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq(GUI,params);     % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 55
%   PDE: x1_{t} = C1*(x1_{(2,0,0)} + x1_{(0,2,0)}   | lam = 1;  C1 = 1;  
%                   + x1_{(0,0,2)}) + lam*x2        | C2 = 0.3;                 (stable iff |lam| < 3*pi^2*sqrt(C1*C2) = 16.2174
%        x2_{t} = C2*(x2_{(2,0,0)} + x2_{(0,2,0)}   | ne = 2 (state size)           for the default parameters)
%                   + x2_{(0,0,2)}) + lam*x1        | 
%   BCs: xi(s1=0) = 0,   xi(s2=0) = 0,   xi(s3=0) = 0,  i = 1,2
%        xi(s1=1) = 0,   xi(s2=1) = 0,   xi(s3=1) = 0;
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq(GUI,params);            % DB, MMP, DJ, 11/21/2026
%
% DB, MMP, DJ, 11/21/2026: BEGIN mixed-domain examples 56-72
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Mixed-domain examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% States on different domains (1D, 2D, 3D, edges, delay as transport),
% coupled in-domain or through boundary conditions. Thresholds marked
% [derived] are computed from the modes or by finite differences, not
% read in the cited paper; [read] and [abstract] mark results taken from
% the cited paper or its abstract.
%==========================================================================
%       Stability analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Delay as transport
%--------------------------------------------------------------------------
    case 56
%   PDE: x_{t} = x_{s1s1} + x_{s2s2} + r*x          | r = 2.9;                  (stable iff tau < 1.034726 for r = 2.9 [derived])
%                 - x(t-tau);  (s1,s2) in [0,pi]^2  | tau = 0.9*1.034726;           Fridman, 2009 [23] (1D version)
%   BCs: x = 0 on all edges.                        |
%   Delay as a 3D transport state:                  |
%        x1 = x;  x2_{t} = -x2_{s3},  s3 in [0,tau] |
%        x2(t,s1,s2,s3=0) = x1(t,s1,s2);            |
%        x1_{t} = ... - x2(t,s1,s2,s3=tau);         |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(GUI,params);    % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 57
%   PDE: z_{tt} = 0.1*z_{s1s1} - 2*z_{t}            | tau = 0.9*2.881401;       (stable iff tau < 2.881401 [derived])
%                 - z(t-tau,s1);  s1 in [0,pi]      |                               Fridman, 2009 [23], Ex. 2
%   BCs: z(t,0) = z(t,pi) = 0;                      |
%   v = z_{t}; y_{t} = -y_{s2}, s2 in [0,tau],      |
%        y(t,s1,s2=0) = z(t,s1);                    |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 58
%   Example 57 in Riemann form on [0,pi/2]:         | tau = 0.9*2.881401;       (stable iff tau < 2.881401 [derived])
%   PDE: e_{t} = v_{s1};                            |                               Fridman, 2009 [23], Ex. 2
%        v_{t} = 0.1*e_{s1} - 2*v                   |
%                 - int_0^s1 y(t,r,tau) dr;         |
%        y_{t} = -y_{s2},  s2 in [0,tau];           |
%   BCs: v(t,0) = 0;  e(t,pi/2) = 0;                |
%        y(t,s1,s2=0) = e(t,s1);                    |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 59
%   PDE: e_{t} = v_{s1};  v_{t} = e_{s1};  s1 in [0,1] | mu = 0.6;  k = 1;      (exponentially stable iff mu > 0.5 for k = 1,
%        y_{t} = -(1/tau)*y_{s2};          s2 in [0,1] | tau = 2;                   tau = 2 [derived]; mu > 1/2 [abstract])
%   BCs: v(t,0) = 0;                                |                               Xu, 2006 [21]
%        e(t,1) = -k*mu*v(t,1) - k*(1-mu)*y(t,1);   |
%        y(t,0) = v(t,1);                           |
%   (e = w_{s1}, v = w_{t} of w_{tt} = w_{s1s1})    |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 60
%   PDE: u_{t} = u_{s1s1};           s1 in [0,pi]   | mu0 = 1;  mu1 = 0.5;      (spectral abscissa -0.762 for the defaults
%        y_{t} = -(1/tau)*y_{s2};    s2 in [0,1]    | tau = 1;                      [finite differences])
%   BCs: u(t,0) = 0;                                |                               Nicaise, 2009 [24]
%        u_{s1}(t,pi) = -mu0*u(t,pi) - mu1*y(t,1);  |
%        y(t,0) = u(t,pi);                          |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
%       Edge and interface couplings
%--------------------------------------------------------------------------
    case 61
%   PDE: u_{t} = u_{s1s1} + u_{s2s2} + lam*u;       | sig = 1;                  (stable iff lam < 0.740174 for sig = 1 [derived])
%        v_{t} = sig*v_{s1s1} - u_{s2}(t,s1,1)      | lam = 0.9*0.740174;           Goldstein, 2006 [25] (box version)
%                 + lam*v;   (s1,s2) in [0,1]^2     |
%   BCs: u_{s1} = v_{s1} = 0 at s1 = 0, 1;          |
%        u(t,s1,0) = 0;  u(t,s1,1) = v(t,s1);       |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge(GUI,params);    % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 62
%   PDE: r_{t} = Dr*r_{s1s1} + nu*f(t,s1,0) - mu*r; | Dr = 2;  d = 1;           (stable iff lam < 13.9198 [derived, finite
%        f_{t} = d*(f_{s1s1} + f_{s2s2}) + lam*f;   | mu = 1;  nu = 1;              differences N = 60])
%   BCs: r = 0 at s1 = 0, 1;                        | lam = 0.9*13.919824;          Berestycki, 2013 [26] (box version)
%        f = 0 at s1 = 0, 1 and s2 = 1;             |
%        -d*f_{s2}(t,s1,0) = mu*r - nu*f(t,s1,0);   |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Road_Field_Berestycki(GUI,params);         % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 63
%   PDE: x1_{t} = (1+s1)*x1_{s1s1} + r1*x1;         | r1 = 0.9*14.507976;       (block thresholds r1* = 14.5080, r2* = 22.2060
%        x2_{t} = x2_{s1s1} + (1+s2^2)*x2_{s2s2}    | r2 = 0.9*22.205980;           [finite differences N = 60]; cascade
%                  + r2*x2;   (s1,s2) in [0,1]^2    |                               stable iff both blocks are [inferred])
%   BCs: x1 = 0 at s1 = 0, 1;  x2 = 0 at s1 = 0, 1; |
%        x2(t,s1,0) = x1(t,s1);  x2(t,s1,1) = 0;    |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 64
%   PDE: ui_{t} = ki*(ui_{s1s1} + ui_{s2s2})        | k1 = 1;  k2 = 2;  H = 1;  (stable iff lam < 13.6869 [derived, finite
%                  + lam*ui,  i = 1,2;  [0,1]^2     | lam = 0.9*13.686938;          differences N = 60])
%   BCs: -k1*u1_{s1}(t,0,s2) = H*(u2-u1)(t,0,s2);   |                               Hickson, 2009 [27] (1D problem)
%         k2*u2_{s1}(t,0,s2) = H*(u2-u1)(t,0,s2);   |
%        u1 = u2 = 0 at s1 = 1, s2 = 0, s2 = 1;     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles(GUI,params);     % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 65
%   PDE: x_{t} = diag(1,d)*(x_{s1s1} + x_{s2s2})    | d = 0.9*4.461640;         (stable iff d < 4.461640 [derived])
%                 + J*x,  J = [20,-25;30,-30];      | bc = 'N';                     Turing, 1952 [28]
%   BCs: bc = 'N': x_{s1} = 0 at s1 = 0, 1,         |
%                  x_{s2} = 0 at s2 = 0, 1;         |
%        bc = 'D': x = 0 on all edges, J -> 2*J;    |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair(GUI,params);         % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 66
%   PDE: q_{t} = -q_{s1} + d*q_{s2s2} - mu*q;       | d = 0.1;  mu = 1;         (stable iff b < 2.302685 [derived])
%                             (s1,s2) in [0,1]^2    | b = 0.9*2.302685;             Kang, 2022 [29]
%   BCs: q(t,0,s2) = b*int_0^1 q(t,s1,s2) ds1;      |
%        q(t,s1,0) = q(t,s1,1) = 0;                 |
%   [PDE,st]: st, the example-specific settings of the example file         % DB, MMP, DJ, 11/21/2026
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 67
%   PDE: e1_{t} = p_{s1};  p_{t} = e1_{s1};         | L = 0.5;                  (exponentially stable for the clamped
%        e2_{t} = q_{s1};  q_{t} = e2_{s1}-th_{s1}; |                               problem [abstract]; free end here)
%        th_{t} = th_{s1s1} - q_{s1};  s1 in [0,L]  |                               Marzocchi, 2002 [30]
%   BCs: e1(t,0) = 0;  q(t,L) = 0;                  |
%        th(t,0) = th(t,L) = 0;                     |
%        p(t,L) = q(t,0);  e1(t,L) = e2(t,0);       |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 68
%   PDE: e_{t} = v_{s1};  v_{t} = e_{s1};           |                           (not exponentially stable; energy
%        w_{t} = w_{s1s1};          s1 in [0,1]     |                               o(t^{-4}) [read])
%   BCs: v(t,1) = 0;  e(t,0) = -w_{s1}(t,0);        |                               Batty, 2016 [31]
%        w(t,1) = 0;  w(t,0) = v(t,0);              |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 69
%   PDE: u_{t} = v;  v_{t} = u_{s1s1} + u_{s2s2};   |                           (not exponentially stable; energy
%        w_{t} = w_{s1s1} + w_{s2s2};  [0,1]^2      |                               o(t^{-2/3}) [read])
%   BCs: u = v = w = 0 at s1 = 1, s2 = 0, s2 = 1;   |                               Batty, 2019 [32]
%        u_{s1}(t,0,s2) = -w_{s1}(t,0,s2);          |
%        w(t,0,s2) = v(t,0,s2);                     |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 70
%   PDE: ph_{t} = ps;  ps_{t} = ph_{s1s1}+ph_{s2s2};|                           (not uniformly stable [abstract])
%        de_{t} = et;                               |                               Beale, 1976 [33]
%        et_{t} = -et - de - ps(t,s1,0);            |
%   BCs: ph = ps = 0 at s1 = 0, 1 and s2 = 1;       |
%        ph_{s2}(t,s1,0) = -et(t,s1);               |
%        de = et = 0 at s1 = 0, 1;                  |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale(GUI,params); % DB, MMP, DJ, 11/21/2026
%--------------------------------------------------------------------------
    case 72
%   PDE: x1_{t} = -x1_{s1};             s1 in [0,1] | r2 = 0.9*pi^2;            (stable iff r2 < pi^2 [derived])
%        x2_{t} = x2_{s2s2} + r2*x2;  [0,1]^2       |
%   BCs: x1(t,0) = 0;                               |
%        x2(t,s1,0) = x1(t,s1);  x2(t,s1,1) = 0;    |
%   (x2 of order 0 in s1)                           |
%   [PDE,st]: st, the example-specific settings of the example file
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: stability, form 'Q'
    [PDE_t,st] = PIETOOLS_PDE_Ex_Transport_Edge_Diffusion(GUI,params);      % DB, MMP, DJ, 11/21/2026
%
%==========================================================================
%       Hinf-gain analysis
%==========================================================================
%--------------------------------------------------------------------------
%       Diffusive/Heat Equation Type Systems
%--------------------------------------------------------------------------
    case 71
%   PDE: x_{t} = x_{s1s1} + x_{s2s2} + w(t,s1);     |                           (Exact L2 gain gam = 0.04571818 [derived])
%   OUT: z(t,s1,s2) = x(t,s1,s2);                   |                               Lidstrom, 2016 [34]
%   BCs: x = 0 on all edges of [0,1]^2;             |
    if BATCH~=0
        disp('No batch input format available for this system, using terms-based format instead.')
        TERM = 1;
        BATCH = 0;
    end
%   Test: H-infinity gain, form 'P'
    [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom(GUI,params); % DB, MMP, DJ, 11/21/2026
% DB, MMP, DJ, 11/21/2026: END mixed-domain examples 56-72
%

%==========================================================================
%       Additional Examples (Undocumented)
%==========================================================================
% Cases 100 and 101 cannot be selected: process_inputs accepts indices up
% to n_examples (72).
case 100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional examples (undocumented) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %------------------------------------------------------------------
% % % % % ODE: xo_{t} = A*xo + B*x(t,0);  PDE: x_{t} = p*x_{ss} + (qc-q)*x;
% % % % % BCs: x_{s}(t,0) = 0,  x(t,1) = C*xo;  batch format only
% % % % %------------------------------------------------------------------

% Parameters and ODE matrices
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
% % % % % ODE: xo_{t} = A*xo + B*x_{s}(t,0);  PDE: x_{t} = p*x_{ss} + (qc-q)*x;
% % % % % BCs: x(t,0) = 0,  x_{s}(t,1) = C*xo;  batch format only
% % % % %------------------------------------------------------------------

% Parameters and ODE matrices
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

% Warn user that executives for 2D examples may take a lot of time.
if index>=33 && index<=48
%   fprintf('\n Examples 33 through 48 correspond to 2D PDEs, that take significant computational time and memory (use ''stripped'' settings).\n'); % DB, MMP, DJ, 11/21/2026 (was)
    fprintf('\n Examples 33 through 48 correspond to 2D PDEs, that take significant computational time and memory (the example file sets its settings, st).\n'); % DB, MMP, DJ, 11/21/2026
end

% Check if the user wants to run the executive
userinp = input('\n Would you like to run the executive associated to this problem? (y/n) \n ---> ','s');
userinp = strrep(userinp,'''','');    % Get rid of potential apostrophes
if strcmpi(userinp,'y') || strcmpi(userinp,'yes')
    PIE = convert_PIETOOLS_PDE(varargout{1});
    assignin('base','PIE',PIE);
    if exist('st','var'),  assignin('base','settings',st);                  % DB, MMP, DJ, 11/21/2026
    else,  evalin('base','clear settings');  end  % no st: run asks a tier  % DB, MMP, DJ, 11/21/2026
%   evalin('base','PIETOOLS_auto_execute');                                 % DB, MMP, DJ, 11/21/2026 (was)
    evalin('base','PIETOOLS_auto_execute_sop');                             % DB, MMP, DJ, 11/21/2026
elseif strcmpi(userinp,'n') || strcmpi(userinp,'no')
    return
else
    userinp = input(' Please answer y (yes) or n (n). \n ---> ','s');
    userinp = strrep(userinp,'''','');    % Get rid of potential apostrophes
    if strcmpi(userinp,'y') || strcmpi(userinp,'yes')
        PIE = convert_PIETOOLS_PDE(varargout{1});
        assignin('base','PIE',PIE);
        if exist('st','var'),  assignin('base','settings',st);              % DB, MMP, DJ, 11/21/2026
        else,  evalin('base','clear settings');  end                        % DB, MMP, DJ, 11/21/2026
%       evalin('base','PIETOOLS_auto_execute');                             % DB, MMP, DJ, 11/21/2026 (was)
        evalin('base','PIETOOLS_auto_execute_sop');                         % DB, MMP, DJ, 11/21/2026
    elseif strcmpi(userinp,'n') || strcmpi(userinp,'no')
        return
    else
       error("Response not recognized; execution of executive aborted.")
    end
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
% % [6] - arXiv version!
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
%   author = {Das, Amritam and Shivakumar, Sachin and Weiland, Siep and Peet, Matthew},
%   year = {2018},
%   pages = {},
%   title = {Representation and Stability Analysis of PDE-ODE Coupled Systems}
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
% % [11] -
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
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
%   title = {State and output feedback boundary control for a coupled PDE–ODE system},
%   journal = {Systems & Control Letters},
%   volume = {60},
%   number = {8},
%   pages = {540-545},
%   year = {2011},
%   issn = {0167-6911},
%   doi = {https://doi.org/10.1016/j.sysconle.2011.04.011},
%   url = {https://www.sciencedirect.com/science/article/pii/S0167691111000922},
%   author = {Shuxia Tang and Chengkang Xie},
%   keywords = {Coupled system, PDEs, Boundary control, Output feedback},
%   abstract = {This note is devoted to stabilizing a coupled PDE–ODE system with interaction at the interface. First, a state feedback boundary controller is designed, and the system is transformed into an exponentially stable PDE–ODE cascade with an invertible integral transformation, where PDE backstepping is employed. Moreover, the solution to the resulting closed-loop system is derived explicitly. Second, an observer is proposed, which is proved to exhibit good performance in estimating the original coupled system, and then an output feedback boundary controller is obtained. For both the state and output feedback boundary controllers, exponential stability analyses in the sense of the corresponding norms for the resulting closed-loop systems are provided. The boundary controller and observer for a scalar coupled PDE–ODE system as well as the solutions to the closed-loop systems are given explicitly.}
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
%   title={Feedback kernel approximations and sensor selection for controlled 2D parabolic PDEs using computational geometry methods},
%   author={Demetriou, Michael A and Hu, Weiwei},
%   booktitle={2019 IEEE 58th Conference on Decision and Control (CDC)},
%   pages={2144--2150},
%   year={2019},
%   organization={IEEE}
% }
%
% % [17] -
% @article{antonelli2021linear,
%   title={Linear stability analysis of the homogeneous Couette flow in a 2D isentropic compressible fluid},
%   author={Antonelli, Paolo and Dolce, Michele and Marcati, Pierangelo},
%   journal={Annals of PDE},
%   volume={7},
%   number={2},
%   pages={1--53},
%   year={2021},
%   publisher={Springer}
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
%   title = {Stability Analysis of the Heat Equation with Time-Delayed Feedback},
%   journal = {IFAC Proceedings Volumes},
%   volume = {42},
%   number = {6},
%   pages = {220-224},
%   year = {2009},
%   note = {6th IFAC Symposium on Robust Control Design},
%   issn = {1474-6670},
%   doi = {https://doi.org/10.3182/20090616-3-IL-2002.00038},
%   url = {https://www.sciencedirect.com/science/article/pii/S1474667015404057},
%   author = {Sina Yamaç çalişkan and Hitay özbay},
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
% % [22] -
% @inproceedings{peet2018new,
%   title={A new state-space representation for coupled PDEs and scalable Lyapunov stability analysis in the SOS framework},
%   author={Peet, Matthew M},
%   booktitle={2018 IEEE Conference on Decision and Control (CDC)},
%   pages={545--550},
%   year={2018},
%   organization={IEEE}
% }
% DB, MMP, DJ, 11/21/2026: BEGIN references [23]-[34]
% % [23] -
% @article{fridman2009exponential,
%   title={Exponential stability of linear distributed parameter systems with time-varying delays},
%   author={Fridman, Emilia and Orlov, Yury},
%   journal={Automatica},
%   volume={45},
%   number={1},
%   pages={194--201},
%   year={2009},
%   doi={10.1016/j.automatica.2008.06.006}
% }
%
% % [24] -
% @article{nicaise2009stability,
%   title={Stability of the heat and of the wave equations with boundary time-varying delays},
%   author={Nicaise, Serge and Valein, Julie and Fridman, Emilia},
%   journal={Discrete and Continuous Dynamical Systems - Series S},
%   volume={2},
%   number={3},
%   pages={559--581},
%   year={2009},
%   doi={10.3934/dcdss.2009.2.559}
% }
%
% % [25] -
% @article{goldstein2006derivation,
%   title={Derivation and physical interpretation of general boundary conditions},
%   author={Goldstein, Gis{\`e}le Ruiz},
%   journal={Advances in Differential Equations},
%   volume={11},
%   number={4},
%   pages={457--480},
%   year={2006}
% }
%
% % [26] -
% @article{berestycki2013influence,
%   title={The influence of a line with fast diffusion on Fisher-KPP propagation},
%   author={Berestycki, Henri and Roquejoffre, Jean-Michel and Rossi, Luca},
%   journal={Journal of Mathematical Biology},
%   volume={66},
%   pages={743--766},
%   year={2013},
%   note={arXiv:1210.3721}
% }
%
% % [27] -
% @article{hickson2009critical,
%   title={Critical times in multilayer diffusion. Part 1: Exact solutions},
%   author={Hickson, Roslyn I and Barry, Steven I and Mercer, Geoffry N},
%   journal={International Journal of Heat and Mass Transfer},
%   volume={52},
%   pages={5776--5783},
%   year={2009},
%   doi={10.1016/j.ijheatmasstransfer.2009.08.013}
% }
%
% % [28] -
% @article{turing1952chemical,
%   title={The chemical basis of morphogenesis},
%   author={Turing, Alan Mathison},
%   journal={Philosophical Transactions of the Royal Society of London. Series B},
%   volume={237},
%   pages={37--72},
%   year={1952}
% }
%
% % [29] -
% @article{kang2022effects,
%   title={Effects of diffusion on the principal eigenvalue for age-structured models with random diffusion},
%   author={Kang, Hao},
%   journal={Proceedings of the Royal Society of Edinburgh Section A},
%   volume={152},
%   number={1},
%   pages={258--280},
%   year={2022},
%   doi={10.1017/prm.2021.10}
% }
%
% % [30] -
% @article{marzocchi2002asymptotic,
%   title={Asymptotic behaviour and exponential stability for a transmission problem in thermoelasticity},
%   author={Marzocchi, Alfredo and Mu{\~n}oz Rivera, Jaime E and Naso, Maria Grazia},
%   journal={Mathematical Methods in the Applied Sciences},
%   volume={25},
%   number={11},
%   pages={955--980},
%   year={2002},
%   doi={10.1002/mma.323}
% }
%
% % [31] -
% @article{batty2016optimal,
%   title={Optimal energy decay in a one-dimensional coupled wave-heat system},
%   author={Batty, Charles J K and Paunonen, Lassi and Seifert, David},
%   journal={Journal of Evolution Equations},
%   volume={16},
%   pages={649--664},
%   year={2016},
%   note={arXiv:1505.04951}
% }
%
% % [32] -
% @article{batty2019optimal,
%   title={Optimal energy decay for the wave-heat system on a rectangular domain},
%   author={Batty, Charles J K and Paunonen, Lassi and Seifert, David},
%   journal={SIAM Journal on Mathematical Analysis},
%   volume={51},
%   number={2},
%   pages={808--819},
%   year={2019},
%   doi={10.1137/18M1195796}
% }
%
% % [33] -
% @article{beale1976spectral,
%   title={Spectral properties of an acoustic boundary condition},
%   author={Beale, J Thomas},
%   journal={Indiana University Mathematics Journal},
%   volume={25},
%   number={9},
%   pages={895--917},
%   year={1976}
% }
%
% % [34] -
% @inproceedings{lidstrom2016hinf,
%   title={H-infinity optimal control for infinite-dimensional systems with strictly negative generator},
%   author={Lidstr{\"o}m, Carolina and Rantzer, Anders and Morris, Kirsten},
%   booktitle={Proceedings of the 55th IEEE Conference on Decision and Control},
%   pages={5275--5280},
%   year={2016},
%   doi={10.1109/CDC.2016.7799077}
% }
% DB, MMP, DJ, 11/21/2026: END references [23]-[34]

function [index,BATCH,TERM,GUI,params] = process_inputs(varargin0,nargin1)
% Subroutine to process the user inputs for the
% examples_PDE_library_PIETOOLS function: the example index, the output
% formats and the parameter strings. It also sets silent_initialize_pde =
% true and resets the executive flags and forms in the base workspace.
% INPUT
% - varargin0: 'cell', the inputs of the library call;
% - nargin1:   'double', their number; 0 prompts for an example index.
% OUTPUT
% - index:     'double', the example index; random when the only input is
%              a format option;
% - BATCH, TERM: 'double', the output positions of the batch and the
%              term-based PDE, 0 for a format not requested (TERM = 1 when
%              no format is requested);
% - GUI:       'double', 1 to load the GUI save file of the example;
% - params:    'cell' of the parameter strings, each ending in ';'.


% Largest example index accepted
n_examples = 72;                                                            % DB, MMP, DJ, 11/21/2026

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
evalin('base','H2_norm = 0;')
evalin('base','H2_norm_dual = 0;')
evalin('base','H2_estimator = 0;')
evalin('base','H2_control = 0;')
% and the forms (empty: 'Q', for H2 'o', unless the example names one)      % DB, MMP, DJ, 11/21/2026
evalin('base','stability_form = '''';');                                    % DB, MMP, DJ, 11/21/2026
evalin('base','Hinf_gain_form = '''';');                                    % DB, MMP, DJ, 11/21/2026
evalin('base','H2_norm_form = '''';');                                      % DB, MMP, DJ, 11/21/2026

params = {};

% Collect the inputs
if nargin1==0 %<-- If there is no input, pause the script and let the user specify an example
    userinp = input('\n Select an example (1 through 50) to convert \n ---> ','s');
    varargin0 = split(userinp,[" ",","]);
    index = str2double(varargin0{1});
    if ~isnan(index) && (index>=0 && index<=n_examples)
        varargin0{1} = index;
    else
        % Invalid index: ask once more; a non-numeric answer is an error
        userinp = input('\n No existing example specified... Please input an integer value 1 through 48 to extract the example \n ---> ','s');
        varargin0 = split(userinp,[" ",","]);
        index = str2double(varargin0{1});
        if isnan(index)
            error('Please specify the desired example as first argument when calling ''examples_PDE_library_PIETOOLS''')
        else
            varargin0 = [str2double(userinp);varargin0];
        end
    end
end
nargin0 = length(varargin0);

if nargin0==1 %<-- If there is one input, this input MUST correspond to the desired PDE
    % An index alone: term-based output; a format option alone: a random
    % example in that format
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
% Two or more inputs: the index first
elseif nargin0>=2
    if isdouble(varargin0{1}) && varargin0{1}<=n_examples
        index = varargin0{1};
    else
        error('First argument must correspond to an example listed in this file')
    end
    % The other inputs: 'batch' (1), 'term' (2) set the output positions in
    % the order given, 'gui' (3) the GUI; other text that evaluates without
    % error is a parameter string; anything else is ignored with a warning
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
% No format requested: term-based output
if BATCH==0 && TERM==0 && GUI==0
    disp('Warning: The input format has not been (appropriately) defined, reverting to terms input format...');
    TERM=1;
end

end