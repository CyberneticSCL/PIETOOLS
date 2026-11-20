%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_NDSDDF_library_PIETOOLS.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLES_NDSDDF_LIBRARY_PIETOOLS is a script that defines, in the calling
% workspace, the system of the one uncommented example block: a neutral
% type system, the 'struct' NDS, with
%   x'(t) - sum_i Ei{i} x'(t-tau(i)) = A0 x(t) + sum_i Ai{i} x(t-tau(i)),
% or a differential difference system, the 'struct' DDF, with
%   x'(t) = A0 x(t) + B1 w(t) + B2 u(t) + Bv v(t),
%   z(t) = C1 x(t) + D11 w(t) + D12 u(t),
%   r_i(t) = Cr{i} x(t) + Drv{i} v(t),  v(t) = sum_i Cv{i} r_i(t-tau(i)),
% (fields not set are zero), together with the executive flags of the
% example, the form of its test (stability_form) where the block names
% one, and st, its recommended settings. The active block is the DDF of
% the multiple showering people (H-infinity control). A block is used by
% removing one comment level from its lines; text that carries two levels
% stays a comment.
%
% INPUT
% none.
% OUTPUT (variables of the calling workspace)
% - NDS or DDF: 'struct', the system of the active block;
% - st:         'struct', its recommended settings (lpisettings_sop);
% - the executive flags and the form variable the active block sets.
%
% This library file contains the definition of some common neutral-type
% (NDS) and DDF systems drawn from the literature. To use one of the examples, simply uncomment
% the variable definitions for that example. Please make sure no other
% example is uncommented, or you will have problems.
%
% The examples are grouped into various problem types. These types are as
% follows:
% 1. Stability Analysis Tests
% 2. Hinf-gain Analysis
% 3. Optimal Control Problems
% 4. Optimal Estimator Design Problems
%
% Most examples in each problem type include the option to solve the LPI
% according to their problem type. Hence, it is recommended that you also
% uncomment the toggle as well. e.g. stability=1, Hinf_gain=1, etc.
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - examples_NDSDDF_library_PIETOOLS
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
% Initial coding MP - 10_01_2020
% MMP, 11/19/2026: Example-specific settings st and the forms of the tests
%                   (stability_form, Hinf_gain_form) in every block with a
%                   certified setting, from the per-example search of
%                   11/19/2026 (maintainer: no default settings);
% MMP, 11/19/2026: the active block (showering people, DDF, Hinf_control)
%                   takes its regulated output, option 2, from eq. (15) of
%                   Peet, Automatica 2021 (arXiv 1910.03881); with every
%                   option commented out, z was absent and the bound fell
%                   to about 1e-3;
%                   the block carries its st;
%                   the delays are tau(i) = i, as in the source (were n);
% MMP, 11/19/2026: Release comment pass (comment-only).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Neutral Type Systems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % This is the version used in [1], as well as Han, 2005; Fridman, 2001; Fridman,Shaked
% % 2002, Han, 2004ab, [6] He et al. 2004; Lien and Chen, 2003[3]; Niculescu, 2000;
% stability=1
% stability_dual=1
% 
% NDS.A0=[-.9 .2;.1 -.9];
% NDS.Ai{1}=[-1.1 -.2;-.1 -1.1];
% NDS.Ei{1}=[-.2 0;.2 -.1];
% NDS.tau(1)= 2.2; 
% % st: the stripped tier, with upper-integral s-length 1; theta-lengths
% %   [2;1;2] (multiplier; upper; lower integral); slack lift degree at its
% %   minimum; stability form 'P' (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 1;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [2;1;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
%%% analytic max is supposedly 2.2255 (Han 2005) although no citation is given
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % This is the version used in [4] - Park and Won, 2000, Ex. 3
% stability=1
% stability_dual=1
% 
% NDS.A0=[-2 0;0 -1];
% NDS.Ai{1}=[0 .5;.5 0];
% NDS.Ei{1}=[.2 0;0 .2];
%  NDS.tau(1)= 20; 
% % st: the stripped tier, with lower-integral s-length 0; upper-integral
% %   s-length 0; theta-lengths [1;2;2] (multiplier; upper; lower integral);
% %   slack lift degree at its minimum; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [1;2;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% supposedly stable for tau<.7516 (verified by pietools. indeed, stable out to at least tau=20)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is the version used in [2] (park and Won, 1999a, Ex. 2.2), as well as [1]; 
% stability=1
% stability_dual=1
% 
% NDS.A0=[-2 .1;0 -2];
% NDS.Ai{1}=[.1 .1;-.1 .2];
% NDS.Ai{2}=[.1 0;0 .1];
% NDS.Ei{1}=[.1 .05;.02 .1];
% NDS.Ei{2}=[.05 0;0 .05];
% NDS.tau(1)= .3; 
% NDS.tau(2)= .6;
% % st: the stripped tier, with lower-integral s-length 0; upper-integral
% %   s-length 0; theta-lengths [1;2;2] (multiplier; upper; lower integral);
% %   slack lift degree at its minimum; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [1;2;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% Said to be stable for any tau1=.3, tau2=.6 (verified by pietools)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % This is the version used in [2](park and Won, 1999a, Ex. 2.3), as well as [1];
% stability=1
% stability_dual=1
% 
% NDS.A0=[-2 1;0 -1];
% NDS.Ai{1}=[0 .3;-.3 0];
% NDS.Ai{2}=[.1 -.05;.05 .1];
% NDS.Ei{1}=[0 -.1;-.1 0];
% NDS.Ei{2}=[.05 0;0 .05];
% NDS.tau(1)= 2.262; 
% NDS.tau(2)= .5;
% % st: the stripped tier, with lower-integral s-length 0; upper-integral
% %   s-length 0; theta-lengths [1;2;2] (multiplier; upper; lower integral);
% %   slack lift degree at its minimum; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [1;2;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% Said to be stable for any h1>0, h2<.552

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is an example from [7], Zhao, 2018
% stability=1
% stability_dual=1
% NDS.A0=[-2 0;0 -.9];
% NDS.Ai{2}=[-1 0;-1 -1];
% NDS.Ei{1}=[.1 0;0 .1];
% h=1.58; %max at 1.5804? (stability verified for h<2.04)
% NDS.tau=[h 3*h];
% % st: the stripped tier, with lower-integral s-length 1; slack lift degree
% %   at its minimum; stability form 'P' (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 1;                                                % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % This is an example from [7], Zhao, 2018 ?????
% % % % stability=1
% % % % stability_dual=1
% % % % NDS.A0=[-2 -1;3.5 0];
% % % % NDS.Ai{2}=[-1.5 0;0 0];
% % % % NDS.Ei{1}=[.9 0;0 0];
% % % % h=.77125;
% % % % %%% max at .7722? 
% % % % NDS.tau=[h 2*h];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % This is an example from [7], Zhao, 2018
% stability=1
% stability_dual=1
% NDS.A0=[-1];
% NDS.Ai{2}=[-2];
% NDS.Ei{1}=[.2];
% h=.43; %max at .5? (stability verified for h<.603)
% NDS.tau=[h 2*h];
% % st: the stripped tier, with upper-integral s-length 1; theta-lengths
% %   [2;1;2] (multiplier; upper; lower integral); slack lift degree at its
% %   minimum; stability form 'P' (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 1;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [2;1;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % An Example from Olgac and Sipahi, 2004
% stability=1
% stability_dual=1
% NDS.A0=[0 1 -1 0; 
% -3.346 -2.715 2.075 -2.007;
% -4 0 2 0;
% -3 0 0 6];
% NDS.Ai{1}=[-1 2 2 -1;
%     3 3 -2 0;
%     1 2 -1 1;
%     2 3 1 -3];
% NDS.Ei{1}=[.2 -.1 .5 -.1 ;
%     -.3 .09 -.15 -.027;
%     -3.333 .1 .2 1;
%     -1 2 .5 1];
%   NDS.tau(1)=.0275; 
% % st: the light tier, with slack lift degree at its minimum; stability
% %   form 'P' (storage <T x, P T x>)
% st = lpisettings_sop('light');                                            % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% stable for tau <.02755 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % An Example from Ramirez and Sipahi, 2020
% stability=1
% stability_dual=1
% NDS.A0=[-2 .2 -.3 0 -.4;
%     .2 -3.8 0 .7 0;
%     .8 0 -1.6 0 0;
%     0 .8 -.6 -2 .3;
%     -1 -.1 -1.5 0 -1.8];
% NDS.Ai{1}=[-2.2 0 0 1 0;
%     1.6 -2.2 1.6 0 0;
%     -0.2 -0.2 -0.2 -0.2 -0.2;
%     0 0.4 -1.4 -3.4 1;
%     -0.2 0.4 -0.1 -1.1 -3.3];
% NDS.Ei{1}=[ 0.40888 0.00888 0.20888 -0.09112 -0.29112;
%     0 0.2 0 0 0.6;
%     -0.1 -0.4 0 -0.8 0;
%     0 0 -0.1 0 0;
%     0 0 0 -0.2 -0.1];
%  NDS.tau(1)=.60; 
% % st: the stripped tier, with theta-lengths [2;2;1] (multiplier; upper;
% %   lower integral); slack lift degree at its minimum; stability form 'P'
% %   (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [2;2;1];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% stable for tau <.603
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Source Uncertain
% stability=1
% stability_dual=1
% NDS.A0=[-4 -1;0 -3];
% NDS.Ai{1}=[2 0;1 1];
% NDS.Ei{1}=+[.2 0;-.1 -.2];
% NDS.tau(1)=.5;
% % st: the stripped tier, with lower-integral s-length 0; upper-integral
% %   s-length 0; theta-lengths [1;2;2] (multiplier; upper; lower integral);
% %   slack lift degree at its minimum; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [1;2;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026

% % 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stability=1
% stability_dual=1
% 
% ki=10;kp=10;
% a=.4;b=50;h=.2;d=.8;sig=.3;
% al1=d+kp;gam1=b*ki*d+a*ki;
% al2=(d-kp)*exp(sig*h);gam2=(b*ki*d-a*ki)*exp(sig*h);
% beta1=(b*kp+a)*d+b*d^2+a*kp+ki;
% beta2=((b*kp+a)*d-b*d^2-a*kp-ki)*exp(sig*h);
% 
% NDS.A0=(1/al1)*[0 al1;-sig^2*al1+sig*beta1-gam1 -beta1+2*sig*al1];
% NDS.Ai{1}=(1/al1)*[0 0;-sig^2*al2+sig*beta2-gam2 -beta2+2*sig*al2];
% NDS.Ei{1}=[0 0;0 -al2/al1];
% NDS.tau(1)=h;
% % st: the stripped tier, with lower-integral s-length 0; upper-integral
% %   s-length 0; slack lift degree at its minimum; stability form 'P'
% %   (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% A Neutral-Type System: 
% %%% A 2-delay Neutral-Type System: 
% %%% \dot x(t)- \sum_i Ei{i}\dot x(t-\tau(i))=A0 x(t)+\sum_i Ai{i}x(t-\tau_i)
% stability=1
% stability_dual=1
% 
% NDS.Ei{1}=.2;NDS.Ei{2}=0;
% NDS.A0=-1;NDS.Ai{1}=0;NDS.Ai{2}=-2;
% h=.4; % max at .4448?
% NDS.tau=[h 2*h];
% % st: the stripped tier, with upper-integral s-length 1; theta-lengths
% %   [2;1;2] (multiplier; upper; lower integral); slack lift degree at its
% %   minimum; stability form 'P' (storage <T x, P T x>)
% st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
% st.lpi.pos.len_com{1} = 1;                                                % MMP, 11/19/2026
% st.lpi.pos.len_sos{1} = [2;1;2];                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% 

% ndim=size(An,1);
% stability=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Use these commands to convert the NDS to a DDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NDS=initialize_PIETOOLS_NDS(NDS);
% DDF=convert_PIETOOLS_NDS(NDS,'pie');
% DDF=minimize_PIETOOLS_DDF(DDF);
% PIE=convert_PIETOOLS_DDF(DDF,'pie');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Differential Difference Equations (DDFs) %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STABILITY TEST EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
% stability=1;
% stability_dual=1;
% DDF.Drv{1}=.5; DDF.Cv{1}=1;DDF.Drv{2}=.25; DDF.Cv{2}=1;
% DDF.tau(1)=1; DDF.tau(2)=2;
% % st: the extreme tier, with slack lift degree at its minimum; slack
% %   integral bases capped at degree Dmax-1; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('extreme');                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% st.lpi.ineq.cap_inc = 'auto';                                             % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%% A pure Difference Equation r_1(t)=Drv{1}*Cv{1}*r_1(t-tau(1))+Drv{2}*Cv{2}*r_2(t-tau(2))
% stability=1;
% stability_dual=1;
% DDF.Drv{1}=.5; DDF.Cv{1}=1;DDF.Drv{2}=.25; DDF.Cv{2}=1;
% DDF.tau(1)=1; DDF.tau(2)=2;
% % st: the extreme tier, with slack lift degree at its minimum; slack
% %   integral bases capped at degree Dmax-1; stability form 'P' (storage
% %   <T x, P T x>)
% st = lpisettings_sop('extreme');                                          % MMP, 11/19/2026
% st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
% st.lpi.ineq.cap_inc = 'auto';                                             % MMP, 11/19/2026
% stability_form = 'P';                                                     % MMP, 11/19/2026
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONTROLLER SYNTHESIS EXAMPLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%% Multiple Showering People - tracking with integral control
% % % This is the DDF implementation
Hinf_control=1
n=2; %This is the number of users and can be changed
shower_ex=2;
nndelay=n;
% the delay for user i is tau(i)=i;
for i=1:nndelay
%   DDF.tau(i)=n;                                                           % MMP, 11/19/2026 (was)
    DDF.tau(i)=i;   % as the comment above and the source, eq. (15)         % MMP, 11/19/2026
end
% set all alphas to 1
alpha=1;
alphav=alpha*ones(n,1);
% set all gammas to 1/n
gamma=1/n;
gammaM=gamma*ones(n,n);
zn=zeros(n);

% State x = [x1; x2], x1, x2 in R^n: x1' = x2 - w, x2' = u - Gam*w + Gam*v
DDF.A0=[zn eye(n);zn zn];

% Coupling: Gam(i,i) = -alpha_i, Gam(i,j) = alpha_j*gamma_ij for i ~= j
Gam=zn;
for i=1:n
    for j=1:n
        if i~=j
            Gam(i,j)=alphav(j)*gammaM(i,j);
        else
            Gam(i,j)=-alphav(i);
        end

    end
end
DDF.B1=[-eye(n);-Gam];
DDF.B2=[zn;eye(n)];

% %There are 3 option for the number of outputs in this example:
% %%%% 1) Full outputs, Full disturbances
% DDF.B1=[-eye(n);-Gam];
% DDF.C1=[eye(n) zn;zn zn];
% DDF.D12=[zn;.1*eye(n)];

% %%%%% 2) 2 outputs, Full disturbances
% DDF.B1=[-eye(n);-Gam];
% Option 2, active: the regulated output of the source, z = [sum_i T1i;     % MMP, 11/19/2026
% .1*sum_i ui], eq. (15) of Peet, Automatica 2021 (arXiv 1910.03881).       % MMP, 11/19/2026
% With every option commented out, Hinf_control has no regulated output.    % MMP, 11/19/2026
% The source sets tau(i) = i, as the loop above does.                       % MMP, 11/19/2026
% DDF.C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];                               % MMP, 11/19/2026 (was)
DDF.C1=[ones(1,n) zeros(1,n);zeros(1,2*n)];                                 % MMP, 11/19/2026
% DDF.D11=[zeros(2,n)];                                                     % MMP, 11/19/2026 (was)
DDF.D11=[zeros(2,n)];                                                       % MMP, 11/19/2026
% DDF.D12=[zeros(1,n);.1*ones(1,n)];                                        % MMP, 11/19/2026 (was)
DDF.D12=[zeros(1,n);.1*ones(1,n)];                                          % MMP, 11/19/2026
%
% % 3) Now for new DDF terms
DDF.Bv=[zn;Gam];

% Delayed channels: r_i = x2_i (Cr{i}), v_i(t) = r_i(t-tau(i)) (Cv{i})
for i=1:nndelay
    DDF.Cr{i}=[zeros(1,n+i-1) 1 zeros(1,n-i)];
DDF.Cv{i}=zeros(n,1);DDF.Cv{i}(i,1)=1;
end
% st: the heavy tier, with theta-lengths [4;2;2] (multiplier; upper; lower
%   integral); slack integral bases capped at degree Dmax-1; free operators
%   of multiplier degree 1, integral degree 2
st = lpisettings_sop('heavy');                                              % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [4;2;2];                                            % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026
st.lpi.var.deg{1} = [1 2 2];                                                % MMP, 11/19/2026
%%%gam_guess; % min in [.32 .38]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] - 
% @article{chen2001criteria,
%   title={Criteria for asymptotic stability of a class of neutral systems via a LMI approach},
%   author={Chen, JD and Lien, CH and Fan, KK and Chou, JH},
%   journal={IEE Proceedings-Control Theory and Applications},
%   volume={148},
%   number={6},
%   pages={442--447},
%   year={2001},
%   publisher={IET}
% }
%
% [2] - 
% @article{park1999asymptotic,
%   title={Asymptotic stability of neutral systems with multiple delays},
%   author={Park, Ju H and Won, S},
%   journal={Journal of Optimization Theory and Applications},
%   volume={103},
%   number={1},
%   pages={183--200},
%   year={1999},
%   publisher={Springer}
% }
% 
% [3] - 
% @article{lien2000stability,
%   title={Stability conditions for a class of neutral systems with multiple time delays},
%   author={Lien, Chang-Hua and Yu, Ker-Wei and Hsieh, Jer-Guang},
%   journal={Journal of Mathematical Analysis and Applications},
%   volume={245},
%   number={1},
%   pages={20--27},
%   year={2000},
%   publisher={Elsevier}
% }
%
% [4] -
% @article{park2000stability,
%   title={Stability analysis for neutral delay-differential systems},
%   author={Park, Ju H and Won, Sangchul},
%   journal={Journal of the Franklin institute},
%   volume={337},
%   number={1},
%   pages={1--9},
%   year={2000},
%   publisher={Elsevier}
% }
%
% [5] - 
% @article{han2005stability,
%   title={On stability of linear neutral systems with mixed time delays: a discretized Lyapunov functional approach},
%   author={Han, Qing-Long},
%   journal={Automatica},
%   volume={41},
%   number={7},
%   pages={1209--1218},
%   year={2005},
%   publisher={Elsevier}
% }
% 
% [6] - 
% @article{he2004delay,
%   title={Delay-dependent robust stability criteria for uncertain neutral systems with mixed delays},
%   author={He, Yong and Wu, Min and She, Jin-Hua and Liu, Guo-Ping},
%   journal={Systems \& Control Letters},
%   volume={51},
%   number={1},
%   pages={57--65},
%   year={2004},
%   publisher={Elsevier}
% }
% 
% [7] - 
% 
% @article{zhao2018necessary,
%   title={Necessary conditions for exponential stability of linear neutral type systems with multiple time delays},
%   author={Zhao, Ning and Zhang, Xian and Xue, Yu and Shi, Peng},
%   journal={Journal of the Franklin Institute},
%   volume={355},
%   number={1},
%   pages={458--473},
%   year={2018},
%   publisher={Elsevier}
% }
% 
% [8]-
% @article{olgac2004practical,
%   title={A practical method for analyzing the stability of neutral type LTI-time delayed systems},
%   author={Olgac, Nejat and Sipahi, Rifat},
%   journal={Automatica},
%   volume={40},
%   number={5},
%   pages={847--853},
%   year={2004},
%   publisher={Elsevier}
% }

