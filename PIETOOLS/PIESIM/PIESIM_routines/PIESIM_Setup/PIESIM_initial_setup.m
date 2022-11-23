%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_initial_setup.m    PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routines performs the following operations:
% 1) Finds time-derivatives of the user-defined boundary inputs and disturbances and
% adds them to the field 'uinput'
% 2) Transforms initial conditions from primary states to fundamental states
% via differentiation
%
% Input:
% 1) uinput - user-defined boundary inputs and disturbances
% 2) psize - size of the problem. Includes nu, nw, nx, nf, N and n
% 3) type - of class ``char'' - type of the problem: 'PDE', 'DDE' or 'PIE'
%
% Output:
% 1) uinput - temporal derivatives of the user-defined boundary inputs and disturbances are added to the 
% field 'input' - these are needed for temporal integration of PIE equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021

function uinput=PIESIM_initial_setup(uinput,psize,type)
syms st sx;

% Compute temporal derivative of boundary inputs

 if isfield(uinput,'w')
             uinput.wdot=diff(uinput.w,st);
 end
 if isfield(uinput,'u')
             uinput.udot=diff(uinput.u,st);
 end

% Compute initial conditions on the fundamental states from initial
% condition on the primary states - only for PDE and DDE problems

if ~strcmp(type,'PIE')

ns=sum(psize);
psize_aux=[1 psize];
nsum=cumsum(psize);
nsump1=cumsum(psize_aux);

% Degree of smoothness 

 for i=1:length(psize)
 p(nsump1(i):nsum(i))=i-1;
 end

     for i=1:ns
      uinput.ic.PDE(i) =  diff(uinput.ic.PDE(i),sx,p(i));
     end

end
