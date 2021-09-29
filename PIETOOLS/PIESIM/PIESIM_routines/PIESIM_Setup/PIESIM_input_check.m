%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_input_check.m     PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PDE, uinput]=PIESIM_input_check(PDE, uinput)
% Check if necessary inputs are defined to start the simulations
% NOTE: All other variables will be checked in PIETOOLS converter
% Inputs: PDE - PDE structure
%         uinput - user's input structure
% Outputs: PDE - updated PDE structure with variables properly defined (if
% previously undefined). 
%          uinput - updated user's input structure with variables properly defined (if
% previously undefined). 
% All properly defined variables are uchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_1_2021


 if ~isfield(PDE,'n0')
 disp('Warning: PDE.n0 is not defined. Defaulting to zero');
 PDE.n0=0;
 end
 if ~isfield(PDE,'n1')
 disp('Warning: PDE.n1 is not defined. Defaulting to zero');
 PDE.n1=0;
 end
 if ~isfield(PDE,'n2')
 disp('Warning: PDE.n2 is not defined. Defaulting to zero');
 PDE.n2=0;
 end
 if ~isfield(PDE,'nf')
 disp('Warning: PDE.nf is not defined. Defaulting to zero');
 PDE.nf=0;
 end
 if ~isfield(PDE,'nw')
 disp('Warning: PDE.nw is not defined. Defaulting to zero');
 PDE.nw=0;
 end
 if ~isfield(PDE,'nu')
 disp('Warning: PDE.nu is not defined. Defaulting to zero');
 PDE.nu=0;
 end
 if ~isfield(PDE,'nx')
 disp('Warning: PDE.nx is not defined. Defaulting to zero');
 PDE.nx=0;
 end
 ns=PDE.n0+PDE.n1+PDE.n2;
 
 if ~isfield(uinput,'ic')
 disp('Warning: PDE initial conditions are not defined. Defaulting to zero');
 disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
 uinput.ic.PDE(1:ns)=0;
 uinput.ic.ODE(1:PDE.nx)=0;
 end
 
 if ~isfield(uinput.ic,'PDE')
 disp('Warning: PDE initial conditions are not defined. Defaulting to zero');
 uinput.ic.PDE(1:ns)=0;
 end
 
 if (PDE.nx>0)
 if ~isfield(uinput.ic,'ODE')
 disp('Warning: ODE initial conditions are not defined. Defaulting to zero');
 uinput.ic.ODE(1:PDE.nx)=0;
 end
 end
 
 if(size(uinput.ic.PDE,2)<ns)
 disp('Warning: Number of initial conditions on PDE states is less than the number of PDE states');
 disp('Defalting the rest to zero');
 uinput.ic.PDE(1,size(uinput.ic.PDE,2)+1:ns)=0;
 end
 
 if(size(uinput.ic.PDE,2)>ns)
 disp('Warning: Number of initial conditions on PDE states is greater than the number of PDE states');
 disp('Defaulting all initial conditions to zero');
 uinput.ic=rmfield(uinput.ic,'PDE');
 uinput.ic.PDE(1:ns)=0.;
 end
 
 if (PDE.nx>0) 
 if(size(uinput.ic.ODE,2)<PDE.nx)
 disp('Warning: Number of initial conditions on ODE states is less than the number of ODE states');
 disp('Defalting the rest to zero');
 uinput.ic.ODE(size(uinput.ic.PDE,2)+1:PDE.nx)=0;
 end
 
 if(size(uinput.ic.ODE,2)>PDE.nx)
 disp('Warning: Number of initial conditions on ODE states is greater than the number of ODE states');
 disp('Defaulting all ODE initial conditions to zero');
 clear uinput.ic.ODE;
 uinput.ic=rmfield(uinput.ic,'ODE');
 uinput.ic.ODE(1:PDE.nx)=0.;
 end
 
 
 end

 
 if (PDE.nu>0)
 if ~isfield(uinput,'u')
 disp('Warning: nu is greater than zero, but user-defiened u inputs are not provided. Defaulting PDE.nu to zero.');
 PDE.nu=0;
 end
 end
 
 if (PDE.nu>0)
 if (size(uinput.u,2)<PDE.nu)
 disp('Warning: Number of provided u inputs is less than nu.');    
 disp('Defalting the rest of u inputs and their time derivatives to zero');
 uinput.u(size(uinput.u,2)+1:PDE.nu)=0;
 uinput.udot(size(uinput.u,2)+1:PDE.nu)=0;
 end
 
 if (size(uinput.u,2)>PDE.nu)
 disp('Warning: Number of provided u inputs is  greater than nu.');    
 disp('Defalting PDE.nu to zero');
 PDE.nu=0;
 end
 end
 
 if (PDE.nw>0)
 if ~isfield(uinput,'w')
 disp('Warning: nw is greater than zero, but user-defiened w inputs are not provided. Defaulting PDE.nw to zero.');
 PDE.nw=0;
 end
 end
 
 if (PDE.nw>0)
 if (size(uinput.w,2)<PDE.nw)
 disp('Warning: Number of provided w inputs is less than nw.');    
 disp('Defalting the rest of w inputs and their time derivatives to zero');
 uinput.w(size(uinput.w,2)+1:PDE.nw)=0;
 uinput.wdot(size(uinput.w,2)+1:PDE.nw)=0;
 end
 
 if (size(uinput.w,2)>PDE.nw)
 disp('Warning: Number of provided w inputs is greater than nw.');    
 disp('Defalting PDE.nw to zero');
 PDE.nw=0;
 end
 end
 
 
 if ~isfield(uinput,'ifexact')   
 disp('Watning: uinput.ifexact is not specified. Defaulting  to false');
 uinput.ifexact=false;
 end
 
 if(uinput.ifexact)
 if ~isfield(uinput,'exact')
 disp('Warning: exact solution is not provided. Defaulting uinput.ifexact to false');
     uinput.ifexact=false;
 end 
 end
 
 if(uinput.ifexact)
 if (size(uinput.exact,2)~=ns)
 disp('Warning: number of exact solutions provided does not match the number of PDE states');
 disp('Defaulting uinput.ifexact to false');
     uinput.ifexact=false;
 end 
 end
 
 if ~isfield(PDE,'dom')
 disp('PDE domain is not defined. Defaulting to [-1, 1]');
 PDE.dom=[-1 1];
 end
 
 if (PDE.dom(1)==PDE.dom(2))
 disp('Warning: left and right ends of the domain are the same. Defaulting domain to [-1, 1]');
 PDE.dom=[-1 1];
 end
 
 

 
 
 
  



