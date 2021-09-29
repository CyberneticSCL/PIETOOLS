%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_Input_domain_transform.m     PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE=PIESIM_domain_transform(PDE)
% This routine transforms the physical PDE domain x\in[a,b] into the
% computational domain xi\in[-1,1]
% Inputs: PDE - PDE structure
% Outputs: PDE - updated PDE structure that uses the change of variables to
% transform the domain x\in[a,b] into xi\in-[1,1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021

pvar s;

a=PDE.dom(1);
b=PDE.dom(2);
nbc=2*(PDE.n1+2*PDE.n2);

% Transform PDE to the computational domain x\in[-1,1]
% Warning: all PIE operators need to be evaluated on the computational
% domain with PDE.dom=-1, PDE.dom=1. a and b signify an original domain and 
% are used for internal purposes

PDE.dom=[-1 1];

% % This transdorms PDE coefficients to account for the PDE transform from physical to
% % computational domain, and takes care of the transformaton of the boundary
% % condition derivatives for n2 states 
% % 

 sp=0.5*(b-a)*s+0.5*(b+a);
 if isfield(PDE,'A0')
 PDE.A0=subs(PDE.A0,s,sp);
 end
 if isfield(PDE,'A1')
 PDE.A1=2*PDE.A1/(b-a);
 PDE.A1=subs(PDE.A1,s,sp);
 end
 if isfield(PDE,'A2')
 PDE.A2=4*PDE.A2/(b-a)^2;
 PDE.A2=subs(PDE.A2,s,sp);
 end
 if isfield(PDE,'B22')
 PDE.B22=subs(PDE.B22,s,sp);
 end
 if isfield(PDE,'B21')
 PDE.B21=subs(PDE.B21,s,sp);
 end
 if isfield(PDE,'E')
 PDE.E=subs(PDE.E,s,sp);
 end
 if isfield(PDE,'B')
 for i=(nbc-2*PDE.n2+1):nbc
 PDE.B(:,i)=PDE.B(:,i)*2/(b-a);
 end
 end
 if isfield(PDE,'E0')
     if (size(PDE.E0,2)==nbc)
 for i=(nbc-2*PDE.n2+1):nbc
 PDE.E0(:,i)=PDE.E0(:,i)*2/(b-a);
 end
 end
 end
 
 
 
 
  



