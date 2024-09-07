%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_discretize_ops_2D.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform discretization of PIE operators with Chebyshev methods  
%
% Inputs:
% 1) PIE - PIE structure of the problem
% 2) psize - size of the PIE problem: all variables defining the size of the PIE problem
%
%
% Outputs:
% 1) Dop - discrete PIE operators containing Chebyshev matrices for
% all PIE operators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 04_16_2024

function Dop=PIESIM_discretize_ops_2D(PIE, psize);

% disp('Setting up Chebyshev matrices for the PIE system');
% 
 Dop.D11cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D11,psize,0);
 Dop.D12cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D12,psize,0);
 Dop.D21cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D21,psize,0);
 Dop.D22cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.D22,psize,0);
 [Dop.Twcheb,Dop.Twcheb_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(PIE.Tw,psize,1);
 [Dop.Tucheb,Dop.Tucheb_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(PIE.Tu,psize,1);
 Dop.B1cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.B1,psize,1);
 Dop.B2cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.B2,psize,1);
 Dop.C1cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.C1,psize,2);
 Dop.C2cheb=PIESIM_fullPI2Mat_cheb_2D(PIE.C2,psize,2);
 Dop.Acheb=PIESIM_fullPI2Mat_cheb_2D(PIE.A,psize,3);
 [Mcheb, Dop.Mcheb_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(PIE.T,psize,3);
 if isfield(PIE,'T0') 
     [Mcheb0, Dop.Mcheb0_nonsquare]=PIESIM_fullPI2Mat_cheb_2D(PIE.T0,psize,3);
 end
%  
  Dop.Mcheb_inv=inv(Mcheb);
  Dop.Atotal=Dop.Mcheb_inv*Dop.Acheb;





