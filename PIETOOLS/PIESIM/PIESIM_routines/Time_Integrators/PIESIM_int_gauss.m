%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_gauss.m    PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE using  Gauss integration
%    of the analytical representation of the discretized
%    in space solution as a convolution in time (similar to a solution of a system of ODEs)
%
% Inputs:
% psize - size of the PIE problem: nw, nu, nf, nx 
% opts - options for temporal scheme parameters
% uinput   - user-defined boundary inputs
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
% Dop - discretized PIE operators
%
% Output:
% solcoeff - contains:
% 1) solcoeff.acheb_f - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_26_2021

function solcoeff=PIESIM_int_gauss(psize, opts, uinput, coeff, Dop)

nw=psize.nw;
nu=psize.nu;
nf=psize.nf;
nx=psize.nx;
tf=opts.tf;
Nint=opts.Nint;
Npgauss=opts.Npgauss;

Twcheb=Dop.Twcheb;
Tucheb=Dop.Tucheb;
B1cheb=Dop.B1cheb;
B2cheb=Dop.B2cheb;
Acheb=Dop.Acheb;
Mcheb_inv=Dop.Mcheb_inv;
Atotal=Dop.Atotal;
Nsize=size(Atotal,1);
acheb_f0=coeff.acheb_f0;
if (nf>0)
acheb_force=coeff.acheb_force;
tforce=uinput.tforce;
end


% %        
%    Calculate size of each interval 
     tint=tf/Nint;
%        
%            
           intsum(1:Nsize,Nint)=0;
           [zgll,w]=zwgll(Npgauss);
                 
%       
        for n=1:Nint
            
%           Uniform intervals are assumed
            t0 = (n-1)*tint;
            t1 = n*tint;
%           
% % % Mapping of points
           z=0.5*(t1-t0)*zgll+0.5*(t0+t1);
%           
           integral_gauss(1:Nsize,1)=0;
           for j=1:Npgauss+1 
            inhom(1:Nsize,1)=0;
            
            
            if (nw>0)
             % Contribution to inhomogeneous term due to boundary disturbances  
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
                 wvec(:,1)=double(subs(uinput.w(:),z(j)));
                 inhom=inhom+Mcheb_inv*B1cheb*wvec;
             end 
              % Contribution to inhomogeneous term due to time derivative of boundary disturbances  
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             wdotvec(:,1)=double(subs(uinput.wdot(:),z(j)));
             inhom=inhom-Mcheb_inv*Twcheb*wdotvec;
             end
            end
            
             if (nu>0)
             % Contribution to inhomogeneous term due to boundary inputs  
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
                 uvec(:,1)=double(subs(uinput.u(:),z(j)));
                 inhom=inhom+Mcheb_inv*B2cheb*uvec;
             end 
              % Contribution to inhomogeneous term due to time derivative of boundary inputs 
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             udotvec(:,1)=double(subs(uinput.udot(:),z(j)));
             inhom=inhom-Mcheb_inv*Tucheb*udotvec;
             end
            end
            

    % Contribution to inhomogeneous term due to distributed forcing 
    % (not necessarily polynomial)
            for i=1:nf
             inhom=inhom+Mcheb_inv*double(subs(tforce(i),z(j)))*acheb_force(:,i);
            end
            integral_gauss=integral_gauss+0.5*(t1-t0)*w(j)*expm(Atotal*(tf-z(j)))*inhom;
           end
           intsum(:,n)=integral_gauss;
           t0=t1;
        end
           integral_gauss=sum(intsum,2);
           
        solcoeff.acheb_f=expm(Atotal*tf)*acheb_f0+integral_gauss;
        solcoeff.tf=tf;
       
       
        
        
       