%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_int_symbolic.m    PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs time-advancement of a discretized PIE using a convolution integral in
% a symbolic form
% Convolution integral is an exact solution in time of a
% discretized in space PIE (similar to an analytical solution of
% system of ODEs)
% Provides the most accurate solution but slow for large N (due to the use
% of symbolic integration)
% For a_t = Atotal*a + inv(T)*f(t), Atotal = inv(T)*A = V*D*inv(V)
% diagonalizable, f = B1*w - Tw*w_t + B2*u - Tu*u_t:
%   a(tf) = V*exp(D*tf)*inv(V)*a(0)
%           + V*int_0^tf exp(D*(tf-s))*inv(V)*inv(T)*f(s) ds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Inputs:
% psize - size of the PIE problem: nw, nu, no 
% opts - options for temporal scheme parameters
%    ('struct' with field tf, the final time)
% uinput   - user-defined boundary inputs
%    ('struct' with fields w, wdot, u, udot, the temporal parts of the
%    inputs, indexed as arrays and symbolic in st)
% coeff - Chebyshev coefficients of initial conditions and forcing terms, if any
%    ('struct' with fields acheb_f0, w, u)
% Dop - discretized PIE operators
%    ('struct' with fields Atotal, V, D (Atotal = V*D*inv(V)), Tcheb_inv,
%    B1cheb, B2cheb, Twcheb, Tucheb)
%
% % Output:
% solcoeff - contains:
% 1) solcoeff.final - Chebyshev coefficients of the final
% solution, ODE+PDE states
% 2) solcoeff.tf - final time of the solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_int_symbolic
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 5_26_2021
% MMP, 11/19/2026: Release comment pass (comment-only).

function solcoeff=PIESIM_int_symbolic(psize, opts, uinput, coeff, Dop)
        % st: time variable of the inputs; tt: upper limit of the
        % convolution integrals
        syms tt intsym st;

        nw=psize.nw;
        nu=psize.nu;
        no=psize.no;        
        tf=opts.tf;
        
        % Discretized operators, eigendecomposition and initial
        % coefficients
        Twcheb=Dop.Twcheb;
        Tucheb=Dop.Tucheb;
        B1cheb=Dop.B1cheb;
        B2cheb=Dop.B2cheb;
        Acheb=Dop.Acheb;
        Tcheb_inv=Dop.Tcheb_inv;
        Atotal=Dop.Atotal;
        V=Dop.V;
        D=Dop.D;
        Nsize=size(Atotal,1);
        acheb_f0=coeff.acheb_f0;
        
           % Row j of inhom: B1*coeff.w*wvec - Tw*coeff.w*wdotvec
           % + B2*coeff.u*uvec - Tu*coeff.u*udotvec, with
           % wvec(i) = int_0^tf exp(lambda_j*(tf-s))*w_i(s) ds, lambda_j =
           % D(j,j) (likewise wdotvec, uvec, udotvec)
           inhom=zeros(Nsize,Nsize);
           for j=1:Nsize
             lambda=D(j,j);
             if (nw>0)
                 
    % Contribution to inhomogeneous term due to boundary disturbances  
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             for i=1:size(uinput.w,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.w(i),st,[0,tt]);
             wvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary disturbances             
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             for i=1:size(uinput.wdot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.wdot(i),st,[0,tt]);
             wdotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B1cheb)==0&any(B1cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B1cheb*coeff.w*wvec).';
             end
             if (isempty(Twcheb)==0&any(Twcheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Twcheb*coeff.w*wdotvec).';
             end
             end % nw>0
             
              if (nu>0)
                 
    % Contribution to inhomogeneous term due to boundary inputs  
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             for i=1:size(uinput.u,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.u(i),st,[0,tt]);
             uvec(i,1)=double(subs(intsym,tt,tf));
             end        
             end 
             
    % Contribution to inhomogeneous term due to time derivative of boundary inputs             
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             for i=1:size(uinput.udot,2);      
             intsym=int(exp(lambda*(tt-st))*uinput.udot(i),st,[0,tt]);
             udotvec(i,1)=double(subs(intsym,tt,tf));
             end 
             end
            
             
             if (isempty(B2cheb)==0&any(B2cheb,'all')~=0)
             inhom(j,:)=inhom(j,:)+(B2cheb*coeff.u*uvec).';
             end
             if (isempty(Tucheb)==0&any(Tucheb,'all')~=0)
             inhom(j,:)=inhom(j,:)-(Tucheb*coeff.u*udotvec).';
             end
             end % nu>0
  
           end    % for j loop
       
   % Adding all contributions to the time integral. Only considering the
   % columns of inhom array that are non-zero
           
       addition_symbolic=0;
       if (isempty(inhom)==0) 
           [row,col]=find(inhom);
           C=unique(col);
       for k=1:size(C,1)
        i=C(k);
       addition_symbolic=addition_symbolic+diag(inhom(:,i))*inv(V)*Tcheb_inv(:,i);
       end
       end

     
        % a(tf) = V*exp(D*tf)/V*a(0) + V*addition_symbolic
        exp_lambda=exp(diag(D*tf));
        solcoeff.final=V*diag(exp_lambda)/V*acheb_f0+V*addition_symbolic;
        solcoeff.tf=tf;
        