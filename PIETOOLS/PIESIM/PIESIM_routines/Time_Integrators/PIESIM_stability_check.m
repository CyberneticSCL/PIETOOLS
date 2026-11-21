%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_stability_check.m    PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs check of a numerical stability of time integration and
% outputs suggestions if unstable
% BDF of order k = opts.Norder: the boundary locus
% z(t) = sum_{j=0}^{k} alpha_j*exp(-i*j*t), t in [0,2*pi], encloses the
% bounded region where the scheme is unstable. An eigenvalue lambda of
% Atotal with lambda*dt inside it, |lambda*dt| < |z| on the ray of
% angle(lambda*dt), gives the amplification max |z|/|lambda*dt| > 1
% (unstable; borderline within 0.05 of 1). Symbolic integration: reports
% whether some eigenvalue has a positive real part.
%
% Inputs:
% opts - options for temporal scheme parameters
%    ('struct' with fields intScheme, Norder, dt, ploteig; ploteig = 'yes'
%    plots the eigenvalues times dt)
% Atotal - matrix operator for the discrete linear system x_t=A_total x
% (achieved upon spatial discretization)
%
% Output:
% none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_stability_check
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
% Initial coding YP  - 6_16_2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.

function PIESIM_stability_check(opts, Atotal);
         
 %   Neigs = size(Atotal,1);
 %   lam = eigs(Atotal, Neigs);

       % Tolerances: borderline band of the amplification, real-part sign
       % test, angle window (rad)
       tol=0.05;
       tol_real=1e-5;
       tol_angle=0.1;

        % Eigenvalues of Atotal
        lam=eig(Atotal);
        Neig=length(lam);
        lamR = real(lam);
        lamI = imag(lam);

      if(strcmp(opts.ploteig,'yes'))
    % Plot eigenvalues
    
    figure("Name","Eigenplot, PIESIM")
    plot(lamR*opts.dt, lamI*opts.dt, 'x','MarkerSize',12,'linewidth',2)
    xlabel("\lambda_{real}dt")
    ylabel("\lambda_{imag}dt")
    xline(0,'linewidth',2);
    yline(0,'linewidth',2);

    ax = gca;
    ax.FontSize = 24;
    H=gca;
    H.LineWidth=3;
    end
    
    switch opts.intScheme
        case 1
    % BDF scheme

    % Boundary locus z(t) of the BDF scheme of order opts.Norder
    t=linspace(0,2*pi);
    switch(opts.Norder)

    case 1
    % First order
    bdf1=[1 -1];
    z=(bdf1(1)*exp(i.*t)+bdf1(2))./exp(i.*t);

    case 2
    % Second order
    bdf2=[3/2 -2 1/2];
    z=(bdf2(1)*exp(2.*i.*t)+bdf2(2)*exp(i.*t)+bdf2(3))./exp(2.*i.*t);

    case 3
    % Third order
    bdf3=[11/6 -3 3/2 -1/3];
    z=(bdf3(1)*exp(3.*i.*t)+bdf3(2)*exp(2.*i.*t)+bdf3(3)*exp(i.*t)+bdf3(4))./exp(3.*i.*t);

    case 4
    % Fourth order
    bdf4=[25/12 -4 3 -4/3 1/4];
    z=(bdf4(1)*exp(4.*i.*t)+bdf4(2)*exp(3.*i.*t)+bdf4(3)*exp(2.*i.*t)+bdf4(4)*exp(i.*t)+bdf4(5))./exp(4.*i.*t); 
    end

    % Per eigenvalue: the locus point of maximal modulus among those within
    % tol_angle of the closest angle to lambda*dt; am(k) = |z|/|lambda*dt|
    % if z and lambda*dt have real parts of one sign (product > tol_real)
    % and imaginary parts of one sign, else 0
    an=angle(z);

    for k=1:Neig
        lambda=lam(k)*opts.dt;

    diff=abs(an-angle(lambda));

    m=min(diff);
    m_max=m+tol_angle;

    ind_permit=find(diff<m_max);
    z_permit=z(ind_permit);

    [zmax,indmax_inter]=max(abs(z_permit));
    indmax=ind_permit(indmax_inter);

    if (real(lambda)*real(z(indmax))>tol_real & imag(lambda)*imag(z(indmax))>=0)
    am(k)=abs(z(indmax))/abs(lambda);
    else
    am(k)=0;
    end

    end

    % Unstable if ampl > 1; borderline if |ampl-1| < tol
    ampl=max(am);
    closeflag=0;

    if(ampl>1) 
        Stable=0;
        if(ampl<1+tol)
            closeflag=1;
        end
    else
        Stable=1;
        if(ampl>1-tol)
            closeflag=1;
        end
    end

    % Report; the time step dt*(ampl+tol) moves lambda*dt out of the
    % unstable region
    if (Stable) 
        if (closeflag)
        disp('Time integration scheme is on the borderline of stability.');
        disp('Numerical simulations might be unstable.');
        formatSpec= 'If unstable, try increasing time step to %10.8f\n';
        fprintf(formatSpec, opts.dt*(ampl+tol));
        else
        disp('Time integration scheme is numerically stable for the given problem.');
        disp('Any observed instabilities must be physical.');
        end
   else
        % Unstable: a time step dt*(ampl+tol) below 1 is suggested, else the
        % instability is attributed to the problem
        if (closeflag)
        disp('Time integration scheme is on the border of instability.');
        disp('Numerical simulations may be unstable.');
        formatSpec= 'If unstable, try increasing time step to %10.8f\n';
        fprintf(formatSpec, opts.dt*(ampl+tol));
        else
        disp('Time integration scheme is numerically unstable for the given problem.');
        if opts.dt*(ampl+tol)<1
        formatSpec= 'Try increasing time step to %10.8f\n';
        fprintf(formatSpec, opts.dt*(ampl+tol));
        if (opts.Norder~=1)
        disp('or decreasing an order of the scheme (opts.Norder).');
        end
        else
        disp('This is likely because the underlying physical problem is unstable.');
        end
        end
     end

        case 2
    % Symbolic integration

    if (lamR<=0) 
        disp('There are no eigenvalues with positive real part. Symbolic integration will produce bounded results.');
    else
        disp('There is an eigenvalue with positive real part. Symbolic integration may produce unbounded results.');
    end

    end
