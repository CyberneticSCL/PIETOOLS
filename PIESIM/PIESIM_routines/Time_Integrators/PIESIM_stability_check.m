%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_stability_check.m    PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs check of a numerical stability of time integration and
% outputs suggestions if unstable
%
% Inputs:
% opts - options for temporal scheme parameters
% Atotal - matrix operator for the discrete linear system x_t=A_total x
% (achieved upon spatial discretization)
%
% Output:
% none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 6_16_2022

function PIESIM_stability_check(opts, Atotal);
         
 %   Neigs = size(Atotal,1);
 %   lam = eigs(Atotal, Neigs);

       tol=0.2;

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

    an=angle(z);

    for k=1:Neig
        lambda=lam(k)*opts.dt;

    diff=abs(an-angle(lambda));
    [m,indmin]=min(diff);
     
    if (real(lambda)*real(z(indmin))>0 & imag(lambda)*imag(z(indmin))>0)
    am(k)=abs(z(indmin))/abs(lambda);
    else
    am(k)=0;
    end
    end

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
        if (closeflag)
        disp('Time integration scheme is on the border of instability.');
        disp('Numerical simulations may be unstable.');
        formatSpec= 'If unstable, try increasing time step to %10.8f\n';
        fprintf(formatSpec, opts.dt*(ampl+tol));
        else
        disp('Time integration scheme is numerically unstable for the given problem.');
        formatSpec= 'Try increasing time step to %10.8f\n';
        fprintf(formatSpec, opts.dt*(ampl+tol));
        if (opts.Norder~=1)
        disp('or decreasing an order of the scheme (opts.Norder).');
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
