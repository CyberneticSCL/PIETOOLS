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
         
    Neigs = size(Atotal,1);
    lam = eigs(Atotal, Neigs);
    
    lamR = zeros(length(lam),1);
    lamI = zeros(length(lam),1);
    
    for i = 1:length(lam)
        lamR(i) = real(lam(i));
        lamI(i) = imag(lam(i));
    end

    
    switch opts.intScheme
        case 1
    % BDF scheme

    % StabSirPar corresponds to a stability region paramerer (radius and
    % center of a circle) for Norder=1, 2, and an estimte stability region
    % property for Norder=3, 4
    
      StabCirPar=[1, 2, 3.4, 5.4];
      StabImPar=[0, 0, 2, 5];

     Stable_circle=false; 
     Stable_im=false; 
     if (opts.Norder<=2) 
         Stable_im=true;
     end
     if (abs(lam*opts.dt-StabCirPar(opts.Norder))>=StabCirPar(opts.Norder)) 
         Stable_circle=true;
     end
     if(opts.Norder>=3)
         if (lamI*opts.dt==0 | abs(lamR*opts.dt)>=1e-2)
             Stable_im=true;
         elseif (abs(lamI*opts.dt)>=StabImPar(opts.Norder)) 
             Stable_im=true;
         end
     end

     Stable=[Stable_circle, Stable_im];

     ampl=0;
     
    if (Stable) 
        disp('Time integration scheme is numerically stable for the given problem.');
        disp('Any observed instabilities must be physical.');
        if(opts.Norder>=3)
            disp('WARNING: Stability check is imprecise for high-order schemes.')
        end
    else
    if (~Stable_circle) 
        ampl=max(2*StabCirPar(opts.Norder)*cos(angle(lam))./abs(lam*opts.dt));
    end
    if (~Stable_im) 
        ampl=max(max(StabImPar(opts.Norder)./abs(lamI*opts.dt)),ampl);
    end
        if (~Stable_circle)
        disp('Time integration scheme is numerically unstable for the given problem.');
        else
        disp('Time integration scheme may be numerically unstable for the given problem.');
        end

        if (opts.Norder<=2)
        formatSpec= 'Try increasing time step to %4.3f\n';
        else
        formatSpec= 'Try increasing time step to (an estimate of) %4.3f\n';
        end
        fprintf(formatSpec, opts.dt*ampl);
        if (opts.Norder~=1)
        disp('or decreasing an order of the scheme (opts.Norder).');
        end
        if (opts.Norder>=3)
        disp('WARNING: time step estimates may be inaccurate for high-order schemes.');
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
