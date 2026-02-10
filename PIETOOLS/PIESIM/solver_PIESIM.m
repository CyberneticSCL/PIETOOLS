<<<<<<< Updated upstream:PIETOOLS/PIESIM/solver_PIESIM.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version of PIETOOLS for 1D and 2D problems: 
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu

% This is the main driver for the PIESIM code if a numerical solution of PDE/ODE problem is the only task required.
% It can be called instead of PIETOOLS_PDE in this case. 
% Examples can be drawn from 'examples_pde_library_PIESIM.m' for 1D or
% 'examples_pde_library_PIESIM_2D.m' for 2D.
%
% This routine performs:
%
% 1) Conversion of a PDE (or PDE+ODE) problem to a PIE representation
% 2) Spatial discretization of the PIE operators (with Chebyshev polynomial approximation - high-order)
% 3) Temporal discretization (up to 4th order - user-defined parameter) and time advancement of the resulting ODE
% system 
% 4) Plotting and output of results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date,
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% YP - Added 2D functionality - 4_16_2024

% YP 12/21/2025 - (1) moved options before example setup, (2) pass options to examples for
% generating data array disturbance, (3) moved uinput.ifexact to
% opts.ifexact

clear all;
clc;
close all;
format long;

%--------------------------------------------------------------
% SETUP OF THE SIMULATIONS USER INPUT BEGINS
%--------------------------------------------------------------

% This script can be used to simulate the following linear systems:
% PDE/ODE or PDE/DDE

% Set up the dimenion of the problem (dim=1 for 1D problems or dim=2 for 2D problems).

dim=1;

% Set up options for the simulations

% If exact solution is known (see examles) and is desired to be provided for testing,
% select opts.ifexact=true.
% NOTE: Only choose opts.ifexact=true if exact solution for all the states is
%       available. Otherwise, choose opts.ifexact=false.

opts.ifexact=true;
%-----------------------------------------------------------
opts.plot='yes';
opts.ploteig='yes';

% Here we define parameters related to simulation.

% Input N - the Chebyshev polynomial discretization order of the
%            distributed states
opts.N=16;

%-----------------------------------------------------------
% Input the desired final time of the solution
opts.tf=1.6;

%-----------------------------------------------------------
% Choose temporal integration scheme 
%  opts.intScheme = 1 - Backward Differentiation Formula (BDF) 
%  opts.intScheme = 2 - Analytical integration in symbolic form 
% Note: opts.intScheme=2 will only work if the boundary and forcing inputs
%       are simple integrable functions of time, and the matrix 
%       Atotal=inv(T)*A is diagonalizable. An error will be issued if 
%       matrix is not diagonalizable, and a default integration scheme 
%       given by opts.intScheme = 1 (BDF) of order 2 (opts.Norder=2) will 
%       be executed.

% Choose opts.intScheme=1 if a
%  temporal history of solution is required (solution history is not
%  provided with analytical integration)
% if opts.intScheme = 1
%  Choose the order of numerical time integration scheme (Norder). A time
%  integration scheme available is Backward Differentiation Formula (BDF).
%  opts.Norder = 1, 2, 3 or 4 can be chosen. Lower order yield more robust
%  schemes, and higher order more accurate schemes. Also input the desired
%  time step (opts.dt)

opts.intScheme=1;

if (opts.intScheme==1)
    opts.Norder = 2;
    opts.dt=0.02;
end


% Set up example to run from the examples library (additional examples can
% be added to the library by the user).
%------------------------------------------------------------------------------
% For 1D problems: example = xxx (between 1 and 41) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_1D.m'
% For 2D problems: example = xxx (between 1 and 19) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_2D.m'
%------------------------------------------------------------------------------

    example=7;

    if (dim==1)
    if (example<1|example>41)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_1D(example);
    else   % dim=2
    if (example<1|example>19)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_2D(example);
    end

%--------------------------------------------------------------------------
%   USER INPUT ENDS
%--------------------------------------------------------------------------

% If disturbances are used while supplying PIE input to PIESIM, 
% user needs to ensure that finite-dimensional distubrances are
% entered before infinite-dimensional disturbances for a proper
% functionality. Disturbances are not reordered automatically with the PIE
% structure as they are with the PDE structure.

% PIE = convert_PIETOOLS_PDE(PDE,'pie');

if exist('PIE','var')
    solution = PIESIM(PIE,opts,uinput,PDE.n.n_pde);
elseif exist('DDE','var')
    solution=PIESIM(DDE,opts,uinput);
else
    solution = PIESIM(PDE,opts,uinput);
end






=======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver_PIESIM.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PIESIM Version of PIETOOLS for 1D and 2D problems: 
% For support, contact Y. Peet, Arizona State University at ypeet@asu.edu

% This is the main driver for the PIESIM code if a numerical solution of PDE/ODE problem is the only task required.
% It can be called instead of PIETOOLS_PDE in this case. 
% Examples can be drawn from 'examples_pde_library_PIESIM.m' for 1D or
% 'examples_pde_library_PIESIM_2D.m' for 2D.
%
% This routine performs:
%
% 1) Conversion of a PDE (or PDE+ODE) problem to a PIE representation
% 2) Spatial discretization of the PIE operators (with Chebyshev polynomial approximation - high-order)
% 3) Temporal discretization (up to 4th order - user-defined parameter) and time advancement of the resulting ODE
% system 
% 4) Plotting and output of results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date,
% authorship, and a brief description of modifications
%
% Initial coding YP  - 5_29_2021
% YP - Added 2D functionality - 4_16_2024

% YP 12/21/2025 - (1) moved options before example setup, (2) pass options to examples for
% generating data array disturbance, (3) moved uinput.ifexact to
% opts.ifexact

clear all;
clc
close all;
format long;

%--------------------------------------------------------------
% SETUP OF THE SIMULATIONS USER INPUT BEGINS
%--------------------------------------------------------------

% This script can be used to simulate the following linear systems:
% PDE/ODE or PDE/DDE

% Set up the dimenion of the problem (dim=1 for 1D problems or dim=2 for 2D problems).

dim=2;

% Set up options for the simulations

% If exact solution is known (see examles) and is desired to be provided for testing,
% select opts.ifexact=true.
% NOTE: Only choose opts.ifexact=true if exact solution for all the states is
%       available. Otherwise, choose opts.ifexact=false.

opts.ifexact=true;
%-----------------------------------------------------------
opts.plot='yes';
opts.ploteig='yes';

% Further custom-made post-processing options
opts.movie='no';
opts.error='no';

% Here we define parameters related to simulation.

% Input opts.N - the Chebyshev polynomial discretization order of the
%            distributed states 

% For 2D: if different order of discretization is required in different directions,
% enter opts.N as an array, such as opts.N=[16,32]

opts.N=16;

%-----------------------------------------------------------
% Input the desired final time of the solution
opts.tf=0.1;

%-----------------------------------------------------------
% Choose temporal integration scheme 
%  opts.intScheme = 1 - Backward Differentiation Formula (BDF) 
%  opts.intScheme = 2 - Analytical integration in symbolic form 
% Note: opts.intScheme=2 will only work if the boundary and forcing inputs
%       are simple integrable functions of time, and the matrix 
%       Atotal=inv(T)*A is diagonalizable. An error will be issued if 
%       matrix is not diagonalizable, and a default integration scheme 
%       given by opts.intScheme = 1 (BDF) of order 2 (opts.Norder=2) will 
%       be executed.

% Choose opts.intScheme=1 if a
%  temporal history of solution is required (solution history is not
%  provided with analytical integration)
% if opts.intScheme = 1
%  Choose the order of numerical time integration scheme (Norder). A time
%  integration scheme available is Backward Differentiation Formula (BDF).
%  opts.Norder = 1, 2, 3 or 4 can be chosen. Lower order yield more robust
%  schemes, and higher order more accurate schemes. Also input the desired
%  time step (opts.dt)

opts.intScheme=1;

if (opts.intScheme==1)
    opts.Norder = 2;
    opts.dt=0.02;
end


% Set up example to run from the examples library (additional examples can
% be added to the library by the user).
%------------------------------------------------------------------------------
% For 1D problems: example = xxx (between 1 and 41) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_1D.m'
% For 2D problems: example = xxx (between 1 and 31) to correspond to an Example number in
% the 'examples_pde_library_PIESIM_2D.m'
%------------------------------------------------------------------------------

    example=26;

    if (dim==1)
    if (example<1|example>41)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_1D(example);
    else   % dim=2
    if (example<1|example>32)
        disp('Warning: Example number is outside of the range. Defaulting to example=1');
        example=1;
    end
    [PDE,uinput]=examples_pde_library_PIESIM_2D(example);
    end

%--------------------------------------------------------------------------
%   USER INPUT ENDS
%--------------------------------------------------------------------------

% If disturbances are used while supplying PIE input to PIESIM, 
% user needs to ensure that finite-dimensional distubrances are
% entered before infinite-dimensional disturbances for a proper
% functionality. Disturbances are not reordered automatically with the PIE
% structure as they are with the PDE structure.

% PIE = convert_PIETOOLS_PDE(PDE,'pie');

if exist('PIE','var')
    solution = PIESIM(PIE,opts,uinput,PDE.n.n_pde);
elseif exist('DDE','var')
    solution=PIESIM(DDE,opts,uinput);
else
    [solution,grid] = PIESIM(PDE,opts,uinput);
end


% 
% % Make a movie for a 2D case (if opts.movie='yes')
% 
if strcmp(opts.movie,'yes') && isfield(solution,'timedep') && dim==2
figure;
X=grid.phys{1};
Y=grid.phys{2};
t=solution.timedep.dtime;
%Z= solution.timedep.pde{3}(:,:,1,1)';
Z=solution.timedep.observed{4}(:,:,1,1)';

% % Fix color scale
minZ = min(Z(:));
maxZ = max(Z(:));
%h = imagesc(X, Y, Z);  % initial frame
%h=contour(X, Y, Z, 100);
h = surf(X, Y, Z,'EdgeColor','none');  % initial frame
axis xy 
axis([min(X) max(X) min(Y) max(Y)])
axis equal
shading interp;  % smooth colors
view(2)          % top-down view
colormap(jet)
colorbar
caxis([minZ maxZ])

for k = 1:length(t)
    %Z= solution.timedep.pde{3}(:,:,1,k)';
    Z=solution.timedep.observed{4}(:,:,1,k)';
    set(h, 'ZData', Z);  % update surface
 %   set(h, 'CData', Z);  % update image
  %  h=contour(X, Y, Z, 100);
    title(['t = ', num2str(t(k))])
    drawnow
end


% % Movie with imagesc and interpolation - only works for uniform grids
% 
% % % Time vector
% 
% t = solution.timedep.dtime;
% 
% % --- Create a fine Cartesian grid for super-resolution ---
% nx = 1000;   % increase for higher smoothness
% ny = 1000;
% 
% xvec = unique(grid.phys{1});
% yvec = unique(grid.phys{2});
% 
% [Xq, Yq] = meshgrid(linspace(min(xvec), max(xvec), nx), ...
%                     linspace(min(yvec), max(yvec), ny));
% 
% % --- Interpolate first frame onto fine grid ---
% %Z = solution.timedep.pde{3}(:,:,1,1);   % [Ny,Nx]
% 
% Z=solution.timedep.observed{4}(:,:,1,1);
% 
% %     % Optional: flip rows if vortex appears mirrored
%       Z = flipud(Z);   % uncomment only if needed
%       Z = fliplr(Z);   % uncomment only if needed
% 
% Zq = interp2(xvec, yvec, Z, Xq, Yq, 'spline');
% 
% % --- Plot using imagesc ---
% figure
% h = imagesc(linspace(min(xvec), max(xvec), nx), ...
%             linspace(min(yvec), max(yvec), ny), Zq);
% axis xy
% axis equal
% colormap(jet)
% shading interp
% colorbar
% 
% % --- Fixed color scale ---
% %Zall = solution.timedep.pde{3}(:,:,1,:);
% Zall=solution.timedep.observed{4}(:,:,1,1);
% caxis([min(Zall(:)) max(Zall(:))])
% 
% % --- Animation loop ---
% for k = 1:length(t)
% %    Z = solution.timedep.pde{3}(:,:,1,k);       % coarse data
%     Z=solution.timedep.observed{4}(:,:,1,k);
% 
% %     % Optional: flip rows if vortex appears mirrored
%       Z = flipud(Z);   % uncomment only if needed
%       Z = fliplr(Z);   % uncomment only if needed
% % 
%     Zq = interp2(xvec, yvec, Z, Xq, Yq, 'spline');  % super-resolve
%     set(h, 'CData', Zq)
%     title(['t = ', num2str(t(k))])
%     drawnow
% end
% 
 end % opts.movie
% 

if strcmp(opts.error,'yes') 

% Compute L2 and Hinfty errors and plot their time evolution
% 
syms sx sy;

Nsteps=size(solution.timedep.dtime,2);

exact_output=diff(uinput.exact(2),sx)-diff(uinput.exact(1),sy);

area=4; % Area of the computattional grid [-1,1]^2

            if (length(opts.N)==1)
                opts.N(2)=opts.N(1);
            end

for n=1:1
       for ntime=1:Nsteps;
         tt=solution.timedep.dtime(ntime);

         % Evaluate exact solution
            exsol_numgrid_time = subs(subs(uinput.exact(n),sx,grid.phys{1}),sy,grid.phys{2}');
            exsol_numgrid = double(subs(exsol_numgrid_time,tt));

            exout_numgrid_time = subs(subs(exact_output,sx,grid.phys{1}),sy,grid.phys{2}');
            exout_numgrid = double(subs(exout_numgrid_time,tt));

         % Evaluate error
            error_state=exsol_numgrid-solution.timedep.pde{3}(:,:,n,ntime);
            error_out=exout_numgrid-solution.timedep.observed{4}(:,:,1,ntime);

        % Hinfty error
            Hinfty_error_state(ntime)=max(abs(error_state),[],'all');
            Hinfty_error_out(ntime)=max(abs(error_out),[],'all');

            % L2 error
            L2_error_state(ntime)=sqrt(sum(error_state.^2,'all')/prod(opts.N+1));
            L2_error_out(ntime)=sqrt(sum(error_out.^2,'all')/prod(opts.N+1));

            % Compute using Clenshaw-Curtis integration and compare
            [~, wx] = clencurt(opts.N(1));   % x ∈ [a,b], wx size 1×(N(1)+1)
            [~, wy] = clencurt(opts.N(2));   % y ∈ [c,d], wy size 1×(N(2)+1)
            W = wx(:) * wy(:).';   % (N(1)+1) × (N(2)+1)
            L2_error_state_CC(ntime)= sqrt(wx * error_state.^2 * wy.'/area);
            L2_error_out_CC(ntime)= sqrt(wx * error_out.^2 * wy.'/area);


       end % for ntime

       figure;
       plot(solution.timedep.dtime,L2_error_state,'b','Linewidth',2);
       hold on;
        plot(solution.timedep.dtime,L2_error_out,'bo');
        hold on;
       plot(solution.timedep.dtime,L2_error_state_CC,'r','Linewidth',2);  
        hold on;
        plot(solution.timedep.dtime,L2_error_out_CC,'ro');
       xlabel('t');ylabel('error');  
       title('Plot of errors versus time');
       ax = gca;   ax.FontSize = 22;

end % opts.error

end

>>>>>>> Stashed changes:PIESIM/solver_PIESIM.m
