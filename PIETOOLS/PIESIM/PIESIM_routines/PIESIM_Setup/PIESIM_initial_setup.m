%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_initial_setup.m    PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routines performs the following operations:
% 1) Finds time-derivatives of the user-defined boundary inputs and disturbances and
% adds them to the field 'uinput'
% 2) Transforms initial conditions from primary states to fundamental states
% via differentiation
%
% Input:
% 1) uinput - user-defined boundary inputs and disturbances
% 2) psize - size of the problem. Includes nu, nw, no, nf, N and n
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
% YP 12/31/2025 - changed disturbance and control inputs specifications to function handles (for
% both symbolic and double)
% YP 1/22/2026 - added differentiation of disturbances and control inputs
% to support inhomogeneous boudnary inputs in 2D

function uinput=PIESIM_initial_setup(uinput,psize,type)
syms st sx sy;

% Setup function handles to boundary inputs and their temporal derivatives

 if isfield(uinput,'w')
     for k = 1:numel(uinput.w)
         % Symbolic disturbances
        if isa(uinput.w{k}, 'sym')
        uinput.w{k}=matlabFunction(uinput.w{k},'Vars', st);
        uinput.wdot{k}=matlabFunction(diff(uinput.w{k},st),'Vars', st);
        else
            % Double type distubrances of size 1 (constant in time) 
            if (length(uinput.w{k})==1)
                uinput.w{k}=@(t) uinput.w{k};
                uinput.wdot{k}=@(t) 0;
            else
            %  Double type distubrances of size other than 1-
           %   Build a cubic spline from the disturbance data
        pp=spline(uinput.w{k}(1,:),uinput.w{k}(2,:));
        tmin=min(uinput.w{k}(1,:));
        tmax=max(uinput.w{k}(1,:));
        uinput.w{k} = @(t) ppval(pp,t).* (t >= tmin & t<= tmax);
        uinput.wdot{k} = @(t) ppval(fnder(pp),t).*(t >= tmin & t<= tmax);
            end
        end


        end % for k
     end % isfield(uinput,'w')

      if isfield(uinput,'u')
     for k = 1:numel(uinput.u)
         % Symbolic disturbances
        if isa(uinput.u{k}, 'sym')
        uinput.u{k}=matlabFunction(uinput.u{k},'Vars', st);
        uinput.udot{k}=matlabFunction(diff(uinput.u{k},st),'Vars', st);
        else
            % Double type distubrances of size 1 (constant in time) 
            if (length(uinput.u{k})==1)
                uinput.u{k}=@(t) uinput.u{k};
                uinput.udot{k}=@(t) 0;
            else
            %  Double type distubrances of size other than 1-
           %   Build a cubic spline from the disturbance data
        pp=spline(uinput.u{k}(1,:),uinput.u{k}(2,:));
        tmin=min(uinput.u{k}(1,:));
        tmax=max(uinput.u{k}(1,:));
        uinput.u{k} = @(t) ppval(pp,t).* (t >= tmin & t<= tmax);
        uinput.udot{k} = @(t) ppval(fnder(pp),t).*(t >= tmin & t<= tmax);
            end
        end

        end % for
     end % isfield(uinput,'u')

% Compute initial conditions on the fundamental states from initial
% condition on the primary states - only for PDE and DDE problems

if ~strcmp(type,'PIE') 

if (psize.dim==1)

    % 1D case
ns=sum(psize.n,'all');

psize_aux=[1 psize.n];
nsum=cumsum(psize.n);
nsump1=cumsum(psize_aux);

% Degree of smoothness 

 for i=1:length(psize.n)
 p(nsump1(i):nsum(i))=i-1;
 end

     for i=1:ns
      uinput.ic.PIE(i) =  diff(uinput.ic.PDE(i),sx,p(i));
     end

else

    % 2D case
% Use x_tab rather than psize to determine degree of smootness
        nsx=sum(psize.nx,'all');
        nsy=sum(psize.ny,'all');
        ns2d=sum(psize.n,'all');
        ns = nsx+nsy+ns2d;

for i=1:ns
    ii=psize.no+i;
    ic_PDE_xder =  diff(uinput.ic.PDE(i),sx,psize.x_tab(ii,end-1));
    uinput.ic.PIE(i) =  diff(ic_PDE_xder,sy,psize.x_tab(ii,end));
end

% Differentiate disturbances and control inputs if needed

for i=1:psize.nw
    w_xder =  diff(uinput.wspace{i},sx,psize.w_tab(i,end-1));
    uinput.wspace{i} =  diff(w_xder,sy,psize.w_tab(i,end));
end

for i=1:psize.nu
    u_xder =  diff(uinput.uspace{i},sx,psize.u_tab(i,end-1));
    uinput.uspace{i} =  diff(u_xder,sy,psize.u_tab(i,end));
end

end


end


