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
% 2) psize - size of the problem
% 3) PIE - PIE structure of the problem to redefine disturbances due to
% corner values if needed
% 4) type - of class ``char'' - type of the problem: 'PDE', 'DDE' or 'PIE'
%
% Output:
% 1) uinput - temporal derivatives of the user-defined boundary inputs and disturbances are added to the 
% field 'input' - these are needed for temporal integration of PIE equations
% 2) psize - can be modified to account for corner values of disturbances
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
% YP 2/6/2026 - added corner values for inhomogeneous boudnary inputs in 2D

function [uinput,psize]=PIESIM_initial_setup(uinput,psize,PIE,type)
syms st sx sy sym_array;

sym_array=[sx,sy];

% Account for multiplicity of disturbances
w_tab = repelem(PIE.w_tab, PIE.w_tab(:,2), 1);
u_tab = repelem(PIE.u_tab, PIE.u_tab(:,2), 1);

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
    ns = sum(psize.nx,'all')+sum(psize.ny,'all')+sum(psize.n,'all');

for i=1:ns
    ii=psize.no+i;
    ic_PDE_xder =  diff(uinput.ic.PDE(i),sx,psize.x_tab(ii,end-1));
    uinput.ic.PIE(i) =  diff(ic_PDE_xder,sy,psize.x_tab(ii,end));
end

% Compare the size of disturbances after PIE conversion to original size
% and renumber if needed

  nw_PDE=psize.nw;
  nw_PIE = sum(w_tab(:,2));
  idx_shift=nw_PIE-nw_PDE;

  if (nw_PIE>nw_PDE)
  uinput.wspace = [cell(1, idx_shift), uinput.wspace];
  uinput.w = [cell(1, idx_shift), uinput.w];

  % We also have to fill the first idx_shift rows with the corresponding
  % corner values
  for i=1:idx_shift
  pos_in_reordered_array = find(uinput.w_order == w_tab(i,1)); % w_tab corresponds to original PDE ordering
  idx_parent=pos_in_reordered_array+idx_shift;
  uinput.w{i}=uinput.w{idx_parent};
  uinput.wspace{i}=uinput.wspace{idx_parent};
  end
  end % if

  psize.nw=nw_PIE;
  psize.nw0=psize.nw0+idx_shift;

 
% Differentiate disturbances if needed (2D only)

for i=1:psize.nw
    w_xder =  diff(uinput.wspace{i},sx,w_tab(i,end-1));
    uinput.wspace{i} =  diff(w_xder,sy,w_tab(i,end));
end


for i=1:idx_shift
    % Do the substitution for corner values
  % If idx is not equal to i, input i must correspond to a corner value
  idx = find(w_tab(:,1)==w_tab(i,1),1,'last');
  % Determine which variables to use for corner values
  is_bval = logical(w_tab(idx,3:4)-w_tab(i,3:4));
  corner_value=double(subs(uinput.wspace{i},sym_array(is_bval),PIE.dom(is_bval,1)));
  uinput.wspace{i}=[];
  if (length(uinput.w{i}==1))
  uinput.w{i}=corner_value*uinput.w{i};
  else
  uinput.w{i}(1,:)=corner_value*uinput.w{i}(1,:);
  end
end

% Dealing with control inputs 

for i=1:psize.nu
    u_xder =  diff(uinput.uspace{i},sx,psize.u_tab(i,end-1));
    uinput.uspace{i} =  diff(u_xder,sy,psize.u_tab(i,end));
end

end % 2D case

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
            %  Double type disturbances of size other than 1-
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


end


