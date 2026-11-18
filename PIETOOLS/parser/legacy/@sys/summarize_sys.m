function summarize_sys(obj,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that displays a summary of sys object explains number of
% states and equations
% Input: 
% obj - sys class object
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
if strcmp(obj.type,'pde')
    statelist = obj.states;
    ode_list = statelist(find(strcmp(statelist.type,'ode')));
    n_ode = sum(ode_list.veclength);
    pde_list = statelist(find(strcmp(statelist.type,'pde')));
    n_pde = sum(pde_list.veclength);
    out_list = statelist(find(strcmp(statelist.type,'out')));
    n_out = sum(out_list.veclength);
    in_list = statelist(find(strcmp(statelist.type,'in')));
    n_in = sum(in_list.veclength);

    disp('----- Summary of identified states/dependent variables -----');
    fprintf('Number of ODE states: %d\n',n_ode);
    fprintf('Number of PDE states: %d\n',n_pde);
    fprintf('Number of regulated outputs: %d\n',n_out-sum(obj.ObservedOutputs));
    fprintf('Number of Observed outputs: %d\n',sum(obj.ObservedOutputs));
    fprintf('Number of Disturbance inputs: %d\n',n_in-sum(obj.ControlledInputs));
    fprintf('Number of Control inputs: %d\n',sum(obj.ControlledInputs));
    
    eqnType = varargin{1};
    
    disp('----- Summary of identified equations -----');
    fprintf('Number of ODE dynamics: %d\n',sum(strcmp(eqnType,'ode')));
    fprintf('Number of PDE dynamics: %d\n',sum(strcmp(eqnType,'pde')));
    fprintf('Number of output equations: %d\n',sum(strcmp(eqnType,'out')));
    fprintf('Number of Boundary conditions: %d\n',sum(strcmp(eqnType,'bc')));
end

end