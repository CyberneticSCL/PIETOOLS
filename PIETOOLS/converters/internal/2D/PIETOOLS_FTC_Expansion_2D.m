function [Pop_f2x, Pop_b2x, bc_state_tab] = PIETOOLS_FTC_Expansion_2D(PDE)
% Construct PI operators Pop_f2x and Pop_b2x such that the PDE state x
% can be expressed in terms of the associated fundamental state xf and a
% set of core boundary values xb as
%   x = Pop_f2x * xf + Pop_b2x * xb;
% 
% INPUTS:
% - PDE:    A struct or pde_struct object defining a PDE in the terms
%           format (see also "@pde_struct/initialize").
%
% OUTPUTS:
% - Pop_f2x:    opvar2d object describing the map from the fundamental
%               state xf to the PDE state x defined by PDE.x_tab.
% - Pop_b2x:    opvar2d object describing the map from the core boundary
%               state xb to the PDE state x defined by PDE.x_tab.
% - bc_state_tab:   A Nx(2+2*nvars) array, defining the core boundary state. 
%               For j=1,...N, element
%               (j,2) should provide the size of the state component;
%               (j,3:2+nvars) should be binary indices, indicating for
%               each of the nvars variables whether the component j depends
%               on this variable. If the associated fundamental state
%               component depends on a variable s, but the boundary state
%               component does not, the boundary state component is
%               evaluated at the lower boundary s=a (for s\in[a,b]).
%               (j,3+nvars:end) should be integer values, indicating for
%               each of the nvars variables to what degree the state is
%               differentiated wrt this variable.
%
% NOTES:
% - For a state variables x(s1,s2), differentiable up to order i wrt s1 and
%   order j wrt s2, the associated fundamental state variable xf(s1,s2) is
%   defined as
%       xf(s1,s2) = \partial_{s1}^{i} \partial_{s2}^{j} x(s1,s2);
%
% - For a PDE state 
%       x = [ x0        ] = [ x0                ]
%           [ --------- ]   [ ----------------- ]
%           [ x1(s1)    ]   [ x1_{0}(s1)        ]
%           [ --------- ]   [ :                 ]
%           [ x2(s2)    ]   [ x1_{K1}(s1)       ]
%           [ --------- ]   [ ----------------- ]
%           [ x3(s1,s2) ]   [ x2_{0}(s2)        ]
%                           [ :                 ]
%                           [ x2_{M1}(s2)       ]
%                           [ ----------------- ]
%                           [ x3_{0,0}(s1,s2)   ]
%                           [ x3_{1,0}(s1,s2)   ]
%                           [ :                 ]
%                           [ x3_{N1,0}(s1,s2)  ]
%                           [ x3_{0,1}(s1,s2)   ]
%                           [ :                 ]
%                           [ x3_{N1,N2}(s1,s2) ],
%
%   where s1\in[a1,b2] and s2\in[a2,b2], and where xk_{i,j} is 
%   differentiable up to degree i with respect to s1, and degree j wrt s2, 
%   we define the associated core boundary state as
%
%  xb = [                                           x1_{1}(a1)          ]
%       [                                           x1_{2}(a1)          ]
%       [ \partial_{s1}                             x1_{2}(a1)          ]
%       [                                           x1_{3}(a1)          ]
%       [                                           :                   ]
%       [ \partial_{s1}^(K-1)                       x1_{K}(a1)          ]
%       [ ------------------------------------------------------------- ]
%       [                                           x2_{1}(a2)          ]
%       [                                           x2_{2}(a2)          ]
%       [ \partial_{s2}                             x2_{2}(a2)          ]
%       [                                           x2_{3}(a2)          ]
%       [                                           :                   ]
%       [ \partial_{s2}^(M-1)                       x2_{M}(a2)          ]
%       [ ------------------------------------------------------------- ]
%       [                                           x3_{1,1}(a1,a2)     ]
%       [                                           x3_{2,1}(a1,a2)     ]
%       [ \partial_{s1}                             x3_{2,1}(a1,a2)     ]
%       [                                           x3_{3,1}(a1,a2)     ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1-1}                      x3_{N1,1}(a1,a2)    ]
%       [                                           x3_{1,2}(a1,a2)     ]
%       [                      \partial_{s2}        x3_{1,2}(a1,a2)     ]
%       [                                           x3_{2,2}(a1,a2)     ]
%       [ \partial_{s1}                             x3_{2,2}(a1,a2)     ]
%       [                      \partial_{s2}        x3_{2,2}(a1,a2)     ]
%       [ \partial_{s1}        \partial_{s2}        x3_{2,2}(a1,a2)     ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1-1} \partial_{s2}^{N2-1} x3_{N1,N2}(a1,a2)   ]
%       [ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ]
%       [                                           x3_{0,1}(s1,a2)     ]
%       [ \partial_{s1}                             x3_{1,1}(s1,a2)     ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1}                        x3_{N1,1}(s1,a2)    ]
%       [                                           x3_{0,2}(s1,a2)     ]
%       [                      \partial_{s2}        x3_{0,2}(s1,a2)     ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1}   \partial_{s2}^{N2-1} x3_{N1,N2}(s1,a2)   ]
%       [ - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ]
%       [                                           x3_{1,0}(a1,s2)     ]
%       [                                           x3_{2,0}(a1,s2)     ]
%       [ \partial_{s1}                             x3_{2,0}(a1,s2)     ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1-1}                      x3_{N1,0}(a1,s2)    ]
%       [                     \partial_{s2}^{N2}    x3_{1,N2}(a1,s2)    ]
%       [                                           :                   ]
%       [ \partial_{s1}^{N1-1} \partial_{s2}^{N2}   x3_{N1,N2}(a1,s2)   ]
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_FTC_Expansion_2D
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/09/2022
% Add support for 1D case, DJ - 10/16/2024.
% DJ, 12/06/2024: Bugfix `ind2sub` for Matlab 2024b;


% Extract PDE information.
dom = PDE.dom;
vars = PDE.vars;
x_tab = PDE.x_tab;
nvars = size(vars,1);


% Establish how many state components there are that vary in each possible
% combination of variables.
dep_tab = x_tab(:,3:2+nvars);
diff_tab = x_tab(:,3+nvars:2*nvars+2);
ncomps = size(x_tab,1);
np_arr = x_tab(:,2);    % Number of state variables in each component
nnp_arr = cumsum([0;np_arr]);   % Row number in the full state associated to the first variable of each component
if nvars==0
    x_tab_rindcs_00 = true(length(np_arr),1);
    x_tab_rindcs_10 = false(length(np_arr),1);
    x_tab_rindcs_01 = false(length(np_arr),1);
    x_tab_rindcs_11 = false(length(np_arr),1);
elseif nvars==1
    x_tab_rindcs_00 = ~dep_tab;                 % Row indices in x_tab associated to finite dimensional states
    x_tab_rindcs_10 = ~x_tab_rindcs_00;         % Row indices in x_tab associated to states that vary just in s1
    x_tab_rindcs_01 = false(length(np_arr),1);
    x_tab_rindcs_11 = false(length(np_arr),1);
elseif nvars==2    
    x_tab_rindcs_00 = ~any(dep_tab,2);                 % Row indices in x_tab associated to finite dimensional states
    x_tab_rindcs_10 = dep_tab(:,1) & ~dep_tab(:,2);   % Row indices in x_tab associated to states that vary just in s1
    x_tab_rindcs_01 = ~dep_tab(:,1) & dep_tab(:,2);   % Row indices in x_tab associated to states that vary just in s2
    x_tab_rindcs_11 = dep_tab(:,1) & dep_tab(:,2);    % Row indices in x_tab associated to states that vary in s1 and s2
else
    error('At most 2 spatial variables are currently supported.')
end

% Set the dimension of the state.
np_00 = sum(x_tab(x_tab_rindcs_00,2));
np_10 = sum(x_tab(x_tab_rindcs_10,2));
np_01 = sum(x_tab(x_tab_rindcs_01,2));
np_11 = sum(x_tab(x_tab_rindcs_11,2));

if nvars<=1
    np_op = [np_00; np_10];
else
    np_op = [np_00; np_10; np_01; np_11];
end
nnp_op = cumsum([0;np_op]);

% Establish the number of core boundary values
nbc_tab = dep_tab.*diff_tab;        % For each order of differentiability, along each spatial direction, we need a BC

if ~any(x_tab_rindcs_10)
    nbc_10_00 = 0;
else
    nbc_10_00 = np_arr(x_tab_rindcs_10)'*(nbc_tab(x_tab_rindcs_10,1));      % Number of finite-dimensional core boundary states for x_{10}(s_1)
end
if ~any(x_tab_rindcs_01)
    nbc_01_00 = 0;
else
    nbc_01_00 = np_arr(x_tab_rindcs_01)'*(nbc_tab(x_tab_rindcs_01,2));      % Number of finite-dimensional core boundary states for x_{01}(s_2)
end
if ~any(x_tab_rindcs_11)
    nbc_11_00 = 0;
    nbc_11_10 = 0;
    nbc_11_01 = 0;
else
    nbc_11_00 = np_arr(x_tab_rindcs_11)'*(prod(nbc_tab(x_tab_rindcs_11,:),2));  % Number of finite-dimensional core boundary states for x_{11}(s_1,s_2)
    nbc_11_10 = np_arr(x_tab_rindcs_11)'*(nbc_tab(x_tab_rindcs_11,2));          % Number of core boundary states for x_{11}(s_1,s_2) that VARY in just s_1 (so are BCs along s_2)
    nbc_11_01 = np_arr(x_tab_rindcs_11)'*(nbc_tab(x_tab_rindcs_11,1));          % Number of core boundary states for x_{11}(s_1,s_2) that VARY in just s_2 (so are BCs along s_1)
end
% Set the dimension of the core boundary state.
if nvars<=1
    nbc_op = [nbc_10_00; 0];
else
    nbc_op = [nbc_10_00 + nbc_01_00 + nbc_11_00;
              nbc_11_10;
              nbc_11_01;
              0];
end

% Initalize a table to keep track of what each column of the BC operator 
% corresponds to.
% First column indicates a state component
% Second column indicates the size of this component
% Next nvar columns are logical indices indicating whether the boundary
% state varies in each of the nvars variables
% The remaining nvar columns are logical indices indicating whether the
% state is differentiated along each variable.
bc_state_tab = zeros(sum(nbc_op),2+2*nvars);

% Initalize the PI operators.
if nvars<=1
    opvar Pop_f2x Pop_b2x;
    Pop_f2x.dim = [np_op,np_op];    Pop_f2x.I = dom;
    Pop_f2x.var1 = vars(1);         Pop_f2x.var2 = vars(2);
    Pop_b2x.dim = [np_op,nbc_op];  Pop_b2x.I = dom;
    Pop_b2x.var1 = vars(1);         Pop_b2x.var2 = vars(2);

    intrr_params = cell(2,2);
    intrr_params{1,1} = {Pop_f2x.P};
    intrr_params{2,2} = {Pop_f2x.R.R0, Pop_f2x.R.R1, Pop_f2x.R.R2};
    
    bndry_params = cell(2,2);
    bndry_params{2,1} = {Pop_b2x.Q2};
else
    Pop_f2x = opvar2d([],[np_op,np_op],dom,vars);
    Pop_b2x = opvar2d([],[np_op,nbc_op],dom,vars);

    intrr_params = cell(4,4);
    intrr_params{1,1} = {Pop_f2x.R00};
    intrr_params{2,2} = Pop_f2x.Rxx;
    intrr_params{3,3} = Pop_f2x.Ryy;
    intrr_params{4,4} = Pop_f2x.R22;
    
    bndry_params = cell(4,4);
    bndry_params{2,1} = {Pop_b2x.Rx0};
    bndry_params{3,1} = {Pop_b2x.Ry0};
    bndry_params{4,1} = {Pop_b2x.R20};
    bndry_params{4,2} = Pop_b2x.R2x;
    bndry_params{4,3} = Pop_b2x.R2y;
end


nBCs_arr = zeros(1,2^nvars - 1);

for comp=1:ncomps
    
    % Establish which element of Rcell we're interested in
    % e.g. if the state component varies only in s1, we should look at Rxx
    nx_comp = np_arr(comp);
    x_rindcs = nnp_arr(comp)+1:nnp_arr(comp+1);   % Row indices in the full state associated to component comp
    Pop_rnum = (x_rindcs(1)>nnp_op(1:end-1) & x_rindcs(1)<=nnp_op(2:end));
    rindcs = x_rindcs - nnp_op(Pop_rnum);     % Row indices in the actual parameters, associated to the current component
    
    isvar = logical(dep_tab(comp,:));
    nvars_cc = sum(isvar);
    
    % Establish which parameter within Rkk we are interested in
    % e.g. if the state component is differentiable in s1 but not in s2, we
    % should look at R22{2,1}
    param_lin_sz = 3.^(0:nvars-1)';
    param_lin_sz = param_lin_sz((1:nvars_cc)');
    Dval = diff_tab(comp,:);
    Dval_cc = Dval(1,isvar);
    param_idx = (Dval_cc>0)*param_lin_sz + 1;
    
    if all(Dval_cc==0)
        % If the PDE state component is not differentiable along any 
        % dimension, then the PDE state component is just equal to the
        % fundamental state component.
        intrr_params{Pop_rnum,Pop_rnum}{param_idx}(rindcs,rindcs) = eye(nx_comp);
        continue
    end
    
    % The highest order derivative along each spatial direction is also
    % scaled with a particular value, which we collect in an array:
    % Establish for each maximal order of the derivative a factor with
    % which to scale the derivative of the PDE state.
    isdiff = Dval>0;
    intrr_fctr_arr = polynomial(ones(1,nvars));
    var1_d = vars(isdiff,1);    % Variables wrt which the considered state component is differentiable
    var2_d = vars(isdiff,2);
    Dval_d = Dval(isdiff);
    intrr_fctr_arr(isdiff) = (1./factorial(Dval_d-1)).*(var1_d'-var2_d').^(Dval_d-1);
    
    % For each spatial direction along which the state is differentiable,
    % we express the PDE state in terms of the highest order derivative of
    % the state along this direction, and all the lower order derivatives
    % of the state at the boundary of the domain. All of these "core
    % boundary" states need to be scaled with a factor that depends on the
    % order of the derivative of the state along each spatial direction.
    % For each direction, we collect these factors for each allowed
    % derivative along this direction in the "bndry_fctr_cell".
    bndry_fctr_cell = cell([isdiff+1,1]);
    bndry_fctr_cell{1} = 1;
    deglist_cell = cell(isdiff+1);
    deglist_cell{1} = zeros(1,0);
    
    fctr_cell = bndry_fctr_cell;
    fctr_cell{1} = prod(intrr_fctr_arr);
    
    % Set the factor associated to the fundamental state.
    intrr_params{Pop_rnum,Pop_rnum}{param_idx} = polynomial(intrr_params{Pop_rnum,Pop_rnum}{param_idx});
    intrr_params{Pop_rnum,Pop_rnum}{param_idx}(rindcs,rindcs) = kron(fctr_cell{1},eye(nx_comp));
    
    
    % We loop over each possible combination of variables wrt which the
    % state is differentiable.
    for kk=2:numel(bndry_fctr_cell)
        % Establish which combination of variables is associated to "kk".
        sub_indx_kk = cell(1,nvars);
        [sub_indx_kk{:}] = ind2sub([isdiff+1,1],kk);                        % DJ, 12/06/2024
        sub_indx_kk = cell2mat(sub_indx_kk);
        log_indx_kk = logical(sub_indx_kk-1);
        
        % Establish which of the variables in the combination is last.
        vnum_kk = find(log_indx_kk,1,'last');
        
        % For the other variables in the combination, we already have a
        % factor associated to their combination.
        log_indx_old = log_indx_kk;
        log_indx_old(vnum_kk) = false;
        if ~any(log_indx_old)
            % % If the new variable is the only one in the combination, we
            % % establish a set of factors associated to all allowed
            % % derivatives of the PDE state wrt this variable.
            
            % Determine the order of differentiability along the direction
            % assocciated to vnum_kk.
            Dval_kk = Dval(vnum_kk);
            % Extract the spatial variable associated to vnum_kk.
            var1_kk = vars(vnum_kk,1);
            % Extract the lower boundary of the domain associated to
            % vnum_kk.
            aval_kk = dom(vnum_kk,1);
            % Establish all allowed derivative orders of the state at the
            % boundary var1_kk=aval_kk.
            deglist_leq_D = (0:Dval_kk-1);
            
            bndry_fctr_cell{kk} = (1./factorial(deglist_leq_D)).*(var1_kk-aval_kk).^(deglist_leq_D);  % Factors mapping core boundary state to PDE state
            deglist_cell{kk} = zeros(Dval_kk,nvars);
            deglist_cell{kk}(:,vnum_kk) = deglist_leq_D;

        else
            % % If the new variable is considered in combination with some
            % % other variables, multiply the factor associated to the new
            % variable, with that associated to the earlier combination.
            
            % Establish the factor associated to the earlier combination.
            lin_indx_old = log_indx_old*cumprod([1;(isdiff(1:end-1)+1)'])+1;
            fctr_kk_old = bndry_fctr_cell{lin_indx_old};
            deglist_old = deglist_cell{lin_indx_old};
            
            % Establish the factor associated to the new variable.
            lin_indx_new = 2^(vnum_kk-1) + 1;
            fctr_kk_new = bndry_fctr_cell{lin_indx_new};
            deglist_new = deglist_cell{lin_indx_new};
        
            % Establish a factor for each allowed combination of derivatives
            % along the considered spatial dimensions.
            deglist_cell{kk} = repmat(deglist_old,[size(deglist_new,1),1]) + kron(deglist_new,ones(size(deglist_old,1),1));
            bndry_fctr_cell{kk} = kron(fctr_kk_new,fctr_kk_old);
        end
        
        % For any variable not included in the combination, we take a
        % derivative of maximal order, and scale accordingly with the
        % interior factor.
        is_intrr = ~log_indx_kk;    % Along which variables are we considering the interior of the domain?
        is_intrr(~isvar) = false;   % Variables that the component does not depend on should not contribute
        fctr_cell_kk = bndry_fctr_cell{kk}*prod(intrr_fctr_arr(is_intrr));
        
        
        % % % Set the factors in the appropriate parameter.
        
        % Establish which parameter in the BC operator we are considering.
        %is_intrr = ~logical(sub_indx_kk-1); 
        Pop_cnum = is_intrr*(2.^(0:nvars-1))'+1;  % Which parameter in the BC operator are we considering?
        
        % Establish which element of the parameter we are considering: we
        % should always perform partial integration along directions in
        % which the core boundary state varies.
        if nvars<=1
            lin_idx_param = 1;
        else
            log_idx_param = (is_intrr & isdiff);
            sz_param = size(bndry_params{Pop_rnum,Pop_cnum});
            lin_idx_param = log_idx_param*[1;cumprod(sz_param(1:end-1))] + 1;
        end
        %lin_idx_param = ceil(numel(bndry_params{Pop_rnum,Pop_cnum})/2); % NOT QUITE... if Dval=[0,1], should be in R2x{1} not R2x{2}
        
        % Establish which columns of the parameter our factors should be
        % placed in.
        nBCs_kk = nx_comp*length(fctr_cell_kk);    % How many core boundary states do we have?
        %lin_indx_nBCs = log_indx_kk*(2.^(0:nvars-1))'+1;
        cindcs = nBCs_arr(Pop_cnum)+1 : nBCs_arr(Pop_cnum)+nBCs_kk;
        nBCs_arr(Pop_cnum) = nBCs_arr(Pop_cnum)+nBCs_kk;
        
        % Set the boundary parameter factors.    
        bndry_params{Pop_rnum,Pop_cnum}{lin_idx_param} = polynomial(bndry_params{Pop_rnum,Pop_cnum}{lin_idx_param});
        bndry_params{Pop_rnum,Pop_cnum}{lin_idx_param}(rindcs,cindcs) = kron(fctr_cell_kk,eye(nx_comp));
        
        % Keep track of which variables the considered core boundary
        % components depends on, and to what degree it is differentiated
        % wrt each of these variables.
        x_cindcs = cindcs + sum(nbc_op(1:Pop_cnum-1));
        bc_state_tab(x_cindcs,1) = x_tab(comp,1);
        bc_state_tab(x_cindcs,2) = nx_comp;
        bc_state_tab(x_cindcs,3:2+nvars) = repmat(is_intrr,[length(x_cindcs),1]);
        deglist_kk = deglist_cell{kk};
        deglist_kk(:,is_intrr) = diff_tab(comp,is_intrr);
        deglist_kk = kron(deglist_kk,ones(nx_comp,1));
        bc_state_tab(x_cindcs,3+nvars:2+2*nvars) = deglist_kk;
        
    end
    
end

% Having determined values for the parameters, now assign the parameters 
% to the pre-initialized operators.
if nvars<=1
    Pop_f2x.P = intrr_params{1,1}{1};
    Pop_f2x.R.R0 = intrr_params{2,2}{1};
    Pop_f2x.R.R1 = intrr_params{2,2}{2};
    Pop_f2x.R.R2 = intrr_params{2,2}{3};

    Pop_b2x.Q2 = bndry_params{2,1}{1};
else
    Pop_f2x.R00 = intrr_params{1,1}{1};
    Pop_f2x.Rxx = intrr_params{2,2};
    Pop_f2x.Ryy = intrr_params{3,3};
    Pop_f2x.R22 = intrr_params{4,4};
    
    Pop_b2x.Rx0 = bndry_params{2,1}{1};
    Pop_b2x.Ry0 = bndry_params{3,1}{1};
    Pop_b2x.R20 = bndry_params{4,1}{1};
    Pop_b2x.R2x = bndry_params{4,2};
    Pop_b2x.R2y = bndry_params{4,3};
end

bc_state_tab = unique(bc_state_tab,'stable','rows');

end