function [PDE] = expand_PIETOOLS_PDE_BCs(PDE)
% Construct PI operators Pop_f2x and Pop_b2x such that the PDE state x
% can be expressed in terms of the associated fundamental state xf and a
% set of core boundary values xb as
%   x = Pop_f2x * xf + Pop_b2x * xb;
% 
% INPUTS:
% - PDE:    A struct or pde_struct object defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - expand_PIETOOLS_PDE_BCs
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
%


% Initialize the PDE, if it seems this has not yet been done.
if ~isa(PDE,'pde_struct') || (size(PDE.BC_tab,1)~=sum(sum(PDE.x_tab(:,3+PDE.dim:2+2*PDE.dim))))
    PDE = initialize_PIETOOLS_PDE(PDE,true);
end

% Extract some important parameters
vars = PDE.vars;
var1 = PDE.vars(:,1);
dom = PDE.dom;
nvars = length(var1);
global_varname_1 = cell(nvars,1);
for kk=1:nvars
    global_varname_1{kk} = var1(kk).varname{1};
end
pvar dumvar

% Set a new PDE structure, in which the BCs are expanded and sorted
PDE_new = PDE;
PDE_new.BC = cell(0,1);
BC_tab = zeros(0,2+2*nvars);

for ii=1:numel(PDE.BC)
    % Loop over all the boundary conditions
    BC_ii = PDE.BC{ii};
    
    % For each boundary function that varies along a spatial dimension, we
    % decompose the function into edge and corner values.
    has_vars_BC_ii = PDE.BC_tab(ii,3:2+nvars);       % Which variables does the BC depend on?
    nBCs_ii = prod(2.^has_vars_BC_ii);               % How many different types of boundaries does the boundary comprise?
    % Decompose e.g. plane into plane+edges+corners
    BC_new_ii = cell(nBCs_ii,1);
    BC_new_ii{1}.size = BC_ii.size;
    BC_new_ii{1}.vars = BC_ii.vars;
    BC_new_ii{1}.dom = BC_ii.dom;
    BC_new_ii{1}.term = cell(1,0);
    BC_new_ii{1}.int_tab = zeros(0,nvars);
    BC_new_ii{1}.is_xcomp = zeros(0,nvars);
    % For each new BC 0=F_kk(s), keep track of which variables the function
    % F_kk depends on using binary indices.
    dep_tab_expanded_ii = zeros(nBCs_ii,nvars);
    % Loop over all possible combinations of the nvars variables.
    for kk=2:nBCs_ii
        % Establish which combination of variables is associated to index
        % kk.
        sub_indx_kk = cell(1,nvars);
        [sub_indx_kk{:}] = ind2sub(has_vars_BC_ii+1,kk);
        sub_indx_kk = cell2mat(sub_indx_kk);
        dep_tab_expanded_ii(kk,:) = sub_indx_kk - 1;
        
        % Initialize a new (empty) boundary condition 0=F_kk(s) that
        % depends on the considered combination of variables.
        BC_new_ii{kk}.size = BC_new_ii{1}.size;
        BC_new_ii{kk}.vars = BC_new_ii{1}.vars;
        BC_new_ii{kk}.dom = BC_new_ii{1}.dom;
        BC_new_ii{kk}.term = cell(1,0);
        
        %BC_new_ii{kk}.isbc_tab = BC_new_ii{1}.isbc_tab;
        BC_new_ii{kk}.int_tab = BC_new_ii{1}.int_tab;
        BC_new_ii{kk}.is_xcomp = BC_new_ii{1}.is_xcomp;
    end
    % Build a new table with BC information associated to the new
    % (expanded) boundary conditions.
    BC_tab_ii = zeros(nBCs_ii,2+2*nvars);
    BC_tab_ii(:,2) = PDE.BC_tab(ii,2);
    BC_tab_ii(:,3:2+nvars) = dep_tab_expanded_ii;
    
    % For some BCs, decomposition into e.g. separate edge and corner BCs is
    % not possible. For this, we keep track of which variables MUST appear
    % in the current BC ii.
    remove_BC_ii = false(nBCs_ii,1);
    
    % Loop over all terms in the BC ii, decomposing them into terms
    % depending on the core boundary and fundamental state components.
    jj = 1;
    while(jj<=numel(BC_ii.term))
        term_jj = BC_ii.term{jj};
        
        % % First check which state component or input is involved.
        if isfield(term_jj,'x')
            Robj = 'x';
            is_x_Robj = true;
        elseif isfield(term_jj,'w')
            Robj = 'w';
            is_x_Robj = false;
        else
            Robj = 'u';
            is_x_Robj = false;
        end
        Rindx = find(PDE.([Robj,'_tab'])(:,1)==term_jj.(Robj));
        
        % Establish which spatial variables the component depends on.
        has_vars_Rcomp = logical(PDE.([Robj,'_tab'])(Rindx(1),3:2+nvars));
        has_vars_BC_ii_Rcomp = has_vars_BC_ii(has_vars_Rcomp);
        nvars_Rcomp = sum(has_vars_Rcomp);
        dom_Rcomp = dom(has_vars_Rcomp,:);
        vars_Rcomp = vars(has_vars_Rcomp,:);
        
        % For now, do not allow terms in the BCs involving
        % infinite-dimensional inputs.
        if (strcmp(Robj,'w') && any(PDE.w_tab(Rindx,3:2+nvars))) || ...
                (strcmp(Robj,'u') && any(PDE.u_tab(Rindx,3:2+nvars)))
            error(['Term "BC{',num2str(ii),'}.term{',num2str(jj),'}" is not appropriate: ',...
                    ' boundary conditions involving infinite-dimenstional inputs are currently not supported.'])            
        elseif ~is_x_Robj
            % If the term involves an input, we add it's contribution to
            % each of the new boundary conditions, and move on to the next
            % term.            
            for kk=1:nBCs_ii
                BC_new_ii{kk}.term = [BC_new_ii{kk}.term, BC_ii.term{jj}];
                var_log_kk = logical(dep_tab_expanded_ii(kk,:));
                BC_new_ii{kk}.term{end}.C = subs(BC_new_ii{kk}.term{end}.C,vars(~var_log_kk,1),dom(~var_log_kk,1));
                BC_new_ii{kk}.is_xcomp = [BC_new_ii{kk}.is_xcomp, false];
            end
            jj = jj+1;
            continue
        end
        
        % % % If the considered term involves a state component, decompose
        % % % it into core boundary values and the fundamental state.
        
        % Extract the order of the derivative and spatial position at which
        % to evaluate the state.
        Dval = term_jj.D;
        Rloc = term_jj.loc;
        Idoms = term_jj.I;
        Cval = polynomial(term_jj.C);
        
        % Determine the maximal order of differentiability of the state
        % component wrt each spatial variable it depends on.        
        x_diff_max = PDE.x_tab(Rindx,3+nvars:2+2*nvars);
        x_diff_max = x_diff_max(has_vars_Rcomp);
        
        % Establish which variables are evaluated at a boundary.
        isbndry_loc = true(1,nvars_Rcomp);
        int_type_list = zeros(1,nvars_Rcomp);
        for ll=1:nvars_Rcomp
            isbndry_loc(ll) = isdouble(Rloc(ll));
            if isempty(Idoms{ll})
                continue
            elseif isa(Idoms{ll},'double') || isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll})
                int_type_list(ll) = 3;  % Integrating over full spatial domain
                if ~has_vars_BC_ii_Rcomp(ll)
                    % For full integrals, kernels may have been specified
                    % using primary variables, rather than dummy vars. For
                    % use in this function, make sure kernels are specified
                    % using dummy vars.
                    Cval = subs(Cval,vars_Rcomp(ll,1),vars_Rcomp(ll,2));
                end
            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(1))
                int_type_list(ll) = 1;  % Integrating over lower "half" of domain
            elseif isa(Idoms{ll},'polynomial') && isdouble(Idoms{ll}(2))
                int_type_list(ll) = 2;  % Integrating over upper "half" of domain
            else
                error(['The proposed integral in term ',num2str(jj),' seems to be over an inappropriate domain.',...
                        ' There must be a bug in the initialization file...'])
            end
        end
        
        
        % % We check now for two cases in which the state needs to be
        % % decomposed further:
        % % - x(s) is evaluated at an upper boundary s=b for s\in[a,b]
        % % - x(s) is evaluated at s=s or s=a, but x is not differentiated
        % %     to the maximal degree as allowed per the continuity
        % %     constraints
        % % First check for the first case, determining whether any
        % % proposed location matches the upper boundary of the interval in
        % % which the corresponding variable exists.
        bndry_locs = double(Rloc(isbndry_loc));
        dom_bndry_vars = dom_Rcomp(isbndry_loc,:);
        is_bval = (bndry_locs==dom_bndry_vars(:,2));
        
        if any(is_bval)
            % % The state is evaluated at some upper boundary.
            % % In this case, only focus on the first variable for which
            % % the state is evaluated on the upper boundary. We deal with
            % % all other variables in a later iteration of the while loop.
            is_bval_indx = find(is_bval,1,'first');
            is_bndry_indcs = find(isbndry_loc,is_bval_indx,'first');
            is_bval_indx = is_bndry_indcs(end); % index of first variable of the state which is evaluated at the upper boundary?
            
            % Extract the spatial variable which is considered at an upper
            % boundary, and the associated domain.
            var11 = vars_Rcomp(is_bval_indx,1);
            var22 = vars_Rcomp(is_bval_indx,2);
            aval = dom_Rcomp(is_bval_indx,1);
            bval = dom_Rcomp(is_bval_indx,2);
            
            % Make sure the integral and order of the derivative along the
            % considered direction make sense.
            if ~isempty(Idoms{is_bval_indx})
                error(['Term "BC{',num2str(ii),'}.term{',num2str(jj),'}" is not appropriate: ',...
                        ' integration along directions perpendicular to the boundary is not supported.'])   
            end
            Dval_mismatch = x_diff_max(is_bval_indx) - Dval(is_bval_indx);
            if Dval_mismatch<=0
                error(['Term "BC{',num2str(ii),'}.term{',num2str(jj),'}" is not appropriate: ',...
                        ' to evaluate a derivative d^n/ds^n x(s) at a position s=a or s=b for s in [a,b], the state must be differentiable in s up to at least order n+1.'])   
            end            

            % % Next, split the state using the FTC:
            % % x(b) = x(a) + (b-a)*d/ds x(a) + (1/2)*(b-a)^2*d^2/ds^2 x(a) +
            % %           ... + int_{a}^{b} 1/(n-1)! * (b-s)^(n-1) * d^n/ds^n x(s) ds
            % % where x is differentiable up to order at most n in s.
            
            % We squeeze the additional terms into BC_ii:
            BC_ii.term = [BC_ii.term(1:jj),repmat(BC_ii.term(jj),[1,Dval_mismatch]),BC_ii.term(jj+1:end)];
            
            % Adjust the first terms to represent the boundary terms 
            % x(a) + (b-a)*d/ds x(a) + (1/2)*(b-a)^2*d^2/ds^2 x(a) + ...
            for kk=0:Dval_mismatch-1
                BC_ii.term{jj+kk}.D(is_bval_indx) = Dval(is_bval_indx) + kk;
                BC_ii.term{jj+kk}.loc(is_bval_indx) = aval;   % Evaluate at lower boundary
                BC_ii.term{jj+kk}.C = Cval * (1./factorial(kk)).*(bval-aval).^(kk);
            end

            % Adjust the last term to represent the interior term
            % ... + int_{a}^{b} 1/(n-1)! * (b-s)^(n-1) * d^n/ds^n x(s) ds
            BC_ii.term{jj+Dval_mismatch}.D(is_bval_indx) = x_diff_max(is_bval_indx);
            BC_ii.term{jj+Dval_mismatch}.loc(is_bval_indx) = var11;   % Evaluate at interior
            BC_ii.term{jj+Dval_mismatch}.I{is_bval_indx} = [aval, bval];
            BC_ii.term{jj+Dval_mismatch}.C = Cval * (1./factorial(Dval_mismatch-1)).*(bval-var22).^(Dval_mismatch-1);            

        else
            % % The state is not evaluated at an upper boundary of the
            % % domain.
        
            % In this case, check if the order of the derivative of the
            % state in this term is as maximally allowed by the continuity
            % constraints.
            Dval_mismatch = zeros(1,nvars_Rcomp);
            Dval_mismatch(~isbndry_loc) = x_diff_max(~isbndry_loc) - Dval(~isbndry_loc);
            if ~any(Dval_mismatch)
                % The state is differentiated up to the maximally allowed
                % degree --> no further expansion is necessary.
                % Add contribution of this term to each of the decomposed
                % boundary conditions, e.g. separate corner, edge, and
                % plane BCs.                
                for kk=1:nBCs_ii
                    if remove_BC_ii(kk)
                        continue
                    end
                    var_log_kk = logical(dep_tab_expanded_ii(kk,:));
                    var_log_kk_jj = var_log_kk(has_vars_Rcomp);
                    if any(int_type_list(~var_log_kk_jj)==1)
                        % Evaluating integrals from a^s at s=a will not
                        % contribute to the BC.
                        continue
                    end
                    % Establish if the considered term depends on any
                    % variable not included in var_log_kk.
                    has_vars_term = false(1,nvars);
                    has_vars_term(has_vars_Rcomp) = ~isbndry_loc;
                    
                    % Replace integrals int_s^b with int_a^b if we're
                    % evaluating s=a.
                    make_full_int = (int_type_list==2 & ~var_log_kk(has_vars_Rcomp));
                    % Any full integral will remove the dependence on a
                    % variable
                    is_full_int = false(1,nvars);
                    is_full_int(has_vars_Rcomp) = (int_type_list==3 | make_full_int);
                    has_vars_term(is_full_int) = false;
                    
                    % For variables on which BC kk does not depend, we
                    % evaluate the variable at the lower boundary
                    Cval_kk = polynomial(subs(Cval,vars(~var_log_kk,1),dom(~var_log_kk,1)));
                    has_vars_Cval = ismember(global_varname_1,Cval_kk.varname)';
                    has_vars_term = has_vars_term | has_vars_Cval;
                    
                    if any(has_vars_term & ~var_log_kk)
                        % If the term depends on variables that do not
                        % appear in BC kk, then a proper decomposition into
                        % e.g. corner and edge BCs is not possible. Any BC
                        % involving this term must be allowed to depend on
                        % the variables appearing in this term.
                        remove_BC_ii(kk) = true;
                        continue
                    elseif max(max(abs(Cval_kk.coeff)))<=1e-13
                        % If the coefficients are zero, the term does not
                        % contribute.
                        continue
                    else    %if ~any(has_vars_term - var_log_kk)
                        % Split contributions to corners, edges, planes, etc.
                        BC_new_ii{kk}.term = [BC_new_ii{kk}.term, BC_ii.term{jj}];
                        BC_new_ii{kk}.term{end}.C = Cval_kk;
                        int_dom_kk = mat2cell(dom(make_full_int,:),ones(sum(make_full_int),1));
                        BC_new_ii{kk}.term{end}.I(make_full_int) = int_dom_kk;
                        
                        %BC_new_ii{kk}.isbc_tab(has_vars_Rcomp) = BC_new_ii{1}.isbc_tab;
                        int_type_list_full = zeros(1,nvars);
                        int_type_list_full(has_vars_Rcomp) = int_type_list;
                        BC_new_ii{kk}.int_tab = [BC_new_ii{kk}.int_tab; int_type_list_full];
                        BC_new_ii{kk}.is_xcomp = [BC_new_ii{kk}.is_xcomp, true];
                        
                        BC_diff_tab_ii_kk = BC_tab_ii(kk,3+nvars:2+2*nvars);
                        BC_diff_tab_ii_kk(has_vars_Rcomp) = max(BC_diff_tab_ii_kk(has_vars_Rcomp), Dval);
                        BC_tab_ii(kk,3+nvars:2+2*nvars) = BC_diff_tab_ii_kk;
                    end
                end                
                jj = jj+1;
                continue
            end
            
            % % % The state is not differentiated up to the maximally
            % % % degree --> we decompose it into a core boundary state and
            % % % fundamental state component:

            % Establish the first spatial variable wrt which the state must
            % still be differentiated.
            must_diff_indx = find(Dval_mismatch,1,'first');
            dom_Rcomp = dom(has_vars_Rcomp,:);
            aval = dom_Rcomp(must_diff_indx,1);
            bval = dom_Rcomp(must_diff_indx,2);
            vars_Rcomp = vars(has_vars_Rcomp,:);
            var11 = vars_Rcomp(must_diff_indx,1);
            var22 = vars_Rcomp(must_diff_indx,2);

            % Split the state using the FTC:
            % x(s) = x(a) + int_{a}^{s} d/dtheta x(theta) dtheta
            BC_ii.term = [BC_ii.term(1:jj),BC_ii.term(jj:end)];
            
            % % Adjust term jj to represent x(a).
            BC_ii.term{jj}.loc(must_diff_indx) = dom_Rcomp(must_diff_indx,1);   % Evaluate at lower boundary
            BC_ii.term{jj}.I{must_diff_indx} = [];              % Perform no integration wrt variable the state does not depend on
            if isempty(Idoms{must_diff_indx})
                BC_ii.term{jj}.C = Cval;
            else
                % If the state was integrated along the considered spatial
                % direction, this integral must be applied to our
                % coefficients.
                BC_ii.term{jj}.C = int(Cval,var22,Idoms{must_diff_indx}(1),Idoms{must_diff_indx}(2));
            end

            % % Adjust term jj+1 to represent
            % %  int_{a}^{s} d/dtheta f(theta) dtheta.
            BC_ii.term{jj+1}.D(must_diff_indx) = Dval(must_diff_indx) + 1;
            BC_ii.term{jj+1}.loc(must_diff_indx) = var22;   % Evaluate at the same position
            
            % Depending on the user-specified integral, we made need to
            % adjust the coefficients...
            if isempty(Idoms{must_diff_indx})
                % If no integral was requested, simply add "int_{a}^{s}"
                BC_ii.term{jj+1}.I{must_diff_indx} = [aval, var11];
                BC_ii.term{jj+1}.C = subs(Cval,var11,var22);    % primary variable becomes dummy variable in integral
            elseif isequal(Idoms{must_diff_indx}(1)-aval,0) && isequal(Idoms{must_diff_indx}(2)-bval,0)
                % If the state is already integrated from a to b, the
                % term is constant, and we can just multiply with (b-s).
                BC_ii.term{jj+1}.I = Idoms;
                pvar dumvar
                BC_ii.term{jj+1}.C = int(subs(Cval,var22,dumvar),dumvar,var22,bval);
            elseif isequal(Idoms{must_diff_indx}(1)-aval,0) && isequal(Idoms{must_diff_indx}(2)-var11,0)
                % If the state is already integrated from a to s, we have
                % to integrate C from theta to s:
                % int_{a}^{s} int_{a}^{theta} C(theta,eta) * s(eta) deta dtheta
                %   = int_{a}^{b} int_{a}^{b} I(s-theta) * I(theta-eta) * C(theta,eta) * s(eta) deta dtheta
                %   = int_{a}^{b} I(s-eta) int_{eta}^{s} C(theta,eta) dtheta s(eta) deta
                %   = int_{a}^{s} int_{theta}^{s} [C(eta,theta)] deta s(theta) dtheta
                BC_ii.term{jj+1}.I = Idoms;
                BC_ii.term{jj+1}.I{must_diff_indx} = [aval,var11];
                BC_ii.term{jj+1}.C = int(subs(Cval,var11,dumvar),dumvar,var22,var11);
            elseif isequal(Idoms{must_diff_indx}(1)-var11,0) && isequal(Idoms{must_diff_indx}(2)-bval,0)
                % If the state is already integrated from s to b, we have
                % to integrate split the term into two:
                % int_{a}^{s} int_{theta}^{b} C(theta,eta) * s(eta) deta dtheta
                %   = int_{a}^{b} int_{a}^{b} I(s-theta) * I(eta-theta) * C(theta,eta) * s(eta) deta dtheta
                %   = int_{a}^{b} int_{a}^{b} I(s-eta) * I(eta-theta) * C(theta,eta) * s(eta) deta dtheta
                %       + int_{a}^{b} int_{a}^{b} I(s-theta) * I(eta-s) * C(theta,eta) * s(eta) deta dtheta
                %   = int_{a}^{s} int_{a}^{theta} C(eta,theta) deta s(theta) dtheta
                %       + int_{s}^{b} int_{a}^{s} C(eta,theta) deta s(theta) dtheta
                BC_ii.term{jj+1}.I = Idoms;
                BC_ii.term = [BC_ii.term(1:jj+1), BC_ii.term(jj+1:end)];

                BC_ii.term{jj+1}.I{must_diff_indx} = [aval,var11];
                BC_ii.term{jj+1}.C = int(subs(Cval,var11,dumvar),dumvar,aval,var22);

                BC_ii.term{jj+2}.I{must_diff_indx} = [var11,bval];
                BC_ii.term{jj+2}.C = int(subs(Cval,var11,dumvar),dumvar,aval,var11);
            else
                % Otherwise, the integral is not appropriate.
                error(['The proposed integral in term ',num2str(jj),' seems to be over an inappropriate domain.',...
                        ' There must be a bug in the initialization file...'])
            end

            % Move on to the next term.
            %jj = jj+1;
        end
    end
    
    % Establish which terms vary along which boundaries
    kk = 1;
    while kk<=numel(BC_new_ii)
        % If the new BC has no terms, there's no point in keeping it.
        if isempty(BC_new_ii{kk}.term) || remove_BC_ii(kk)
            BC_new_ii = [BC_new_ii(1:kk-1,:); BC_new_ii(kk+1:end,:)];
            BC_tab_ii = [BC_tab_ii(1:kk-1,:); BC_tab_ii(kk+1:end,:)];
            remove_BC_ii = [remove_BC_ii(1:kk-1); remove_BC_ii(kk+1:end)];
            continue
        end
        % If all terms are integrated over the same indefinite domain in a
        % particular spatial direction, the integrand itself must be zero:
        % 0 = int_{a}^s f(tt) dtt --> 0 = f(s)
        is_xcomp_kk = logical(BC_new_ii{kk}.is_xcomp);
        int_tab_kk = BC_new_ii{kk}.int_tab(is_xcomp_kk,:);
        remove_int = any(int_tab_kk,1) & any(int_tab_kk<3,1) & ~any(int_tab_kk(2:end,:) - int_tab_kk(1:end-1,:),1);
        if any(remove_int)
            for jj=1:numel(BC_new_ii{kk}.term)
                if ~is_xcomp_kk(jj)
                    continue
                end
                Rindx = find(PDE.([Robj,'_tab'])(:,1)==BC_new_ii{kk}.term{jj}.x);
                % Establish which spatial variables each component depends on.
                has_vars_Rcomp = logical(PDE.x_tab(Rindx(1),3:2+nvars));
                BC_new_ii{kk}.term{jj}.I{remove_int(has_vars_Rcomp)} = [];
            end
        end
        BC_new_ii{kk} = rmfield(BC_new_ii{kk},'int_tab');
        BC_new_ii{kk} = rmfield(BC_new_ii{kk},'is_xcomp');
        kk = kk+1;
    end
    
    % Update the BCs in the PDE structure, and the BC table.
    PDE_new.BC = [PDE_new.BC; BC_new_ii];    
    BC_tab = [BC_tab; BC_tab_ii];
    
end


% Re-order the BCs based on their dependence on each of the spatial
% variables.
BC_tab(:,1) = 1:size(BC_tab,1);
dep_tab = BC_tab(:,3:2+nvars);
diff_tab = BC_tab(:,3+nvars:2*nvars+2);
BC_tab_alt = [dep_tab(:,end:-1:1), diff_tab(:,end:-1:1), BC_tab(:,1), BC_tab(:,2)];
[BC_tab_alt,new_order] = sortrows_integerTable(BC_tab_alt);
BC_tab_new = [BC_tab_alt(:,end-1:end), BC_tab_alt(:,nvars:-1:1), BC_tab_alt(:,2*nvars:-1:nvars+1)];
BC_tab_new(:,1) = 1:size(BC_tab_new,1);
BC_tab = BC_tab_new;

PDE_new.BC = PDE_new.BC(new_order);
PDE_new.BC_tab = BC_tab;
PDE = PDE_new;



end