function [PDE_new,old_ito_new] = combine_vars(PDE,dom,silent_merge)
% [PDE_new,old_ito_new] = combine_vars(PDE,dom,silent_merge)
% Takes a pde_struct and represents it using fewer spatial variables, all
% on the same domain dom.
%
% INPUTS:
% PDE:  A pde_struct class object, defining a PDE in a manner as
%       outlined in "initialize_PIETOOLS_PDE".
% dom:  (optional) A 1x2 array speicfying the new interval on which all
%       spatial variables must exist. Defaults to [-1,1];
% silent_merge: logical value indicating whether to suppress the summary of
%               the merging process. Defaults to false.
% 
% OUTPUTS:
% PDE_new:  A pde_struct object describing the same system as the input
%           PDE, but with fewer spatial variables (if possible), and with
%           all variables existing on the same interval dom. For any pair
%           of variables, the function checks whether any state component,
%           input or output depends on both variables. If so, these
%           variables cannot be "merged". If not, any instance of the
%           second variable is replaced with the first variable, discarding
%           the second variable. The variables are also scaled and shifted
%           to exist on the same domain.
% old_ito_new:  A nvarsx2 polynomial object providing an expression for the
%               variables in the original PDE structure in terms of those
%               in the new PDE structure.
%
% NOTES:
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
%
% Initial coding DJ - 10/11/2022


% % % Check the inputs
if nargin<1
    error('At least one argument must be passed.')
elseif nargin==1
    dom = [-1,1];
    silent_merge = false;
elseif nargin==2
    if islogical(dom)
        silent_merge = dom;
        dom = [-1,1];
    else
        silent_merge = false;
    end
elseif nargin>3
    error('At most three arguments are accepted.')
end

% Check that the specified domain makes sense
if any(size(dom)~=[1,2])
    error('New global domain should be specified as 1x2 array [a,b].');
elseif dom(2)<=dom(1)
    error('Upper boundary b of domain [a,b] should exceed lower boundary a.');
end


% Initialize the PDE, and extract the old variables and their domains.
if ~PDE.is_initialized
    [PDE,var_order] = initialize(PDE,true); % PDE_out.vars = PDE_in.(var_order,:);
else
    var_order = 1:size(PDE.vars,1);
end
old_vars = PDE.vars;
old_dom = PDE.dom;
nvars = size(old_vars,1);



% % % Determine which variables can be merged

% Establish on which variables each component depends
dep_tab = logical([PDE.x_tab(:,3:2+nvars);
                   PDE.w_tab(:,3:2+nvars);
                   PDE.u_tab(:,3:2+nvars);
                   PDE.y_tab(:,3:2+nvars);
                   PDE.z_tab(:,3:2+nvars);
                   PDE.BC_tab(:,3:2+nvars)]);

% Reorder variables from most to fewest dependent components
n_deps = sum(dep_tab,1);
[~,new_var_order] = sort(n_deps,'descend');
dep_tab = dep_tab(:,new_var_order);
%[~,org_var_order] = sort(new_var_order);      % keep track of old order

% For each variable 1:nvars-1, determine a variable with which it can be
% merged.
var_indcs_old = (1:nvars);
var_indcs_new = (1:nvars);
dep_tab_new = dep_tab;
for j=1:nvars-1
    % For all components depending on variable j, determine on which other
    % variables these components also depend. If there exists even one
    % component depending on both variable j and i, we cannot merge
    % variables j and i.
    dep_tab_j = dep_tab_new(dep_tab_new(:,j),j+1:end);  % which components depend on variable j
    can_merge_j = ~any(dep_tab_j,1);                    % with which variables could j be merged?
    new_var_indx_j = find(can_merge_j,1,'first')+j;        % first var i with which var j can be merged
    if ~isempty(new_var_indx_j)
        % Replace the index of var j with the new index
        var_indcs_new(var_indcs_new==new_var_order(j)) = new_var_order(new_var_indx_j);
        % Update dependency to indicate that var j = var new_var_indx_j
        dep_tab_new(:,new_var_indx_j) = dep_tab_new(:,j) | dep_tab_new(:,new_var_indx_j);
    end
end



% % % Scale all variables to exist on same interval dom

% Determine which of the old variables are retained.
[~,var_indcs_retain,var_indcs_o2n] = unique(var_indcs_new,'stable');
nvars_new = length(var_indcs_retain);

% Build a new set of variables (s_i,theta_i).
new_varnames = cell(nvars_new,2);
for jj=1:nvars_new
    new_varnames{jj,1} = ['s',num2str(jj)];
    new_varnames{jj,2} = ['th',num2str(jj)];
end
new_vars = polynomial(new_varnames);
new_dom = repmat(dom,[nvars_new,1]);

% Assign to each old variable and domain a new variable and domain.
new_vars_aug = new_vars(var_indcs_o2n,:);
new_dom_aug = new_dom(var_indcs_o2n,:);

% Express old vars in terms of new vars.
ai = (old_dom(:,2)-old_dom(:,1))./(new_dom_aug(:,2)-new_dom_aug(:,1));
bi = old_dom(:,1) - ai.*new_dom_aug(:,1);
old_ito_new = repmat(ai,[1,2]).*new_vars_aug + repmat(bi,[1,2]);


% % % Build the new PDE structure

% Initialize the structure as just the old PDE.
PDE_new = PDE;

% Update the global variables and their domains
PDE_new.vars = new_vars;
PDE_new.dom = new_dom;

% % % Update the object tables x_tab, etc. to match the new dependence.
PDE_new.x_tab = PDE.x_tab(:,[1, 2, 2+var_indcs_retain', 2+nvars_new+var_indcs_retain']);
PDE_new.w_tab = PDE.w_tab(:,[1, 2, 2+var_indcs_retain']);
PDE_new.z_tab = PDE.z_tab(:,[1, 2, 2+var_indcs_retain']);
PDE_new.u_tab = PDE.u_tab(:,[1, 2, 2+var_indcs_retain']);
PDE_new.y_tab = PDE.y_tab(:,[1, 2, 2+var_indcs_retain']);
PDE_new.BC_tab = PDE.BC_tab(:,[1, 2, 2+var_indcs_retain', 2+nvars_new+var_indcs_retain']);

% % % Replace the variables in each term.
PDE_new = replace_vars_terms(PDE,PDE_new,'x',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);
PDE_new = replace_vars_terms(PDE,PDE_new,'w',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);
PDE_new = replace_vars_terms(PDE,PDE_new,'u',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);
PDE_new = replace_vars_terms(PDE,PDE_new,'z',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);
PDE_new = replace_vars_terms(PDE,PDE_new,'y',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);
PDE_new = replace_vars_terms(PDE,PDE_new,'BC',old_vars,new_vars_aug,old_dom,new_dom_aug,old_ito_new,ai);

% % % Re-initialize the PDE
PDE_new = initialize(PDE_new,true);


% Indicate in the command line window which variables have been changed.
ismerge = var_indcs_old~=var_indcs_new;
if ~silent_merge
    if ~any(ismerge)
        fprintf('\n','No variables have been merged.\n')
    else
        old_vars_merge = old_vars(var_indcs_old(ismerge),1);
        new_vars_merge = old_vars(var_indcs_new(ismerge),1);
        old_var_names = old_vars_merge(1).varname{1};
        new_var_names = new_vars_merge(1).varname{1};
        if sum(ismerge)==1
            fprintf(['\n','Variable ',old_var_names,' has been merged with variable ',new_var_names,'.\n'])  
        else
            for kk=2:size(old_vars_merge,1)
                old_var_names = [old_var_names,',',old_vars_merge(kk).varname{1}];
                new_var_names = [new_var_names,',',new_vars_merge(kk).varname{1}];
            end
            fprintf(['\n','Variables (',old_var_names,') have been merged with variables (',new_var_names,') respectively.\n'])  
        end
    end
    fprintf(['\n','All spatial variables have been rescaled to exist on the interval [',num2str(dom(1)),',',num2str(dom(2)),'].\n'])
end

if nargout>1
    % By intializing the PDE, the order of the variables may have changed.
    % If the user requests an expression for the old variables in terms of
    % the new ones, we'll have to reorder our expression "old_ito_new" to
    % match the order of the variables in the input PDE.
    [~,old_var_order] = sort(var_order);
    old_ito_new = old_ito_new(old_var_order,:);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function PDE_new = replace_vars_terms(PDE,PDE_new,Lobj,old_vars,new_vars,old_dom,new_dom,new2old,fctr_new2old)
% Loop over all components of the object Lobj, replacing the old variables
% in their PDEs with the new variables.
nvars = size(old_vars,1);
ncomps = numel(PDE.(Lobj));
for ii=1:ncomps    
    % Establish which old variables s the component depends on
    has_vars_Lcomp = logical(PDE.([Lobj,'_tab'])(ii,3:2+nvars));
    % Update to the new variables t and their associated domain 
    old_vars_Lcomp = old_vars(has_vars_Lcomp,:);    % old vars [s,theta]
    new_vars_Lcomp = new_vars(has_vars_Lcomp,:);
    new2old_Lcomp = new2old(has_vars_Lcomp,:);      % map s(t) = var_fctr*t + b;
    PDE_new.(Lobj){ii}.vars = new_vars_Lcomp;
    PDE_new.(Lobj){ii}.dom = new_dom(has_vars_Lcomp,:);
    
    % Loop over all terms in the PDE for the component, replacing the old
    % variables by the new variables.
    if strcmp(Lobj,'w') || strcmp(Lobj,'u')
        break
    end
    nterms = numel(PDE.(Lobj){ii}.term);
    for jj=1:nterms
        % Extract the term;
        term_jj = PDE.(Lobj){ii}.term{jj};
        
        % Establish whether the term involves an input or state.
        if isfield(term_jj,'x')
            Robj = 'x';
            is_x_Rcomp = true;            
        elseif isfield(term_jj,'w')
            Robj = 'w';
            is_x_Rcomp = false;
        else
            Robj = 'u';
            is_x_Rcomp = false;
        end
        Rindx = term_jj.(Robj);
        % Determine on which global variables the RHS component depends.
        has_vars_Rcomp = logical(PDE.([Robj,'_tab'])(Rindx,3:2+nvars));
        nvars_Rcomp = sum(has_vars_Rcomp);
        
        % Determine the old and new vars associated to the RHS component.
        old_vars_Rcomp = old_vars(has_vars_Rcomp,:);    % old vars [s,theta]
        new_vars_Rcomp = new_vars(has_vars_Rcomp,:);    % new vars [t,nu];
        old_dom_Rcomp = old_dom(has_vars_Rcomp,:);      % old domains
        new_dom_Rcomp = new_dom(has_vars_Rcomp,:);      % new domains
        new2old_Rcomp = new2old(has_vars_Rcomp,:);      % map s(t) = var_fctr*t + b;
        var_fctr_Rcomp = fctr_new2old(has_vars_Rcomp);  % factor var_fctr;
        
        
        % % Update the derivative: dR/ds = dR/dt * dt/ds = (1/var_fctr)*dR/dt
        if is_x_Rcomp
            Dval = term_jj.D;
            fctr = prod((1./var_fctr_Rcomp').^Dval);
        else
            fctr = 1;
        end
        
        % % Update the spatial position: x(s) --> x(t), x(L(s)) --> x(L(t))
        if is_x_Rcomp
            locval = term_jj.loc;
            for kk=1:nvars_Rcomp
                if isa(locval(kk),'double') || isdouble(locval(kk))
                    locval_kk = double(locval(kk));
                    if locval_kk == old_dom_Rcomp(kk,1)
                        locval(kk) = new_dom_Rcomp(kk,1);
                    elseif locval_kk == old_dom_Rcomp(kk,2)
                        locval(kk) = new_dom_Rcomp(kk,2);
                    else
                        error(['The position "',Lobj,'{',num2str(ii),'}.term{',num2str(jj),'}.loc" is not properly specified...']);
                    end
                elseif strcmp(locval(kk).varname{1},old_vars_Rcomp(kk,1).varname{1})
                    locval(kk) = new_vars_Rcomp(kk,1);
                elseif strcmp(locval(kk).varname{1},old_vars_Rcomp(kk,2).varname{1})
                    locval(kk) = new_vars_Rcomp(kk,2);
                else
                    error(['The position "',Lobj,'{',num2str(ii),'}.term{',num2str(jj),'}.loc" is not properly specified...']);
                end                
            end
            term_jj.loc = locval;
        end
        
        % % Update the integral: 
        % % int_{L(s)}^{U(s)} ... ds = int_{L(t)}^{U(t)} ...*(ds/dt) dt
        is_partial_int = false(1,nvars_Rcomp);
        for kk=1:nvars_Rcomp
            Idom_kk = term_jj.I{kk};
            if ~isempty(Idom_kk)
                is_partial_int(kk) = any(ismember(new_vars_Rcomp(kk,1).varname{1},new_vars_Lcomp(:,1).varname));
                %is_partial_int(kk) = is_partial_int(kk) && replace_var_Rcomp(kk);
                if isa(Idom_kk,'double') || isdouble(Idom_kk)
                    % Full integral: just have to replace with new domain.
                    Idom_kk = new_dom_Rcomp(kk,:);
                else
                    % Partial integral: will have to change the variable.
                    if isdouble(Idom_kk(1))
                        Idom_kk(1) = new_dom_Rcomp(kk,1);
                    elseif strcmp(Idom_kk(1).varname{1},old_vars_Rcomp(kk,1).varname{1})
                        Idom_kk(1) = new_vars_Rcomp(kk,1);
                    elseif strcmp(Idom_kk(1).varname{1},old_vars_Rcomp(kk,2).varname{1})
                        Idom_kk(1) = new_vars_Rcomp(kk,2);
                    else
                        error(['The integral "',Lobj,'{',num2str(ii),'}.term{',num2str(jj),'}.I" is not properly specified...']);
                    end
                    if isdouble(Idom_kk(2))
                        Idom_kk(2) = new_dom_Rcomp(kk,2);
                    elseif strcmp(Idom_kk(2).varname{1},old_vars_Rcomp(kk,1).varname{1})
                        Idom_kk(2) = new_vars_Rcomp(kk,1);
                    elseif strcmp(Idom_kk(2).varname{1},old_vars_Rcomp(kk,2).varname{1})
                        Idom_kk(2) = new_vars_Rcomp(kk,2);
                    else
                        error(['The integral "',Lobj,'{',num2str(ii),'}.term{',num2str(jj),'}.I" is not properly specified...']);
                    end
                end                        
                term_jj.I{kk} = Idom_kk;  
                fctr = fctr * var_fctr_Rcomp(kk);
            end
        end
            
        % % Finally, update the coefficients: C(s) --> C(t)*fctr        
        Cval = term_jj.C;
        if isa(Cval,'polynomial')
            % Substitute in expression for s in terms of t
            if size(old_vars_Lcomp,1)>0
                Cval = subs(Cval,old_vars_Lcomp(:,1),new2old_Lcomp(:,1));
            end
            if any(is_partial_int)
                new2old_Rcomp = subs(new2old_Rcomp,new_vars_Rcomp(is_partial_int,1),new_vars_Rcomp(is_partial_int,2));
            end
            if size(old_vars_Rcomp,1)>0
                Cval = subs(Cval,old_vars_Rcomp(:),new2old_Rcomp(:));
            end
        end
        % Multiply with scaling due to integration/differentiation.
        term_jj.C = Cval*fctr;
        
        % % Update the term in the new PDE, and continue to the next term.
        PDE_new.(Lobj){ii}.term{jj} = term_jj;
    end
end

end