function PDE = combine_terms(PDE,objs)
% This function combines terms in the equations of "objs" in the PDE
% structure "PDE"
%
% INPUTS:
% - PDE:    A pde_struct type object, defining a PDE in the terms
%           format (see also the "@pde_struct/initialize" function).
% - obj:    'x', 'y', 'z', or 'BC', indicating for what object to
%           combine terms in the associated equations Defaults to 'all', 
%           in which case terms in all equations are combined where 
%           possible.
%
%
% OUTPUTS:
% - PDE:    A structure of the same type as the input, with terms in the
%           specified equations merged where possible. That is, if two
%           terms involve the same state component or input, same delay,
%           same spatial derivative, same spatial position, and same
%           integral, these terms will have been combined in the new PDE
%           structure, with the coefficients in the two terms added to form
%           the coefficients of the new term.
%
% NOTES:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 11/28/2022
%


% % % Check the input arguments
if nargin==1
    % If no particular object is specified, reorder all components
    objs = 'all'; 
end

% If multiple objects are specified, loop over each object in the function.
if isa(objs,'char') && strcmp(objs,'all')
    objs = {'x','z','y','BC'};
elseif isa(objs,'char') && ~ismember(obj,{'x','z','y','BC'})
    error('The second argument must be one of {''x'',''z'',''y'',''BC'',''all''}.')
elseif isa(objs,'char')
    objs = {objs};
end

% Extract number of variables.
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
nvars = PDE.dim;

% % % Loop over equations for each of the specified objects, combining
% % % terms in the equations wherever possible.
for kk=1:numel(objs)

obj = objs{kk};

for ii=1:numel(PDE.(obj))
    % % Use integer array "term_tab" to keep track of
    % - the state component or input in each term (column 1);
    % - the time delay in each term (column 2);
    % - the order of the derivative wrt each spatial variable 
    %   (columns 3 through 2+nvars);
    % - the spatial position at which to evaluate each variable
    %   (columns 3+nvars through 2+2*nvars);
    % - the integral taken of the term
    %   (columns 3+2*nvars through 2+3*nvars).
    term_tab = zeros(numel(PDE.(obj){ii}.term),1+1+3*nvars);
    
    % Loop over the terms, checking what combination of component, delay,
    % derivative, position, and integral it involves.
    for jj=1:numel(PDE.(obj){ii}.term)
        term_jj = PDE.(obj){ii}.term{jj};
        % First check what state or input component the term involves.
        if isfield(term_jj,'x')
            is_x_Robj = true;
            term_tab(jj,1) = term_jj.x;
            has_vars_Rcomp = logical(PDE.x_tab(term_jj.x,3:2+nvars));
        elseif isfield(term_jj,'w')
            is_x_Robj = false;
            term_tab(jj,1) = numel(PDE.x)+term_jj.w;
            has_vars_Rcomp = logical(PDE.w_tab(term_jj.w,3:2+nvars));
        else
            is_x_Robj = false;
            term_tab(jj,1) = numel(PDE.x)+numel(PDE.w)+term_jj.u;
            has_vars_Rcomp = logical(PDE.u_tab(term_jj.u,3:2+nvars));
        end
        dom_Rcomp = PDE.dom(has_vars_Rcomp,:);
        nvars_Rcomp = sum(has_vars_Rcomp);
        
        % Next, determine what (constant) delay the term involves.
        if ~isa(term_jj.delay,'double') && ~isdouble(term_jj.delay)
            error('Currently only constant delays are supported.')
        end
        term_tab(jj,2) = double(term_jj.delay);
        
        % Determine what derivative is taken of the state,
        % and at which position the state is evaluated.
        Dval = zeros(1,nvars);      % degrees of derivatives in this term
        locval = zeros(1,nvars);    % -1 for lower B, 0 for interior, 1 for upper B
        if is_x_Robj
            Dval(has_vars_Rcomp) = double(term_jj.D);
            locval_R = zeros(1,nvars_Rcomp);
            for kk=1:nvars_Rcomp
                if isa(term_jj.loc(kk),'double') || isdouble(term_jj.loc(kk))
                    if double(term_jj.loc(kk))==dom_Rcomp(kk,1)
                        % State is evaluated at lower boundary.
                        locval_R(kk) = -1;
                    else
                        % State is evaluated at upper boundary.
                        locval_R(kk) = 1;
                    end
                end
            end
            locval(has_vars_Rcomp) = locval_R;            
        end
        
        % Determine what integral is taken of the term.
        Ival = zeros(1,nvars);
        Ival_R = zeros(1,nvars_Rcomp);
        for kk=1:nvars_Rcomp
            if ~isempty(term_jj.I{kk})
                if isa(term_jj.I{kk},'double') || isdouble(term_jj.I{kk})
                    % Integrate over full domain of variable kk.
                    Ival_R(kk) = 3;
                elseif isa(term_jj.I{kk}(1),'double') || isdouble(term_jj.I{kk}(1))
                    % Integrate over lower "half" of domain of variable kk.
                    Ival_R(kk) = 1;
                else
                    % Integrate over upper "half" of domain of variable kk.
                    Ival_R(kk) = 2;
                end
            end
        end
        Ival(has_vars_Rcomp) = Ival_R;
        
        % Keep track of the derivative, position, and integral associated
        % to this term jj.
        term_tab(jj,3:2+3*nvars) = [Dval,locval,Ival];
    end
    
    % Check which terms involve the same component, delay, derivative,
    % position, and integral.
    [~,retain_terms,new_indcs] = unique(term_tab,'rows');
    % Keep only one term for each combination of component, ..., integral.
    terms_new = PDE.(obj){ii}.term(retain_terms');
    for jj=1:numel(PDE.(obj){ii}.term)
        % Add coefficients of different terms associated to same
        % combination of component, ..., integral.
        if ismember(jj,retain_terms)
            continue
        else
            terms_new{new_indcs(jj)}.C = terms_new{new_indcs(jj)}.C + PDE.(obj){ii}.term{jj}.C;
        end
    end
end

end

% Re-initialize the PDE structure.
PDE = initialize(PDE,true);

end