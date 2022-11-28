function [logval] = isequal_PDEs_terms(PDE1,PDE2,no_reorder)
%ISEQUAL_PDES_TERMS Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    error('At least two PDE structures must be specified to compare.')
elseif nargin==2
    no_reorder = false;
end

if ~isa(PDE1,'pde_struct')
    try PDE1 = pde_struct(PDE1);
    catch
        error('PDEs to compare must be specified as ''pde_struct'' class objects.')
    end
end
if ~isa(PDE2,'pde_struct')
    try PDE2 = pde_struct(PDE2);
    catch
        error('PDEs to compare must be specified as ''pde_struct'' class objects.')
    end
end

% Express PDEs in a standardized format for comparison.
if ~no_reorder
    % % First PDE:
    % Initialize.
    if ~PDE1.is_initialized
        PDE1 = initialize(PDE1,true);
    end
    % Reorder components.
    PDE1 = reorder_comps(PDE1,'all',true);
    % Scale to standard domain [0,1].
    dom = [zeros(PDE1.dim,1),ones(PDE1.dim,1)];
    PDE1 = combine_vars(PDE1,dom,true);
    % Expand higher-order temporal derivatives.
    if PDE1.has_hotd
        PDE1 = expand_tderivatives(PDE1,true);
    end
    % Expand delays.
    if PDE1.has_delay
        PDE1 = expand_delays(PDE1,true);
    end
    % Reorder components again.
    PDE1 = reorder_comps(PDE1,'all',true);
    
    
    % % Second PDE:
    % Initialize.
    if ~PDE2.is_initialized
        PDE2 = initialize(PDE2,true);
    end
    % Reorder components.
    PDE2 = reorder_comps(PDE2,'all',true);
    % Scale to standard domain [0,1].
    dom = [zeros(PDE2.dim,1),ones(PDE2.dim,1)];
    PDE2 = combine_vars(PDE2,dom,true);
    % Expand higher-order temporal derivatives.
    if PDE2.has_hotd
        PDE2 = expand_tderivatives(PDE2,true);
    end
    % Expand delays.
    if PDE2.has_delay
        PDE2 = expand_delays(PDE2,true);
    end
    % Reorder components again.
    PDE2 = reorder_comps(PDE2,'all',true);
end
    
% Assume innocence until proven guilty.
logval = true;

% First check if the dimension, variables, and domain match
if PDE1.dim~=PDE2.dim
    disp('Spatial dimensionality of the presented PDEs does not match.')
    logval = false;
    return
end
if ~all(all(isequal(PDE1.vars,PDE2.vars)))
    disp('Spatial variables of the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.dom==PDE2.dom))
    disp('Spatial domains of the presented PDEs do not match.')
    logval = false;
    return
end

% % Check if the number of state components, inputs, outputs, and BCs 
% % matches
if numel(PDE1.x)~=numel(PDE2.x)
    disp('Number of state components in the presented PDEs does not match.')
    logval = false;
    return
end
if numel(PDE1.w)~=numel(PDE2.w)
    disp('Number of exogenous inputs in the presented PDEs does not match.')
    logval = false;
    return
end
if numel(PDE1.u)~=numel(PDE2.u)
    disp('Number of controlled inputs in the presented PDEs does not match.')
    logval = false;
    return
end
if numel(PDE1.z)~=numel(PDE2.z)
    disp('Number of regulated outputs in the presented PDEs does not match.')
    logval = false;
    return
end
if numel(PDE1.y)~=numel(PDE2.y)
    disp('Number of observed outputs in the presented PDEs does not match.')
    logval = false;
    return
end
if numel(PDE1.BC)~=numel(PDE2.BC)
    disp('Number of BCs in the presented PDEs does not match.')
    logval = false;
    return
end

% % Check if the components depend on the same numbers of variables
if size(PDE1.x_tab,2)~=size(PDE2.x_tab,2)
    disp('The number of spatial variables on which the state components depend does not match.')
    logval = false;
    return
end
if size(PDE1.w_tab,2)~=size(PDE2.w_tab,2)
    disp('The number of spatial variables on which the exogenous inputs depend does not match.')
    logval = false;
    return
end
if size(PDE1.u_tab,2)~=size(PDE2.u_tab,2)
    disp('The number of spatial variables on which the controlled inputs depend does not match.')
    logval = false;
    return
end
if size(PDE1.z_tab,2)~=size(PDE2.z_tab,2)
    disp('The number of spatial variables on which the regulated outputs depend does not match.')
    logval = false;
    return
end
if size(PDE1.y_tab,2)~=size(PDE2.y_tab,2)
    disp('The number of spatial variables on which the observed outputs depend does not match.')
    logval = false;
    return
end
if size(PDE1.BC_tab,2)~=size(PDE2.BC_tab,2)
    disp('The number of spatial variables on which the BCs depend does not match.')
    logval = false;
    return
end

% % Check that the specifications (size etc.) of the components match.
if ~all(all(PDE1.x_tab==PDE2.x_tab))
    disp('Properties of state components in the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.w_tab==PDE2.w_tab))
    disp('Properties of exogenous inputs in the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.u_tab==PDE2.u_tab))
    disp('Properties of controlled inputs in the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.z_tab==PDE2.z_tab))
    disp('Properties of regulated outputs in the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.y_tab==PDE2.y_tab))
    disp('Properties of observed outputs in the presented PDEs do not match.')
    logval = false;
    return
end
if ~all(all(PDE1.BC_tab==PDE2.BC_tab))
    disp('Properties of BCs in the presented PDEs do not match.')
    logval = false;
    return
end



% % % Finally, check if the terms match
logval_x = ismatch_terms(PDE1,PDE2,'x');
logval_y = ismatch_terms(PDE1,PDE2,'y');
logval_z = ismatch_terms(PDE1,PDE2,'z');
logval_BC = ismatch_terms(PDE1,PDE2,'BC');

logval = logval && logval_x && logval_y && logval_z && logval_BC;

end




function logval = ismatch_terms(PDE1,PDE2,obj)
% Check if the terms in PDE1.(obj) and PDE2.(obj) match.

logval = true;
for ii=1:numel(PDE1.(obj))
    % First check if the number of terms in the equation matches.
    if numel(PDE1.(obj){ii}.term) ~= numel(PDE2.(obj){ii}.term)
        disp('Number of terms in equation ',obj,'{',num2str(ii),'} does not match.')
        logval = false;
        continue
    end
    % Then check if the terms actually match.
    Cfctr = 1;  % For BCs, we can scale with e.g. -1 to get same conditions
    for jj=1:numel(PDE1.(obj){ii}.term)
        term1_jj = PDE1.(obj){ii}.term{jj};
        term2_jj = PDE2.(obj){ii}.term{jj};
        trm_name = [obj,'{',num2str(ii),'}.term{',num2str(jj),'}'];
        
        % Check if the considered component matches
        if isfield(term1_jj,'x')
            if ~isfield(term2_jj,'x') || term1_jj.x~=term2_jj.x
                disp('The considered components in ',trm_name,' do not match.');
                logval = false;
                continue
            end            
        elseif isfield(term1_jj,'w')
            if ~isfield(term2_jj,'w') || term1_jj.w~=term2_jj.w
                disp('The considered components in ',trm_name,' do not match.');
                logval = false;
                continue
            end
        else
            if ~isfield(term2_jj,'u') || term1_jj.u~=term2_jj.u
                disp('The considered components in ',trm_name,' do not match.');
                logval = false;
                continue
            end
        end
        
        if isfield(term1_jj,'x')
            % Check if the derivatives match
            if ~all(term1_jj.D==term2_jj.D)
                disp('The specified derivatives ',trm_name,'.D do not match.')
                logval = false;
            end
            % Check that the spatial positions match
            if ~all(isequal(term1_jj.loc,term2_jj.loc))
                disp('The specified spatial positions ',trm_name,'.loc do not match.')
                logval = false;
            end
        end
        % Check if the delays match
        if ~all(isequal(term1_jj.delay,term2_jj.delay))
            disp('The specified delays ',trm_name,'.delay do not match.')
            logval = false;
        end
        % Check if the integrals match
        if numel(term1_jj.I)~=numel(term2_jj.I)
            disp('The number of integrals in ',trm_name,' does not match.')
            logval = false;
        else
            for ll=1:numel(term1_jj.I)
                if ~all(isequal(term1_jj.I{ll},term2_jj.I{ll}))
                    disp('The domains of integration ',trm_name,'.I{',num2str(ll),'} do not match.')
                    logval = false;
                end
            end
        end
        % Check if the coefficients match
        if ~all(all(isequal(term1_jj.C,Cfctr*term2_jj.C)))
            if jj==1 && strcmp(obj,'BC')
                % If this is the first term in a BC, we allow the term in
                % PDE1 be a scaled version of that in PDE2, so long as all
                % other terms are scaled with this same factor
                if all(all(isequal(term1_jj.C,-term2_jj.C)))
                    % For polynomials, allow only scaling with -1
                    Cfctr = -1;
                elseif isa(term1_jj.C,'double') || isdouble(term1_jj.C) && ...
                        isa(term2_jj.C,'double') || isdouble(term2_jj.C)
                    % For constant coefficients, allow scaling with scalar
                    % factor.
                    CFCTR = double(term1_jj.C./term2_jj.C);
                    if numel(CFCTR)==1 || all(CFCTR(2:end)-CFCTR(1:end-1)==0)
                        Cfctr = CFCTR(1);
                    else
                        disp('The specified coefficients ',trm_name,'.C do not match.')
                        logval = false;
                    end
                else
                    disp('The specified coefficients ',trm_name,'.C do not match.')
                logval = false;
                end
            else
                disp('The specified coefficients ',trm_name,'.C do not match.')
                logval = false;
            end
        end  
    end
end

end