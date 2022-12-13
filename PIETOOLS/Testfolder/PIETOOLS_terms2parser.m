function [PDE_p] = PIETOOLS_terms2parser(PDE_t)
% Takes
if ~PDE_t.is_initialized
    PDE_t = initialize(PDE_t);
end
if PDE_t.dim>1
    error('The function only allows 1D ODE-PDE systems at this time.')
end

pvar t s theta

vars_t = PDE_t.vars;
vars_p = [s,theta];

% Initialize the system structure
PDE_p = sys();
x = cell(numel(PDE_t.x),1);
w = cell(numel(PDE_t.w),1);
u = cell(numel(PDE_t.u),1);
z = cell(numel(PDE_t.z),1);
y = cell(numel(PDE_t.y),1);

% Declare the ODE and PDE state components
for jj=1:numel(PDE_t.x)
    if size(PDE_t.x{jj}.vars,1)==0
        x{jj} = state('ode',PDE_t.x{jj}.size);
    else
        x{jj} = state('pde',PDE_t.x{jj}.size);
    end
end
% Declare the inputs
for jj=1:numel(PDE_t.w)
    w{jj} = state('in',PDE_t.w{jj}.size);
end
for jj=1:numel(PDE_t.u)
    u{jj} = state('in',PDE_t.u{jj}.size);
end
% Declare the outputs
for jj=1:numel(PDE_t.z)
    z{jj} = state('out',PDE_t.z{jj}.size);
end
for jj=1:numel(PDE_t.y)
    y{jj} = state('out',PDE_t.y{jj}.size);
end

% % % Build the equations
% % First the ODE-PDE equations
eqns_x = build_eqns(PDE_t,'x',x,w,u,z,y,vars_t,vars_p,t);
% % Next the output equations
eqns_z = build_eqns(PDE_t,'z',x,w,u,z,y,vars_t,vars_p,t);
eqns_y = build_eqns(PDE_t,'y',x,w,u,z,y,vars_t,vars_p,t);
% % Finally the BCs
eqns_BC = build_eqns(PDE_t,'BC',x,w,u,z,y,vars_t,vars_p,t);

% eqns = [eqns_x; eqns_z; eqns_y; eqns_BC];
% PDE_p = addequation(PDE_p,eqns);
PDE_p = addequation(PDE_p,eqns_x);
if numel(PDE_t.z)>0
    PDE_p = addequation(PDE_p,eqns_z);
end
if numel(PDE_t.y)>0
    PDE_p = addequation(PDE_p,eqns_y);
end
if numel(PDE_t.BC)>0
    PDE_p = addequation(PDE_p,eqns_BC);
end

% Set controlled inputs.
for ii=1:numel(PDE_t.u)
    PDE_p = setControl(PDE_p,u{ii});
end
% Set observed outputs
for ii=1:numel(PDE_t.y)
    PDE_p = setObserve(PDE_p,y{ii});
end

end



function eqns = build_eqns(PDE_t,obj,x,w,u,z,y,vars_t,vars_p,t)
% % % Build the eqns defined by PDE_t.(obj) in parser format.

eqns = [];
for ii=1:numel(PDE_t.(obj))
    eqn = 0;
    for jj=1:numel(PDE_t.(obj){ii}.term)
        term_jj = PDE_t.(obj){ii}.term{jj};
        % Initialize the term: state var, disturbance, or input
        if isfield(term_jj,'x')
            is_x_Robj = true;
            trm = x{term_jj.x};
        elseif isfield(term_jj,'w')
            is_x_Robj = false;
            trm = w{term_jj.w};
        else
            is_x_Robj = false;
            trm = u{term_jj.u};
        end
        if is_x_Robj
            % Take desired derivative of the state
            if ~isempty(term_jj.D) && term_jj.D>0
                trm = diff(trm,vars_p(1),term_jj.D);
            end
            % Evaluate state at desired location
            if ~isempty(term_jj.loc) && (isa(term_jj.loc,'double') || isdouble(term_jj.loc))
                trm = subs(trm,vars_p(1),double(term_jj.loc));
            end
            % If taking partial integral, use dummy variable
            if ~isempty(term_jj.I) && ~isempty(term_jj.I{1}) && ~(isa(term_jj.I{1},'double') || isdouble(term_jj.I{1}))
                trm = subs(trm,vars_p(1),vars_p(2));
            end
        end
        % Account for the delay
        if term_jj.delay>0
            trm = subs(trm,t,t-term_jj.delay);
        end
        % Multiply with coefficient matrix
        if isa(term_jj.C,'polynomial')
            trm = subs(term_jj.C,vars_t',vars_p') * trm;
        else
            trm = term_jj.C * trm;
        end
        % Take the desired integral
        if ~isempty(term_jj.I) && ~isempty(term_jj.I{1})
            trm = int(trm,vars_p(1),term_jj.I{1});
        end
        
        % Add the newly constructed term to the equation
        if isa(eqn,'double') && eqn==0
            eqn = trm;
        else
            eqn = eqn + trm;
        end
    end
    
    % Add the equation to the list of eqs.
    if strcmp(obj,'x')
        eqns = [eqns;
            diff(x{ii},t)==eqn];
    elseif strcmp(obj,'z')
        eqns = [eqns;
                z{ii}==eqn];
    elseif strcmp(obj,'y')
        eqns = [eqns;
                y{ii}==eqn];
    elseif strcmp(obj,'BC')
        eqns = [eqns; -eqn];
    end
    
end

end
