function out = getPDEparams(pdeObj)
equations = pdeObj.equation;
statelist = pdeObj.states;
eqnNum = length(equations);
odeeqns = {}; pdeeqns = {}; outeqns = {}; bceqns = {};
for i=1:eqnNum
    eqnType{i} = identifyEqnType(equations(i));
    if strcmp(eqnType{i},'ode')
        odeeqns{end+1} = equations(i);
    elseif strcmp(eqnType{i},'pde')
        pdeeqns{end+1} = equations(i);
    elseif strcmp(eqnType{i},'out')
        outeqns{end+1} = equations(i);
    elseif strcmp(eqnType{i},'bc')
        bceqns{end+1} = equations(i);
    end
end

[~]=summarize_sys(obj);

out = pde_struct();

out.dim = 1;
out.domain = [0,1];
out.vars = [pvar('s'),pvar('theta')];

% extract ODE params

% extract output params

% extract BC params

% extract PDE params


end
