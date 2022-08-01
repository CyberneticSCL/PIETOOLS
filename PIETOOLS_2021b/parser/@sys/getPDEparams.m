function out = getPDEparams(pdeObj)
equations = pdeObj.equation;
statelist = pdeObj.states;
eqnNum = length(equations);

% build a dictionary of state names for future reference
odeNames = statelist(find(strcmp(statelist.type,'ode'))).statename;
pdeNames = statelist(find(strcmp(statelist.type,'pde'))).statename;
xNames = [odeNames; pdeNames];
outNames = statelist(find(strcmp(statelist.type,'out'))).statename;
inNames = statelist(find(strcmp(statelist.type,'in'))).statename;

% find the term location of states in the equation
veclen_sum = [0;cumsum(equations.statevec.veclength)]+1;
odeidx = find(strcmp(equations.statevec.type,'ode')); % find all terms with ode states
pdeidx = find(strcmp(equations.statevec.type,'pde')); % find all terms with pde states
outidx = find(strcmp(equations.statevec.type,'out')); % find all terms with output 
inidx = find(strcmp(equations.statevec.type,'in')); % find all terms with input



% start parsing the system equations
out = pde_struct();

% set known global properties 
out.dim = 1;
out.domain = [0,1];
out.vars = [pvar('s'),pvar('theta')];

% set state properties
out.x = cell(length(xNames),1);
out.w = cell(length(inNames)-sum(pdeObj.ControlledInputs),1);
out.u = cell(sum(pdeObj.ControlledInputs),1);
out.z = cell(length(outNames)-sum(pdeObj.ObservedOutputs),1);
out.y = cell(sum(pdeObj.ObservedOutputs),1);

isdot_A = []; isout_A=[]; 
for i=1:length(equations.statevec)
    isdot_A = [isdot_A; equations.statevec(i).diff_order(1)*ones(equations.statevec(i).veclength,1)];
end



for i=1:eqnNum
    row = equations(i);
    if ~isequal(strcmp(equations.statevec.type,'out').*row.operator.R0,0) % equation has outputs
        
    elseif ~isequal(isdot_A.*any(strcmp(equations.statevec.type,{'ode','pde'})).*row.operator.R0,0)% equation has dynamics
        
    else % boundary conditions
        
    end
end

end
