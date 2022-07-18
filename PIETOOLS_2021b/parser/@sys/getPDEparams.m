function out = getPDEparams(pdeObj)
equations = pdeObj.equation;
statelist = pdeObj.states;
eqnNum = length(equations);
odeeqns = {}; pdeeqns = {}; outeqns = {}; bceqns = {};
for i=1:eqnNum
    eqnType{i} = identifyEqnType(equations(i));
end

[~]=summarize_sys(obj,eqnType);

veclen_sum = [0;cumsum(equations.statevec.veclength)]+1;
odeidx = find(strcmp(equations.statevec.type,'ode')); % find all terms with ode states
pdeidx = find(strcmp(equations.statevec.type,'pde')); % find all terms with pde states
outidx = find(strcmp(equations.statevec.type,'out')); % find all terms with output 
inidx = find(strcmp(equations.statevec.type,'in')); % find all terms with input


% start parsing the system equations
out = pde_struct();

out.dim = 1;
out.domain = [0,1];
out.vars = [pvar('s'),pvar('theta')];

odeNames = statelist(find(strcmp(statelist.type,'ode'))).statename;
pdeNames = statelist(find(strcmp(statelist.type,'pde'))).statename;
xXDict = containers.map([odeNames;pdeNames],1:length(odeNames)+length(pdeNames));
outNames = statelist(find(strcmp(statelist.type,'out'))).statename;
YZDict = containers.map(outNames,1:length(outNames));
inNames = statelist(find(strcmp(statelist.type,'in'))).statename;
WUDict = containers.map(inNames,1:length(inNames));

% extract ODE params
tmpodeidx = odeidx;
for i=1:length(equations)
    eqnTemp = equations(i);
    if strcmp(eqnType{i},'ode') % if ode equation extract ODE parameters
        for j=1:length(odeidx)% find the time derivative term
            if equations.statevec(odeidx(j)).diff_order{1} 
                derivLoc = odeidx(j);
                derivName = equations.statevec(derivLoc).statename;
            end
        end
        % Now go through each term in the equation to populate parameters
        tempOp = equations.operator.R.R0; % ODEs can only have opvar with P terms 
        tempOq = equations.operator.R.R1; % ODEs can only have opvar with Q1 terms,  R1=R2
        for k=1:length(veclen_sum)-1
            if derivLoc==k % dont add anything since its the derivative term
                out.x{xXDict(derivName)}.term{k}.C = [];
                out.x{xXDict(derivName)}.term{k}.x = [];
                continue
            end
            if strcmp(equations.statevec(k).type,'ode') % contribution from ode state
                out.x{xXDict(derivName)}.term{k}.x = equations.statevec(k).statename;
                out.x{xXDict(derivName)}.term{k}.C = [out.x{xXDict(j)}.term{k}.C;tempOp(1,veclen_sum(k):veclen_sum(k))];
            end
            if strcmp(equations.statevec(k).type,'in') % contribution from inputs
                if obj.ControlledInputs(something) % separate between w and u
                    out.x{xXDict(derivName)}.term{k}.u = equations.statevec(k).statename;
                    out.x{xXDict(derivName)}.term{k}.C = [out.x{xXDict(j)}.term{k}.C;tempOp(1,veclen_sum(k):veclen_sum(k))];
                else
                    out.x{xXDict(derivName)}.term{k}.w = equations.statevec(k).statename;
                    out.x{xXDict(derivName)}.term{k}.C = [out.x{xXDict(j)}.term{k}.C;tempOp(1,veclen_sum(k):veclen_sum(k))];
                end
            end
            if strcmp(equations.statevec(k).type,'pde') % contribution from pde state (integral or boundary term)
            end
        end
    end
    if strcmp(eqnType{i},'pde') % if pde equation extract PDE parameters
    end
    if strcmp(eqnType{i},'out') % if Output equation extract Z/Y parameters
    end
    if strcmp(eqnType{i},'bc') % if BC extract BC parameters
    end
end


end
