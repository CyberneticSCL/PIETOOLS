function summarize_sys(obj)

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
    
    equations = pdeObj.equation;
    eqnNum = length(equations);
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
    
    disp('----- Summary of identified equations -----');
    fprintf('Number of ODE dynamics: %d\n',length(odeeqns));
    fprintf('Number of PDE dynamics: %d\n',length(pdeeqns));
    fprintf('Number of output equations: %d\n',length(outeqns));
    fprintf('Number of Boundary conditions: %d\n',length(bceqns));
end

end