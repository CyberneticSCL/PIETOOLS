function summarize_sys(obj)

if strcmp(obj.type,'pde')
    statelist = obj.states;
    ode_list = statelist(find(strcmp(statelist.type,'ode')));
    n_ode = sum(ode_list.veclength);
    pde_list = statelist(find(strcmp(statelist.type,'pde')));
    n_pde = sum(pde_list.veclength);

    disp('----- Summary of identified states/dependent variables -----');
    fprintf('Number of ODE states: %d\n',n_ode);
    fprintf('Number of PDE states: %d\n',n_pde);
    fprintf('Number of regulated outputs: %d\n',length(obj.ObservedOutputs));
    fprintf('Number of Observed outputs: %d\n',sum(obj.ObservedOutputs));
    fprintf('Number of Disturbance inputs: %d\n',length(obj.ControlledInputs));
    fprintf('Number of Control inputs: %d\n',sum(obj.ControlledInputs));
    
    disp('----- Summary of identified equations -----');
end

end