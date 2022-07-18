function out = identifyEqnType(equation)
statevec = equation.statevec;
classifiedeqns = zeros(length(equation),1);
for i=1:length(statevec)
    if strcmp(statevec(i).type,'ode')
    elseif strcmp(statevec(i).type,'pde')
    elseif strcmp(statevec(i).type,'out')
    end
end

if any(strcmp(statevec.state.type,'out')) %if equation has output
    out = 'out';
elseif any(cellfun(@(x) x(1)~=0, statevec.diff,'un',1))%if term has derivative of time
    if strcmp(statevec.state(find(cellfun(@(x) x(1)~=0, statevec.diff,'un',1))).type,'ode')% ode with derivative of time
        out = 'ode';
    else %pde
        out = 'pde';
    end
else % boundary condition
    out = 'bc';
end
end