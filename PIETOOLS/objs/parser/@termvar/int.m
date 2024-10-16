function out = int(objA, var, lim)
if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

try
    intterms = double(objA.operator.R.R2)|double(objA.operator.R.R1); %R1 and R2 must be zero
    if any(intterms(:))
        error('Double integrals are currently not supported')
    end
catch
    error('Double integrals are currently not supported')
end

T = buildopvar('kernel',eye(length(objA)),'lim',lim,'var',var,'dom',objA.operator.dom);

out = termvar(T*objA.operator,objA.state);
end