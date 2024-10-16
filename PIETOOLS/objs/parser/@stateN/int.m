function obj = int(obj,var,lim)
if (length(lim)~=2)||(~isa(lim,'double')&&~isa(lim,'polynomial'))
    error('Third argument limits must be a 1x2 vector of double or pvars');
end
if isequal(var,pvar('t'))
    error('Integration with respect to time is not currently supported');
end
for i=1:length(obj.len)
    tmpObj = obj(i);
    [loc,logval] = ismember(var,tmpObj.var); % find if var is present in obj(i)
    if logval
        if all(~isequal(tmpObj.dom,lim))
            error("Cannot integrate "+i+"th component of state since the integral is neither an upper nor lower integral on state");
        end
        tmpObj.intLim = {lim};
    else % integration variable not present, multiply by interval length
        tmpMul = lim(2)-lim(1);
        tmpObj.multiplier = tmpMul*tmpObj.multiplier;
    end
    obj(i) = tmpObj;
end
end