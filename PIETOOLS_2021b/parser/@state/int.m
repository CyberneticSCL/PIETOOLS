function out = int(objC, var, lim)

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

if ~poly2double(lim(2))&&~poly2double(lim(1))
    error('Integral should have at least one limit of integration as 0 or 1');
end
if poly2double(lim(2))&&(double(lim(2))~=1)
    error('Upper limit of integration can only be a "pvar" variable or 1');
end
if poly2double(lim(1))&&(double(lim(1))~=0)
    error('Lower limit of integration can only be a "pvar" variable or 0');
end
isdot_C = []; isout_C=[]; 
for i=1:length(objA)
    isdot_C = [isdot_C; objC(i).diff_order(1)*ones(subsref(objC(i),s),1)];
    isout_C = [isout_C; strcmp(objC(i).type,'out')*ones(subsref(objC(i),s),1)];
end
if any((isdot_C|isout_C))
    error("Integration of vectors with outputs or time-derivative of state is not allowed");
end



out = [];
for i=1:length(objC)
idx = find(isequal(objC(i).var,var));
if isempty(idx) % state drops out of integration
    K = int(1,var,lim(1),lim(2));
    out = [out;mtimes(K,objC(i))];
else
    opvar T; s.type = '.'; s.subs = 'veclength';
    if poly2double(lim(2))
            T.R.R2 = eye(subsref(objC(i),s));
    end
    if poly2double(lim(1))
            T.R.R1 = eye(subsref(objC(i),s));
    end
    objC(i).var(idx) = T.var1;
    out = [out;terms(T,objC)];
end
end
end