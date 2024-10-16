function out = int(objC, var, lim)
dom = objC.dom;

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end
if ~poly2double(lim(2))&&~poly2double(lim(1))
    error("Integral should have at least one limit of integration as "+num2str(dom(1))+" or "+num2str(dom(2)));
end
if poly2double(lim(2))&&(double(lim(2))~=dom(2))
    error("Upper limit of integration can only be a pvar variable or "+num2str(dom(2)));
end
if poly2double(lim(1))&&(double(lim(1))~=dom(1))
    error("Lower limit of integration can only be a pvar variable or "+num2str(dom(1)));
end

out = [];
for i=1:length(objC)
idx = find(isequal(objC(i).var,var));
if isempty(idx) % state drops out of integration
    K = int(1,var,lim(1),lim(2));
    out = [out;mtimes(K,objC(i))];
else
    mi T; s.type = '.'; s.subs = 'len'; T.dom = dom;
    T.kernel = eye(subsref(objC(i),s));
    T.lim{idx} = lim;
    objC(i).var(idx) = T.var1;
    out = [out;termvar(T,objC)];
end
end
end