function out = int(objC, var, lim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that integrates a vector of states, x, in variable
% 'var' by 'lim' limit
% Input: 
% obj - state class object, x 
% var - pvar w.r.t. which differentiation is done
% lim - a 1x2 pvar object (or double object) 
% Output:
% out - state class object, \int_{lim{1}}^{lim{2}} x d(var)

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
% <<<<<<< Updated upstream
isdot_C = isdot(objC); isout_C=isout(objC); 
% for i=1:length(objA)
%     isdot_C = [isdot_C; objC(i).diff_order(1)*ones(subsref(objC(i),s),1)];
%     isout_C = [isout_C; strcmp(objC(i).type,'out')*ones(subsref(objC(i),s),1)];
% end
% =======
s.type = '.'; s.subs = 'veclength';
% isdot_C = []; isout_C=[]; 
% for i=1:length(objC)
%     isdot_C = [isdot_C; objC(i).diff_order(1)*ones(subsref(objC(i),s),1)];
%     isout_C = [isout_C; strcmp(objC(i).type,'out')*ones(subsref(objC(i),s),1)];
% end
% isdot_C = boolean(isdot_C); isout_C = boolean(isout_C);
% >>>>>>> Stashed changes
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
    opvar T; s.type = '.'; s.subs = 'veclength'; T.I = dom;
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