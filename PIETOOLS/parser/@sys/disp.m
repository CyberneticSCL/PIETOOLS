function disp(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that displays equations in sys
% Input: 
% obj - sys class object
% Output:
% command line display of obj.params

if strcmp(obj.type,'pde')
    display_PDE(obj.params);
elseif strcmp(obj.type,'pie')
    disp_pie(obj.params);
else
    error('Display is supported only for sys objects of type PDE or PIE');
end
end