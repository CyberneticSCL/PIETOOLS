function display(sys)
if strcmp(sys.type,'pde')
    display_PDE(sys.params);
elseif strcmp(sys.type,'pie')
    display(sys.params);
else
    error('Display is supported only for sys objects of type PDE or PIE');
end
end