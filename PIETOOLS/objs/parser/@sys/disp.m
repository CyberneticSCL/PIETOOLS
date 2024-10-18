function disp(obj)
if strcmp(obj.type,'pde')||strcmp(obj.type,'pie')
    display(obj.params);
else
    error('Display is supported only for sys objects of type PDE or PIE');
end
end