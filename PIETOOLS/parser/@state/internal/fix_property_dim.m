function obj = fix_property_dim(obj)
diffs = obj.diff_order{:}';
vars = obj.var_indep{:}';
if size(diffs,2)>size(vars,2)
    warning('Length of diff_order > length of var_indep. Appending new pvar variables to var_indep');
    n_vars = size(diffs);
    vars(1) = pvar('t');
    for i=size(obj.var_indep,2):n_vars
        vars(i) = pvar(['s',num2str(i-1)]);
    end
elseif size(diffs,2)<size(vars,2)
    warning('Length of diff_order < length of var_indep. Appending zeros to diff_order to match length');
    n_vars = size(diffs,2);
    for i=size(diffs,2):n_vars
        diffs(i) = 0;
    end
end
obj.diff_order = {diffs};
obj.var_indep = {vars};
end