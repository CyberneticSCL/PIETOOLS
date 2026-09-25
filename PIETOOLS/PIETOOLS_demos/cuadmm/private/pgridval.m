function v = pgridval(p,vars,vals)
% Evaluate a scalar polynomial at arrays of values, reading the multipoly
% internals directly rather than calling subs() once per quadrature node.  NOT named peval: the polynomial class already defines that method, which would win dispatch.
%   vars - cell of variable names OR pvar objects
%   vals - cell of arrays, same shape, one per variable
if isa(p,'double')
    v = p*ones(size(vals{1}));  return
end
% Callers pass pvar OBJECTS as the variable list; degmat is indexed by NAME,
% so normalise first.  Without this, lookup silently fails for any kernel that
% actually depends on a variable (a constant kernel returns early and hides it).
for i = 1:numel(vars)
    if isa(vars{i},'polynomial'), vars{i} = vars{i}.varname{1}; end
end
nm = p.varname;  dm = p.degmat;  cf = p.coefficient;
if isempty(nm) || isempty(dm)
    v = full(sum(cf))*ones(size(vals{1}));  return
end
sz = size(vals{1});
v  = zeros(sz);
for k = 1:size(dm,1)
    term = full(cf(k))*ones(sz);
    for j = 1:numel(nm)
        e = full(dm(k,j));
        if e == 0, continue; end
        idx = find(strcmp(nm{j},vars),1);
        if isempty(idx)
            error('pgridval:unbound','polynomial depends on unbound variable %s',nm{j});
        end
        term = term .* (vals{idx}.^e);
    end
    v = v + term;
end
end
