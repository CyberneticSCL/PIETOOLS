function d3_ws(f)
% D3_WS(F), 11/27/2026: saves the variables of the base workspace to the
% .mat file F, class-free: every object becomes a struct of its properties,
% recursively (so that a file written with one class name loads in a tree
% with the other); graphics objects, function handles and LPI programs (a
% struct with field decvartable) are left out.
w = warning('off','all');
cl = onCleanup(@() warning(w));
names = evalin('base','who');
S = struct();
for k = 1:numel(names)
    v = evalin('base',names{k});
    if isa(v,'matlab.graphics.Graphics') || isa(v,'function_handle') || ...
            (isstruct(v) && isfield(v,'decvartable'))
        continue
    end
    S.(names{k}) = plain(v,0);
end
save(f,'-struct','S','-v7.3');
end

function v = plain(v,depth)
% V = PLAIN(V,DEPTH): objects to structs of their properties, recursively,
% to depth 12; the class name is kept in a field 'class_'.
if depth>12,    v = '<depth>';  return,     end
if isa(v,'matlab.graphics.Graphics') || isa(v,'function_handle') || isa(v,'handle')
    v = ['<',class(v),'>'];     return
end
if isobject(v) && ~isstring(v) && ~isdatetime(v) && ~isduration(v) && ~istable(v) && ~iscategorical(v)
    % builtin numel/subsref: the containers overload both (size counts
    % components; P(k) is a block slice), 11/27/2026
    c = class(v);   n = builtin('numel',v);
    if n==1
        v = plain(struct(v),depth+1);
    else
        s = cell(builtin('size',v));
        for k = 1:n,    s{k} = plain(struct(builtin('subsref',v,substruct('()',{k}))),depth+1);  end
        v = s;
    end
    if isstruct(v),     v.class_ = c;   end
elseif isstruct(v)
    f = fieldnames(v);
    for k = 1:numel(v)
        for j = 1:numel(f),     v(k).(f{j}) = plain(v(k).(f{j}),depth+1);   end
    end
elseif iscell(v)
    for k = 1:numel(v),     v{k} = plain(v{k},depth+1);     end
end
end
