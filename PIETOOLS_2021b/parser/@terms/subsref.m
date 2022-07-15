function out = subsref(obj,s)
switch s(1).type
    case '.'
        out = builtin('subsref',obj,s);
    case '()'
        out = terms(obj.operator(s(1).subs{1}),obj.statevec);
    otherwise
        error('Not a valid indexing expression');
end
end