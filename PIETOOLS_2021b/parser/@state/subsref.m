function out = subsref(obj,s)
switch s(1).type
    case '.'
        out = builtin('subsref',obj,s);
    case '()'
        if length(s)==1
            out = state(obj.type(s(1).subs{1}),obj.veclength(s(1).subs{1}),obj.var(s(1).subs{1}),obj.statename(s(1).subs{1}));
        else
            out = state(obj.type(s(1).subs{1}),obj.veclength(s(1).subs{1}),obj.var(s(1).subs{1}),obj.statename(s(1).subs{1}));
            out = out.(s(2).subs);
        end
    case '{}'
        error('Not a valid indexing expression')
end
end