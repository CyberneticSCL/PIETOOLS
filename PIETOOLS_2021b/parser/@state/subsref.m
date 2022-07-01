function out = subsref(obj,s)
switch s(1).type
    case '.'
        if strcmp(s(1).subs,'type')
            for i=1:length(obj)
                out{i,1} = obj(i).type;
            end
        elseif strcmp(s(1).subs,'var')
            for i=1:length(obj)
                out{i,1} = obj(i).var;
            end
        elseif strcmp(s(1).subs,'veclength')
            for i=1:length(obj)
                out(i,1) = obj(i).veclength;
            end
        elseif strcmp(s(1).subs,'delta_val')
            for i=1:length(obj)
                out{i,1} = obj(i).delta_val;
            end
        elseif strcmp(s(1).subs,'diff_order')
            for i=1:length(obj)
                out{i,1} = obj(i).diff_order;
            end
        elseif strcmp(s(1).subs,'statename')
            for i=1:length(obj)
                out(i,1) = obj(i).statename;
            end
        else
            error('Inaccessible/unknown property for the state object');
        end
    case '()'
        out = subsref('builtin',obj,s);
    case '{}'
        error('Not a valid indexing expression')
end
end