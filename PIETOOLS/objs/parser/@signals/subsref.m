function out = subsref(obj,s)
switch s(1).type
    case '.'
        if length(s) == 1
            % Implement obj.PropertyName
            prop = s(1).subs;
            n = length(obj);
            out = [];
            for k = 1:n
                out = [out; obj(k).(prop)];
            end
        else
            error('Invalid indexing expression for signals');
        end
    case '()'
        if length(s) == 1
            % Implement obj(indices)
            idx = s(1).subs{1};
            n = length(idx);
            out = [];
            for k = 1:n
                out = [out; obj(k)];
            end
        elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement obj(ind).PropertyName
            tmp = subsref(obj,s(1));
            out = subsref(tmp,s(2));
        elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
        % Implement obj(indices).PropertyName(indices)
            tmp = subsref(obj,s(1));
            tmp = subsref(tmp,s(2));
            out = subsref(tmp,s(3));
        else
            % Use built-in for any other expression
            error('Invalid indexing expression for signals');
        end
    otherwise
        error('Invalid indexing expression for signals');
end