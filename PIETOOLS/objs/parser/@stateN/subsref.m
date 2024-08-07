function out = subsref(obj,ref)
switch s(1).type
    case '.'
        if length(s) == 1
            % Implement obj.PropertyName
            if strcmp(s(1).subs,'type')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).type;
                end
            elseif strcmp(s(1).subs,'var')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).var;
                end
            elseif strcmp(s(1).subs,'len')
                for i=1:length(obj.len)
                    out(i,1) = obj(i).len;
                end
            elseif strcmp(s(1).subs,'diffOrder')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).diffOrder;
                end
            elseif strcmp(s(1).subs,'maxdiff')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).maxdiff;
                end
            elseif strcmp(s(1).subs,'dom')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).dom;
                end
            elseif strcmp(s(1).subs,'intLim')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).intLim;
                end
            elseif strcmp(s(1).subs,'multiplier')
                for i=1:length(obj.len)
                    out{i,1} = obj(i).multiplier;
                end
            elseif strcmp(s(1).subs,'statename')
                for i=1:length(obj.len)
                    out(i,1) = obj(i).statename;
                end
            else
                error('Inaccessible/unknown property for the state object');
            end
%         elseif length(s) == 2 && strcmp(s(2).type,'()')
%             % Implement obj.PropertyName(indices)
%             if strcmp(s(1).subs,'var')
%                 for i=1:length(obj.len)
%                     out(i,1) = obj(i).var(s(2).subs);
%                 end
%             elseif strcmp(s(1).subs,'diff_order')
%                 for i=1:length(obj.len)
%                     out(i,1) = obj(i).diff_order(s(2).subs);
%                 end
%             else
%                 error('Inaccessible/unknown property for the state object');
%             end
        else
        % Use built-in for any other expression
            error('Invalid indexing operation on state objec');
        end
    case '()'
        if isempty(s(1).subs{1})
                out= [];
                return;
        end
        if length(s) == 1
            % Implement obj(indices)
            out = obj(s(1).subs{1});
        elseif length(s) == 2 && strcmp(s(2).type,'.')
        % Implement obj(ind).PropertyName
            if strcmp(s(2).subs,'type')
               out = obj(s(1).subs{1}).type;
            elseif strcmp(s(2).subs,'var')
                out = obj(s(1).subs{1}).var;
            elseif strcmp(s(2).subs,'len')
                out = obj(s(1).subs{1}).len;
            elseif strcmp(s(2).subs,'diffOrder')
                out = obj(s(1).subs{1}).diffOrder;
            elseif strcmp(s(2).subs,'maxdiff')
                out = obj(s(1).subs{1}).maxdiff;
            elseif strcmp(s(2).subs,'dom')
                out = obj(s(1).subs{1}).dom;
            elseif strcmp(s(2).subs,'intLim')
                out = obj(s(1).subs{1}).intLim;
            elseif strcmp(s(2).subs,'multiplier')
                out = obj(s(1).subs{1}).multiplier;
            elseif strcmp(s(2).subs,'statename')
                out = obj(s(1).subs{1}).statename;
            else
                error('Inaccessible/unknown property for the state object');
            end
%         elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
%         % Implement obj(indices).PropertyName(indices)
%             out = obj(s(1).subs{1}).(s(2).subs(:));
%             out = out(s(3).subs{1});
        else
        % Use built-in for any other expression
            error('Invalid indexing operation on state objec');
        end
    otherwise
        error('Invalid indexing operation on state objec');
end
if isempty(out)
    out = 0;
end