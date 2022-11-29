function out = subsref(obj,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method for sub indexing of obj. Allowed subindexing:
% obj.(property), obj(a:b)
% Input: 
% obj - state class object
% s - indexing structure
% Output:
% out - obj(s) or obj.s, state class object
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
switch s(1).type
    case '.'
        if length(s) == 1
            % Implement obj.PropertyName
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
        elseif length(s) == 2 && strcmp(s(2).type,'()')
            % Implement obj.PropertyName(indices)
            if strcmp(s(1).subs,'var')
                for i=1:length(obj)
                    out(i,1) = obj(i).var(s(2).subs);
                end
            elseif strcmp(s(1).subs,'diff_order')
                for i=1:length(obj)
                    out(i,1) = obj(i).diff_order(s(2).subs);
                end
            else
                error('Inaccessible/unknown property for the state object');
            end
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
            elseif strcmp(s(2).subs,'veclength')
                out = obj(s(1).subs{1}).veclength;
            elseif strcmp(s(2).subs,'diff_order')
                out = obj(s(1).subs{1}).diff_order;
            elseif strcmp(s(2).subs,'statename')
                out = obj(s(1).subs{1}).statename;
            else
                error('Inaccessible/unknown property for the state object');
            end
        elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
        % Implement obj(indices).PropertyName(indices)
            out = obj(s(1).subs{1}).(s(2).subs(:));
            out = out(s(3).subs{1});
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