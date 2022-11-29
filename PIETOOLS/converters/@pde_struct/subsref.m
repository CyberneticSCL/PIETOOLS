function val = subsref(PDE,prop)
% SUBSREF allows properties and elements of the "pde_struct" object PDE to
% be retrieved. The function is automatically called when calling
% e.g. PDE(i,j) or PDE.obj for a "pde_struct" PDE.
%
% INPUTS:
% - PDE:    A pde_struct object of which to extract a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PDE
%           structure to retrieve.
%
% OUTPUTS:
% - val:    The value of the field of "PDE" specified by "prop".
%
% EXAMPLE: calling PDE.x{1}.term{2}.C(3), we have
%   PDE_in = PDE;
%   prop(1).type = '.',     prop(1).subs = 'x';
%   prop(2).type = '{}',    prop(2).subs = {[1]};
%   prop(3).type = '.',     prop(3).subs = 'term';
%   prop(4).type = '{}',    prop(4).subs = {[2]};
%   prop(5).type = '.',     prop(5).subs = 'C';
%   prop(6).type = '()',    prop(6).subs = {[3]};
%
% NOTES:
% For the most part, this function just calls the Matlab built-in subsref
% function. However, in the display, certain coefficients are shown as
% C_{ij}, so we added a feature to extract these coefficients by calling
% PDE.C{i,j}.
%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 11/29/2022
%

% At this point, we just use the built-in subsasgn function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'()') || strcmp(prop(1).type,'{}')
    error(['''()'' and ''{}'' type subsref are currently not supported for "pie_struct" class objects.'])
else
    % Allow coefficients to be retrieved by calling PDE.C{i,j}
    if strcmp(prop(1).type,'.') && strcmpi(prop(1).subs,'C')
        if ~strcmp(prop(2).type,'{}')
            error('To extract the coefficients C_{ij} appearing in term j of equation i, call "PDE.C{i,j}".')
        elseif numel(prop(2).subs)<2
            error('Both an equation number i and term number j must be specified to extract coefficients "PDE.C{i,j}".')
        elseif numel(prop(2).subs)>2
            error('Only an equation number i and term number j can be specified when extracting coefficients "PDE.C{i,j}".')
        else
            ii = prop(2).subs{1};
            jj = prop(2).subs{2};
        end
        % Distinguish between dynamics, outputs, or BCs
        if ii <= 0 || ii~=round(ii) || ii~=real(ii)
            error('Equation number i must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        elseif ii <= numel(PDE.x)
            % Coefficients in PDE dynamics are requested.
            eq = PDE.x{ii};
            eqname = ['x{',num2str(ii),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y)
            % Coefficients in observed output equations are requested.
            eqnum = ii - numel(PDE.x);
            eq = PDE.y{eqnum};
            eqname = ['y{',num2str(eqnum),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y) + numel(PDE.z)
            % Coefficients in regulated output equations are requested.
            eqnum = ii - numel(PDE.x) - numel(PDE.y);
            eq = PDE.z{eqnum};
            eqname = ['z{',num2str(eqnum),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y) + numel(PDE.z) + numel(PDE.BC)
            % Coefficients in BCs are requested.
            eqnum = ii - numel(PDE.x) - numel(PDE.y) - numel(PDE.z);
            eq = PDE.BC{eqnum};
            eqname = ['BC{',num2str(eqnum),'}'];
        else
            error('The proposed equation index exceeds the number of equations in the system.')
        end
        % Extract coefficients from the desired term.
        if jj <= 0 || jj~=round(jj) || jj~=real(jj)
            error('Term number j must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        elseif jj>numel(eq.term)
            error(['The proposed term index exceeds the number of terms in the equation "PDE.',eqname,'".'])
        else
            val = eq.term{jj}.C;
        end
    else
        % Otherwise, we have nothing fancy implemented at this time.
        val = builtin('subsref',PDE,prop);
    end
end

end