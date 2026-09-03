function info = verify(Mop)
% true = verify(MOP) determines whether the dimensions, domains, and
% variable names in a mopvar object are consistent.
%
% INPUTS
% - Mop:    A 'mopvar' object
%
% OUTPUTS
% - info.true = 1 if the Mopvar object is consistent, 0 otherwise   
%   info.flags={} if mopvar is consistent and a cell list of inconsistent values otherwise 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - verify(mopvar)
%
% Copyright (C) 2026 PIETOOLS Team
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
% MP, 01/19/2026: Initial coding

[M,N]=size(Mop);
ref = cellfun(@(f) C{1}.(f), fields, 'UniformOutput', false);
for i=1:M % check row consistency
    row=Mop{i,:};
    fields = {'dim(1)','vars(1)'};
    ref = cellfun(@(f) r{1}.(f), fields, 'UniformOutput', false);
    tfr{i} = all(cellfun(@(s) ...
    isequal(cellfun(@(f) s.(f), fields, 'UniformOutput', false), ref), row));
    % Note, check for and eliminate should eliminate empty elements
    % Cdom=cellfun(@(s) s.dom, row);
    % Cdim=cellfun(@(s) s.dim(1), row);
    % Cvars=cellfun(@(s) s.vars(2), row);
end
for i=1:N % check column consistency
    col=Mop{:,i};
    fields = {'dim(2)','vars(2)'};
    ref = cellfun(@(f) r{1}.(f), fields, 'UniformOutput', false);
    tfc{i} = all(cellfun(@(s) ...
    isequal(cellfun(@(f) s.(f), fields, 'UniformOutput', false), ref), col));
    % Note, check for and eliminate should eliminate empty elements
end


end