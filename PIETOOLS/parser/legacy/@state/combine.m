function [out, varargout] = combine(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that combines a vector of states, x, to get (x_min,A),
% where x_min is the unique set of vectors in, x, and x = A*x_min
% Input: 
% varargin - states to combine, x{1}...x{N}
% Output:
% out - vector of unique states, x{1}...x{M}
% varargout - Matrices A{1}..A{M}, where x{i} = A{i}*out
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
% MMP, 11/09/2026: Inlining pass (maintainer, 11/09/2026): the local
%                  findPermutation (one caller) inlined into the output loop.
out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp)
        s.type = '()'; s.subs = {j};
        temp = subsref(tmp,s); 
        if ~ismember(temp,out)
            out = [out; temp];
        end
    end
end

varargout = cell(1,nargin);

for i=1:nargin
%   varargout{i} = findPermutation(out,varargin{i});                        % MMP, 11/09/2026 (was)
    % the permutation Pi with varargin{i} = Pi*out (the former local        % MMP, 11/09/2026
    % FINDPERMUTATION): one identity block per state of varargin{i}, placed % MMP, 11/09/2026
    % at the position of that state in out                                  % MMP, 11/09/2026
    ss.type = '.'; ss.subs = 'veclength';                                   % MMP, 11/09/2026
    Pi = zeros(sum(subsref(varargin{i},ss)),sum(subsref(out,ss)));          % MMP, 11/09/2026
    blen = [0;cumsum(subsref(varargin{i},ss))]+1; alen = [0;cumsum(subsref(out,ss))]+1; % MMP, 11/09/2026
    [~,pidx] = ismember(varargin{i},out);                                   % MMP, 11/09/2026
    for ib=1:length(blen)-1                                                 % MMP, 11/09/2026
        vlen = subsref(out,ss);                                             % MMP, 11/09/2026
        Pi(blen(ib):blen(ib+1)-1,alen(pidx(ib)):alen(pidx(ib)+1)-1) = eye(vlen(pidx(ib))); % MMP, 11/09/2026
    end                                                                     % MMP, 11/09/2026
    varargout{i} = Pi;                                                      % MMP, 11/09/2026
end
end
