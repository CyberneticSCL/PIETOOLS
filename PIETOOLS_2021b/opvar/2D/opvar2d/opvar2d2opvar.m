function Pop_1d = opvar2d2opvar(Pop_2d,use_space_idx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop_1d = opvar2d2opvar(Pop_2d,use_space_idx) converts an opvar2d object 
% to a corresponding opvar object if possible.
% 
% INPUT
% - Pop_2d: opvar2d object, representing an operator mapping
%           R^n0 x L2^n1[s1] x L2^n2[s2] x L2^n3[s1,s2]
%               --> R^m0 x L2^m1[s1] x L2^m2[s2] x L2^m3[s1,s2].
%           Conversion is only possible if n3 and m3 are zero, and either
%           n1 and m1 or n2 and m2 are zero.
% - use_space_idx: char 'x' or 'y', or scalar 1 or 2, indicating whether
%           the output opvar should map L2^n1[s1] --> L2^m1[s2] or 
%           L2^n2[s2] --> L2^m2[s2]. Only really useful if both n1 and m1,
%           and n2 and m2 are zero, in which case it is unclear which space
%           the output opvar object should map to/from.
% 
% OUTPUT
% - Pop_1d: opvar object, representing an operator mapping
%           R^n0 x L2^n1[s1] --> R^m0 x L2^m1[s1] if n2=m2=0, or
%           R^n0 x L2^n2[s2] --> R^m0 x L2^m2[s2] if n1=m1=0.
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, - 08/12/2022.


% Conversion from 2d opvar to 1d opvar is only supported if the 2d operator
% only maps functions in at most one spatial variable.
Pdim = Pop_2d.dim;
use_space = [any(Pdim(2,:)) ; any(Pdim(3,:))];  % Does the operator map functions in x or in y?
if any(Pdim(4,:)) || any(Pdim(2,:)) && any(Pdim(3,:))
    error(['Conversion of the 2D operator to a 1d operator is not supported:',...
            ' the 2d operator maps functions in more than 1 spatial variable.'])
elseif nargin==2
    if ischar(use_space_idx)
        if strcmp(use_space_idx,'x') || strcmp(use_space_idx,Pop_2d.var1(1).varname{1})
            use_space(1) = true;
        elseif strcmp(use_space_idx,'y') || strcmp(use_space_idx,Pop_2d.var1(2).varname{1})
            use_space(2) = true;
        else
            error(['The second argument to the function is not appropriate.',...
                    ' Please specify on which of the spatial variables of the 2D operator should carry over to the 1D operator.'])
        end
    else
        use_space(use_space_idx) = true;
    end
    if all(use_space)
        error(['The specified variable does not match that of the functions the 2D operator maps (to).'])
    end
end
if ~any(use_space)
    warning(['It is unclear what function space the 1D operator should map to/from;',...
                ' Continuing with the first 1D function space the 2D operator maps.'])
    use_space(1) = true;
end

% Initialize a 1D operator on the desired domain.
Pop_1d = opvar();
Pop_1d.I = Pop_2d.I(use_space,:);
Pop_1d.var1 = Pop_2d.var1(use_space,:);
Pop_1d.var2 = Pop_2d.var2(use_space,:);

% Set the dimensions and parameters of the operator.
Pop_1d.dim = Pdim([true; use_space; false],:);
if use_space(1)
    Pop_1d.P = Pop_2d.R00;
    Pop_1d.Q1 = Pop_2d.R0x;
    Pop_1d.Q2 = Pop_2d.Rx0;
    Pop_1d.R.R0 = Pop_2d.Rxx{1};
    Pop_1d.R.R1 = Pop_2d.Rxx{2};
    Pop_1d.R.R2 = Pop_2d.Rxx{3};
else
    Pop_1d.P = Pop_2d.R00;
    Pop_1d.Q1 = Pop_2d.R0y;
    Pop_1d.Q2 = Pop_2d.Ry0;
    Pop_1d.R.R0 = Pop_2d.Ryy{1};
    Pop_1d.R.R1 = Pop_2d.Ryy{2};
    Pop_1d.R.R2 = Pop_2d.Ryy{3};
end

end