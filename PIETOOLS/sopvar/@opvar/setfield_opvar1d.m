function P = setfield_opvar1d(P,fname,val)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = setfield_opvar1d(P,fname,val) sets the field fname of the PI operator
%   P: R^m x L2^n[a,b] -> R^p x L2^q[a,b]
% as if P were a 2025 'opvar' object, acting as
%   [ P    Q1 ] [ u ]       (P u + int_a^b Q1(s) x(s) ds )
%   [ Q2   R  ] [ x ] (s) = (Q2(s) u + R0(s) x(s) + int_a^s R1(s,t) x(t) dt
%                                    + int_s^b R2(s,t) x(t) dt ),
% so that P.(fname) then returns val.
%
% INPUT
% P:      opvar class object with the 1-D field layout: at most
%         one spatial variable, and every block row (column) of R before
%         every block row (column) of L2
% fname:  char or string, a field name of the 1-D layout: 'P', 'Q1', 'Q2',
%         'R', 'R.R0', 'R.R1', 'R.R2', 'I', 'var1' or 'var2'
% val:    the value of the field:
%         'P':      p x m 'double', constant;
%         'Q1':     p x n 'double' or 'polynomial' in var1 (or in var2);
%         'Q2':     q x m 'double' or 'polynomial' in var1;
%         'R':      struct with fields R0, R1, R2, each q x n, R0 in var1,
%                   R1 and R2 in var1 and var2;
%         'R.R0', 'R.R1', 'R.R2': one of these;
%         [] for a zero component;
%         'I':      1 x 2 array [a,b];
%         'var1':   1 x 1 'polynomial', the spatial variable s;
%         'var2':   1 x 1 'polynomial', only s_dum (var1 with '_dum') is
%                   accepted, the name a 'sopvar' gives the dummy variable;
%
% OUTPUT
% P:      the opvar with the field set
%
% NOTES:
% - The fields are those of the 2-D layout on the one variable x = s: P is R00,
%   Q1 is R0x, Q2 is Rx0, and R is Rxx with R.R0-R.R2 its parameters 1-3,
%   and each is set as that field of the 2-D layout. A component on a space of
%   dimension 0 adds a block row (column) for it, the other components in
%   it being zero.
% - Without a spatial variable, s is the variable the kernel names, else
%   's', on [0,1], as in the 2025 class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - setfield_opvar1d (opvar)
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/15/2026: Initial coding;

% Check the structure of P and the field name; default variable name s
if ~is_opvar1d(P)
    error('opvar:notOpvar1d',['The operator does not have the 1-D field layout (P, Q1, Q2, R) of an operator '...
        'R^m x L2^n[a,b] -> R^p x L2^q[a,b]; see ''is_opvar1d''.'])
end
if ~(ischar(fname) || (isstring(fname) && isscalar(fname)))
    error('opvar:badField','The field name must be a char or string.')
end
fname = char(fname);
opts = struct('vars',{{'s'}},'label',fname);

% The field of the 2-D layout on x = s; one variable for I, var1, var2
switch fname
    case 'P',       P = setfield_opvar2d(P,'R00',val,[],opts);
    case 'Q1',      P = setfield_opvar2d(P,'R0x',val,[],opts);
    case 'Q2',      P = setfield_opvar2d(P,'Rx0',val,[],opts);
    case 'R'
        % The three kernels, as the parameters {R0;R1;R2} of Rxx
        if ~isstruct(val) || ~isscalar(val) || ~isempty(setxor(fieldnames(val),{'R0';'R1';'R2'}))
            error('opvar:badValue','''R'' must be a struct with fields R0, R1 and R2.')
        end
        P = setfield_opvar2d(P,'Rxx',{val.R0;val.R1;val.R2},[],opts);
    case {'R.R0','R.R1','R.R2'}
        % One kernel: parameter 1, 2 or 3 of Rxx
        P = setfield_opvar2d(P,'Rxx',val,str2double(fname(end))+1,opts);
    case {'I','var1','var2'}
        % One spatial variable: one row of I, one variable
        if (strcmp(fname,'I') && size(val,1)~=1) || (~strcmp(fname,'I') && numel(val)~=1)
            error('opvar:badValue','The 1-D field layout has one spatial variable: ''%s'' must have one row.',fname)
        end
        P = setfield_opvar2d(P,fname,val,[],opts);
    case 'dim'
        error('opvar:setDim','''dim'' follows from the components; set those instead.')
    otherwise
        error('opvar:badField',['''%s'' is not in the 1-D field layout; use one of '...
            '''P'', ''Q1'', ''Q2'', ''R'', ''R.R0'', ''R.R1'', ''R.R2'', ''I'', '...
            '''var1'', ''var2''.'],fname)
end

end
