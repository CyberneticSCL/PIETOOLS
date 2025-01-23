function eq_opts = get_eq_opts_2D(Qop,eq_opts,ztol)
% EQ_OPTS = GET_EQ_OPTS_2D(QOP,EQ_OPTS,ZTOL) takes an operator
% "Qop" for which to declare an equality constraint "Qop==Qeop", and
% determines what options to use in declaring the object "Qeop" using
% 'poslpivar2d'.
%
% INPUT
% - Qop:    'opvar2d' or 'dopvar2d' object for which to enforce an
%           equality constraint Qop==Qeop for an operator decision variable
%           Qeop;
% - eq_opts:    'struct' specifying desired settings for declaring an 
%               'dopvar2d' object Qeop to enforce the constraint Qop==Qeop.
%               The struct should include a field 'exclude'
%               specifying whether any components of the operator should be
%               excluded, and a field 'sep' specifying whether the operator
%               should be separable (see the 2D settings files for more
%               details).
% - ztol:       Scalar value of type 'double' specifying a tolerance below
%               which coefficients are assumed to be zero, in checking
%               whether the operator is separable.
%
% OUTPUT
% - eq_opts:    'struct' with updated settings for declaring a `dopvar2d'
%               object Qeop for enforcing the equality constraint Qop==Qeop.
%               Specifically, if the function finds certain components of
%               Qop to be zero, the associated fields in eq_opts.exclude
%               will be set to true. Similarly, if the function finds Qop
%               to be seperable, the associated field eq_opts.sep will be
%               set to true accordingly. 
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
% Initial coding DJ - 07/31/2024

% Determine the required dimension of the operator Qeop in the equality
%       Qop==Qeop;
Qdim = Qop.dim(:,1);

% % Determine which components of Qeop can be excluded.
eq_opts_exc = eq_opts.exclude;
eq_opts_sep = eq_opts.sep;
% First check the map R-->R
if Qdim(1)
    % Check if Qop.R00 is zero.
    if isempty(Qop.R00) || (isdouble(Qop.R00) && all(all(double(Qop.R00)==0))) || all(all(Qop.R00.C==0))
        eq_opts_exc(1) = 1;
    end
end
% Then the map L2[x] --> L2[x]
if Qdim(2)
    % Check if Qop.Rxx{1} is zero.
    if isempty(Qop.Rxx{1}) || (isdouble(Qop.Rxx{1}) && all(all(double(Qop.Rxx{1})==0))) || all(all(Qop.Rxx{1}.C==0))
        eq_opts_exc(2) = 1;
    end
    % Check if Qop.Rxx{2} is zero.
    if isempty(Qop.Rxx{1}) || (isdouble(Qop.Rxx{2}) && all(all(double(Qop.Rxx{2})==0))) || all(all(Qop.Rxx{2}.C==0))
        eq_opts_exc(3) = 1;
    end
    % Check if Qop.Rxx{3} is zero.
    if isempty(Qop.Rxx{1}) || (isdouble(Qop.Rxx{3}) && all(all(double(Qop.Rxx{3})==0))) || all(all(Qop.Rxx{3}.C==0))
        eq_opts_exc(4) = 1;
    elseif ~eq_opts_exc(3)
        % Check if Qop.Rxx{2}=Qop.Rxx{3}.
        Rxx_diff = Qop.Rxx{2}-Qop.Rxx{3};
        if all(all(Rxx_diff.C<=ztol))
            eq_opts_sep(1) = 1;
        end
    end
end
% Then the map L2[y] --> L2[y]
if Qdim(3)
    % Check if Qop.Ryy{1} is zero.
    if isempty(Qop.Ryy{1}) || (isdouble(Qop.Ryy{1}) && all(all(double(Qop.Ryy{1})==0))) || all(all(Qop.Ryy{1}.C==0))
        eq_opts_exc(5) = 1;
    end
    % Check if Qop.Ryy{2} is zero.
    if isempty(Qop.Ryy{2}) || (isdouble(Qop.Ryy{2}) && all(all(double(Qop.Ryy{2})==0))) || all(all(Qop.Ryy{2}.C==0))
        eq_opts_exc(6) = 1;
    end
    % Check if Qop.Ryy{3} is zero.
    if isempty(Qop.Ryy{3}) || (isdouble(Qop.Ryy{3}) && all(all(double(Qop.Ryy{3})==0))) || all(all(Qop.Ryy{3}.C==0))
        eq_opts_exc(7) = 1;
    elseif ~eq_opts_exc(6)
        % Check if Qop.Rxx{2}=Qop.Rxx{3}.
        Ryy_diff = Qop.Ryy{2}-Qop.Ryy{3};
        if all(all(Ryy_diff.C<=ztol))
            eq_opts_sep(2) = 1;
        end
    end
end
% Finally the map L2[x,y]-->L2[x,y];
if Qdim(4)
    % Loop over all parameters, checking if they are zero.
    param_idcs = [8,9,10,11,13,14,12,15,16];
    for ii=1:9
        if isempty(Qop.R22{ii}) || (isdouble(Qop.R22{ii}) && all(all(double(Qop.R22{ii})==0))) || all(all(Qop.R22{ii}.C==0))
            eq_opts_exc(param_idcs(ii)) = 1;
        end
    end
    % Check if Qop.R22{3,1}==Qop.R22{2,1}
    if ~eq_opts_exc(9) && ~eq_opts_exc(10)
        R22_diff = Qop.R22{2,1}-Qop.R22{3,1};
        if all(all(R22_diff.C<=ztol))
            eq_opts_sep(3) = 1;
        end
    end
    % Check if Qop.R22{1,3}==Qop.R22{1,2}
    if ~eq_opts_exc(11) && ~eq_opts_exc(12)
        R22_diff = Qop.R22{1,2}-Qop.R22{1,3};
        if all(all(R22_diff.C<=ztol))
            eq_opts_sep(4) = 1;
        end
    end
    % Check if Qop.R22{3,2}==Qop.R22{2,2} and Qop.R22{3,3}==Qop.R22{2,3}
    if (~eq_opts_exc(13) && ~eq_opts_exc(14)) || (~eq_opts_exc(15) && ~eq_opts_exc(16))
        R22_diff_1 = Qop.R22{2,2}-Qop.R22{3,2};
        R22_diff_2 = Qop.R22{2,3}-Qop.R22{3,3};
        if all(all(R22_diff_1.C<=ztol)) && all(all(R22_diff_2.C<=ztol))
            eq_opts_sep(5) = 1;
        end
    end
    % Check if Qop.R22{2,3}==Qop.R22{2,2} and Qop.R22{3,3}==Qop.R22{3,2}
    if (~eq_opts_exc(13) && ~eq_opts_exc(15)) || (~eq_opts_exc(14) && ~eq_opts_exc(16))
        R22_diff_1 = Qop.R22{2,2}-Qop.R22{2,3};
        R22_diff_2 = Qop.R22{3,2}-Qop.R22{3,3};
        if all(all(R22_diff_1.C<=ztol)) && all(all(R22_diff_2.C<=ztol))
            eq_opts_sep(6) = 1;
        end
    end
end
eq_opts.exclude = eq_opts_exc;
eq_opts.sep = eq_opts_sep;

end