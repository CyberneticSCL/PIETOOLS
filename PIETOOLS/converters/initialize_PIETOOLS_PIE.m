function PIE = initialize_PIETOOLS_PIE(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = initialize_PIETOOLS_PIE(PIE) takes a PIE data structure and
% checks that all the necessary fields are appropriately specified, and
% assigns a default value to all the optional fields that have not been
% specified.
% 
% INPUT
%   PIE: "struct" or "pie_struct" class object.
%
% OUTPUT
%   PIE: "pie_struct" class object describing the same system as the input.
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Initial coding SS - 08/01/2022

if nargin==0
    PIE = pie_struct();
    return
end
if isa(PIE,'pie_struct')
    PIE = initialize(PIE);
    return
elseif ~isa(PIE,'struct')
    error('Input must be a ''struct'' or ''pie_struct'' class object.')
end

if ~isfield(PIE,'dom')
    PIE.dom = [0,1];
end
if ~isfield(PIE,'vars')
    PIE.vars = [pvar('s'),pvar('theta')];
end

% store dimension information in the structure
dimstruct = struct('nx',0,'nw',0,'nu',0,'nz',0,'ny',0);

% mandatory fields
if ~isfield(PIE,'T')
    error("PIE.T is a mandatory element of PIE structure and cannot be defaulted to zero");
else
    if any(PIE.T.dim(:,1)~=PIE.T.dim(:,2))
        error("PIE.T operator must be a map between symmetric spaces");
    end
    dimstruct.nx = PIE.T.dim(:,1);
end
if ~isfield(PIE,'A')
    error("PIE.A is a mandatory element of PIE structure and cannot be defaulted to zero");
else
    if any(PIE.A.dim(:,1)~=PIE.A.dim(:,2))||any(PIE.A.dim(:,2)~=PIE.T.dim(:,2))
        error("PIE.A operator must be a map between symmetric spaces");
    end
end

% find remaining dimensions
dimstruct.nw = find_nw(PIE);
dimstruct.nu = find_nu(PIE);
dimstruct.nz = find_nz(PIE);
dimstruct.ny = find_ny(PIE);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% optional fields, separate based on 1d and 2d cases
if ~isfield(PIE,'C1')
    if PIE.dim==1
        opvar tmp; tmp.dim = [[dimstruct.nz;0],PIE.T.dim(:,1)];
        PIE.C1 = tmp;
    else % 2d
        opvar2d tmp; tmp.dim = [[dimstruct.nz;0;0;0],PIE.T.dim(:,1)];
        PIE.C1 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.C1.dim(:)~=[[dimstruct.nz;0];PIE.T.dim(:,1)])
            warning('PIE.C1, PIE.T have incompatible dimensions. Defaulting C1 to zero operator');
            opvar tmp; tmp.dim = [[dimstruct.nz;0],PIE.T.dim(:,1)];
            PIE.C1 = tmp;
        end
    else
        if any(PIE.C1.dim(:)~=[[dimstruct.nz;0;0;0];PIE.T.dim(:,1)])
            warning('PIE.C1, PIE.T have incompatible dimensions. Defaulting C1 to zero operator');
            opvar2d tmp; tmp.dim = [[dimstruct.nz;0;0;0],PIE.T.dim(:,1)];
            PIE.C1 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'C2')
    if PIE.dim==1
        opvar tmp; tmp.dim = [[dimstruct.ny;0],PIE.T.dim(:,1)];
        PIE.C2 = tmp;
    else % 2d
        opvar2d tmp; tmp.dim = [[dimstruct.ny;0;0;0],PIE.T.dim(:,1)];
        PIE.C2 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.C2.dim(:)~=[[dimstruct.ny;0];PIE.T.dim(:,1)])
            warning('PIE.C2, PIE.T have incompatible dimensions. Defaulting C2 to zero operator');
            opvar tmp; tmp.dim = [[dimstruct.ny;0],PIE.T.dim(:,1)];
            PIE.C2 = tmp;
        end
    else
        if any(PIE.C2.dim(:)~=[[dimstruct.ny;0;0;0];PIE.T.dim(:,1)])
            warning('PIE.C2, PIE.T have incompatible dimensions. Defaulting C2 to zero operator');
            opvar2d tmp; tmp.dim = [[dimstruct.ny;0;0;0],PIE.T.dim(:,1)];
            PIE.C2 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'Tw')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0]];
        PIE.Tw = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0;0;0]];
        PIE.Tw = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.T.dim(:,1)~=PIE.Tw.dim(:,1))||any(PIE.Tw.dim(:,2)~=[dimstruct.nw;0])
            warning('PIE.Tw and PIE.T map to different spaces. Defaulting Tw to zero operator');
            opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0]];
            PIE.Tw = tmp;
        end
    else
        if any(PIE.T.dim(:,1)~=PIE.Tw.dim(:,1))||any(PIE.Tw.dim(:,2)~=[dimstruct.nw;0;0;0])
            warning('PIE.Tw and PIE.T map to different spaces. Defaulting Tw to zero operator');
            opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0;0;0]];
            PIE.Tw = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'Tu')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0]];
        PIE.Tu = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0;0;0]];
        PIE.Tu = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.T.dim(:,1)~=PIE.Tu.dim(:,1))||any(PIE.Tu.dim(:,2)~=[dimstruct.nu;0])
            warning('PIE.Tu and PIE.T map to different spaces. Defaulting Tu to zero operator');
            opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0]];
            PIE.Tu = tmp;
        end
    else
        if any(PIE.T.dim(:,1)~=PIE.Tu.dim(:,1))||any(PIE.Tu.dim(:,2)~=[dimstruct.nu;0;0;0])
            warning('PIE.Tu and PIE.T map to different spaces. Defaulting Tu to zero operator');
            opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0;0;0]];
            PIE.Tu = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'B1')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0]];
        PIE.B1 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0;0;0]];
        PIE.B1 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.B1.dim(:)~=PIE.Tw.dim(:))
            warning('PIE.Tw, PIE.B1 map to different spaces. Defaulting B1 to zero operator');
            opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0]];
            PIE.B1 = tmp;
        end
    else
        if any(PIE.B1.dim(:)~=PIE.Tw.dim(:))
            warning('PIE.Tw, PIE.B1 map to different spaces. Defaulting B1 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nw;0;0;0]];
            PIE.B1 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'B2')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0]];
        PIE.B2 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0;0;0]];
        PIE.B2 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.B2.dim(:)~=PIE.Tu.dim(:))
            warning('PIE.Tu, PIE.B2 map to different spaces. Defaulting B2 to zero operator');
            opvar tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0]];
            PIE.B2 = tmp;
        end
    else
        if any(PIE.B2.dim(:)~=PIE.Tu.dim(:))
            warning('PIE.Tu, PIE.B2 map to different spaces. Defaulting B2 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.T.dim(:,1),[dimstruct.nu;0;0;0]];
            PIE.B2 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'D11')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nw;0]];
        PIE.D11 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nw;0;0;0]];
        PIE.D11 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.D11.dim(:)~=[PIE.C1.dim(:,1);[dimstruct.nw;0]])
            warning('PIE.C1, PIE.D11 have incompatible dimensions. Defaulting D11 to zero operator');
            opvar tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nw;0]];
            PIE.D11 = tmp;
        end
    else
        if any(PIE.D11.dim(:)~=[PIE.C1.dim(:,1);[dimstruct.nw;0;0;0]])
            warning('PIE.C1, PIE.D11 have incompatible dimensions. Defaulting D11 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nw;0;0;0]];
            PIE.D11 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'D12')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nu;0]];
        PIE.D12 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nu;0;0;0]];
        PIE.D12 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.D12.dim(:)~=[PIE.C1.dim(:,1);[dimstruct.nu;0]])
            warning('PIE.C1, PIE.D12 have incompatible dimensions. Defaulting D12 to zero operator');
            opvar tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nu;0]];
            PIE.D12 = tmp;
        end
    else
        if any(PIE.D12.dim(:)~=[PIE.C1.dim(:,1);[dimstruct.nu;0;0;0]])
            warning('PIE.C1, PIE.D12 have incompatible dimensions. Defaulting D12 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.C1.dim(:,1),[dimstruct.nu;0;0;0]];
            PIE.D12 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'D21')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nw;0]];
        PIE.D21 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nw;0;0;0]];
        PIE.D21 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.D21.dim(:)~=[PIE.C2.dim(:,1);[dimstruct.nw;0]])
            warning('PIE.C2, PIE.D21 have incompatible dimensions. Defaulting D21 to zero operator');
            opvar tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nw;0]];
            PIE.D21 = tmp;
        end
    else
        if any(PIE.D21.dim(:)~=[PIE.C2.dim(:,1);[dimstruct.nw;0;0;0]])
            warning('PIE.C2, PIE.D21 have incompatible dimensions. Defaulting D21 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nw;0;0;0]];
            PIE.D21 = tmp;
        end
    end
end
%--------------------------------------------------------------------------
if ~isfield(PIE,'D22')
    if PIE.dim==1
        opvar tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nu;0]];
        PIE.D22 = tmp;
    else
        opvar2d tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nu;0;0;0]];
        PIE.D22 = tmp;
    end
else
    if PIE.dim==1
        if any(PIE.D22.dim(:)~=[PIE.C2.dim(:,1);[dimstruct.nu;0]])
            warning('PIE.C2, PIE.D22 have incompatible dimensions. Defaulting D22 to zero operator');
            opvar tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nu;0]];
            PIE.D22 = tmp;
        end
    else
        if any(PIE.D22.dim(:)~=[PIE.C2.dim(:,1);[dimstruct.nu;0;0;0]])
            warning('PIE.C2, PIE.D22 have incompatible dimensions. Defaulting D22 to zero operator');
            opvar2d tmp; tmp.dim = [PIE.C2.dim(:,1),[dimstruct.nu;0;0;0]];
            PIE.D22 = tmp;
        end
    end
end
end


function val = find_nw(PIE)
flag =1;
if isfield(PIE,'Tw')
    val = PIE.Tw.dim(1,2);
    flag = 0;
end
if isfield(PIE,'B1')
    if ~flag && (val~=PIE.B1.dim(1,2))
        flag=1;
    else
        val = PIE.B1.dim(1,2);
        flag = 0;
    end
end
if isfield(PIE,'D11')
    if ~flag && (val~=PIE.D11.dim(1,2))
        flag=1;
    else
        val = PIE.D11.dim(1,2);
        flag = 0;
    end
end
if isfield(PIE,'D21')
    if ~flag && (val~=PIE.D21.dim(1,2))
        flag=1;
    else
        val = PIE.D21.dim(1,2);
        flag = 0;
    end
end
if flag
    val =0;
end
end
function val = find_nu(PIE)
flag =1;
if isfield(PIE,'Tu')
    val = PIE.Tu.dim(1,2);
    flag = 0;
end
if isfield(PIE,'B2')
    if ~flag && (val~=PIE.B2.dim(1,2))
        flag=1;
    else
        val = PIE.B2.dim(1,2);
        flag = 0;
    end
end
if isfield(PIE,'D12')
    if ~flag && (val~=PIE.D12.dim(1,2))
        flag=1;
    else
        val = PIE.D12.dim(1,2);
        flag = 0;
    end
end
if isfield(PIE,'D22')
    if ~flag && (val~=PIE.D22.dim(1,2))
        flag=1;
    else
        val = PIE.D22.dim(1,2);
        flag = 0;
    end
end
if flag
    val =0;
end
end
function val = find_nz(PIE)
flag =1;
if isfield(PIE,'C1')
    val = PIE.C1.dim(1,1);
    flag = 0;
end
if isfield(PIE,'D11')
    if ~flag && (val~=PIE.D11.dim(1,1))
        flag=1;
    else
        val = PIE.D11.dim(1,1);
        flag = 0;
    end
end
if isfield(PIE,'D12')
    if ~flag && (val~=PIE.D12.dim(1,1))
        flag=1;
    else
        val = PIE.D12.dim(1,1);
        flag = 0;
    end
end
if flag
    val =0;
end
end
function val = find_ny(PIE)
flag =1;
if isfield(PIE,'C2')
    val = PIE.C2.dim(1,1);
    flag = 0;
end
if isfield(PIE,'D21')
    if ~flag && (val~=PIE.D21.dim(1,1))
        flag=1;
    else
        val = PIE.D21.dim(1,1);
        flag = 0;
    end
end
if isfield(PIE,'D22')
    if ~flag && (val~=PIE.D22.dim(1,1))
        flag=1;
    else
        val = PIE.D22.dim(1,1);
        flag = 0;
    end
end
if flag
    val =0;
end
end