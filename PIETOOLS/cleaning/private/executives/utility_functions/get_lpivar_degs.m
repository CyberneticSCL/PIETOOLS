function Qdeg = get_lpivar_degs(Rop,Top)
% QDEG = GET_LPIVAR_DEGS(ROP,TOP) takes an opvar or dopvar object ROP, and
% returns an array QDEG of degrees to pass to "lpivar" for generating an
% operator variable QOP such that TOP'*QOP can match all monomials
% appearing in ROP, i.e. we can enforce TOP'*QOP==ROP.
%
% INPUT
% - Rop:    'opvar', 'dopvar', 'opvar2d' or 'dopvar2d' object;
% - Top:    'opvar', 'dopvar', 'opvar2d' or 'dopvar2d' object;
%
% OUTPUT
% - Qdeg:   1x3 array of nonnegative integers if input is 1D ('opvar'),
%           or struct with fields 'dx', 'dy', and 'd2' if input is 2D
%           specifying maximal degrees of monomials in call to "lpivar". 
%           In particular, calling Qop = lpivar(Top.dim,Qdeg), for 1D case,
%           Qdeg(1) will determine the degrees of Qop.Q1, Qop.Q2, and 
%           Qop.R.R0 in s, Qdeg(2) the degrees of s_dum in Qop.R.R1 and 
%           Qop.R.R2, and Qdeg(3) the degrees of s in Qop.R.R2.
%           For 2D case, Qdeg.dx{1} will determine degrees of Qop.Rxx{1} in
%           s1, Qdeg.dx{2} degrees of Qop.Rxx{2} in s1 and s1_dum,
%           Qdeg.dx{3} degrees of Qop.Rxx{3} in s1 and s1_dum, Qdeg.dy
%           determines degrees of Qop.Ryy, Qdeg.d2{i,j} determines degrees
%           of Qop.R22{i,j};
%
% NOTES
% The current implementation is very safe but likely also unnecessarily
% expensive. It just checks the maximal degrees of all monomials defining
% Rop, not at all accounting for the operator Top, or distinguishing 
% degrees in s and s_dum. Future updates may improve upon this.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getdegs_lpivar
%
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 01/06/2025: Initial coding;
% DJ, 01/17/2025: Update to support 2D operators;
%

if (isa(Rop,'opvar') || isa(Rop,'dopvar')) && (isa(Top,'opvar') || isa(Top,'dopvar'))
    % % 1D case
    % Check maximal degrees of monomials in Rop.Q1, Rop.Q2, and Rop.R.R0.
    Rop_deg1 = 0;
    if ~isempty(Rop.Q1) && (isa(Rop.Q1,'polynomial')||isa(Rop.Q1,'dpvar')) && ~isempty(Rop.Q1.degmat)
        Rop_deg1 = max(Rop_deg1,max(Rop.Q1.degmat));
    end
    if ~isempty(Rop.Q2) && (isa(Rop.Q2,'polynomial')||isa(Rop.Q2,'dpvar')) && ~isempty(Rop.Q2.degmat)
        Rop_deg1 = max(Rop_deg1,max(Rop.Q2.degmat));
    end
    if ~isempty(Rop.R.R0) && (isa(Rop.R.R0,'polynomial')||isa(Rop.R.R0,'dpvar')) && ~isempty(Rop.R.R0.degmat)
        Rop_deg1 = max(Rop_deg1,max(Rop.R.R0.degmat));
    end
    % Check maximal degrees of monomials in Rop.R.R1, Rop.R.R2.
    Rop_deg2 = 0;
    if ~isempty(Rop.R.R1) && (isa(Rop.R.R1,'polynomial')||isa(Rop.R.R1,'dpvar')) && ~isempty(Rop.R.R1.degmat)
        Rop_deg2 = max(Rop_deg2,max(max(Rop.R.R1.degmat)));
    end
    if ~isempty(Rop.R.R2) && (isa(Rop.R.R2,'polynomial')||isa(Rop.R.R2,'dpvar')) && ~isempty(Rop.R.R2.degmat)
        Rop_deg2 = max(Rop_deg2,max(max(Rop.R.R2.degmat)));
    end
    Qdeg = [Rop_deg1,Rop_deg2,Rop_deg2-1];  % subtract 1 to account for multiplication with Top

elseif (isa(Rop,'opvar2d') || isa(Rop,'dopvar2d')) && (isa(Top,'opvar2d') || isa(Top,'dopvar2d'))
    % % 2D case
    [~,Rmaxdegs,~] = getdeg(Rop);
    Qdeg = struct();
    Qdeg.dx = cell(3,1);
    Qdeg.dy = cell(1,3);
    Qdeg.d2 = cell(3,3);

    % Degrees in just s1 and s1_dum
    Qdeg.dx{1} = max([Rmaxdegs.Rx0(2,1),Rmaxdegs.R0x(2,1),Rmaxdegs.Rxx{1}(2,1)]);
    Qdeg.dx{2} = Rmaxdegs.Rxx{2}([2;3;4],1);
    Qdeg.dx{3} = Rmaxdegs.Rxx{3}([2;3;4],1);

    % Degrees in just s2 and s2_dum
    Qdeg.dy{1} = max([Rmaxdegs.Ry0(1,2),Rmaxdegs.R0y(1,2),Rmaxdegs.Ryy{1}(1,2)]);
    Qdeg.dy{2} = Rmaxdegs.Ryy{2}(1,[2,3,4]);
    Qdeg.dy{3} = Rmaxdegs.Ryy{3}(1,[2,3,4]);

    % Degrees in (s1,s2) and (s1_dum,s2_dum);
    Qdeg.d2{1,1} = max(cat(3,Rmaxdegs.R22{1,1}(1:2,1:2),Rmaxdegs.Rxy(1:2,1:2),Rmaxdegs.Ryx(1:2,1:2),...
                            Rmaxdegs.R02(1:2,1:2),Rmaxdegs.R20(1:2,1:2),...
                            Rmaxdegs.Rx2{1}(1:2,1:2),Rmaxdegs.Ry2{1}(1:2,1:2), ...
                            Rmaxdegs.R2x{1}(1:2,1:2),Rmaxdegs.R2y{1}(1:2,1:2)),[],3);
    Qdeg.d2{2,1} = max(cat(3,Rmaxdegs.R22{2,1}(1:4,1:2),Rmaxdegs.Rx2{2}(1:4,1:2),Rmaxdegs.R2x{2}(1:4,1:2)),[],3);
    Qdeg.d2{3,1} = max(cat(3,Rmaxdegs.R22{3,1}(1:4,1:2),Rmaxdegs.Rx2{3}(1:4,1:2),Rmaxdegs.R2x{3}(1:4,1:2)),[],3);
    Qdeg.d2{1,2} = max(cat(3,Rmaxdegs.R22{1,2}(1:2,1:4),Rmaxdegs.Ry2{2}(1:2,1:4),Rmaxdegs.R2y{2}(1:2,1:4)),[],3);
    Qdeg.d2{1,3} = max(cat(3,Rmaxdegs.R22{1,3}(1:2,1:4),Rmaxdegs.Ry2{3}(1:2,1:4),Rmaxdegs.R2y{3}(1:2,1:4)),[],3);
    Qdeg.d2{2,2} = Rmaxdegs.R22{2,2};
    Qdeg.d2{3,2} = Rmaxdegs.R22{3,2};
    Qdeg.d2{2,3} = Rmaxdegs.R22{2,3};
    Qdeg.d2{3,3} = Rmaxdegs.R22{3,3};
else
    error("Inputs should be of type 'opvar', 'dopvar', 'opvar2d', or 'dopvar2d'.")
end

end