function Pop = sopvar2opvar2d(Psop)
% POP = SOPVAR2OPVAR2D(PSOP) takes a sopvar object representing a 2D PI
% operator component and returns an opvar2d object representing the same
% operator.
%
% INPUTS
% - Psop:   'sopvar' object representing a 2D PI operator. It can map
%           between different function space, but the input and output
%           domains cannot be coupled domains of different dimensionality
%           (i.e. L2^n and R^m are supported, but R^m x L2^n is not).
%
% OUTPUTS
% - obj:    'opvar2d' object representing the same operator as the input;
%
% NOTES
% At this time, the output operator is always expressed in terms of
% variables [s1;s2] (primary) and [s1_dum;s2_dum] (dummy), independent of 
% the variables specified in Psop.vars
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sopvar2opvar2d
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
% DJ, 05/27/2026: Initial coding

% Extract the dimension, variables and domain of the operator
dims = Psop.dims;
vars = Psop.vars;
dom = Psop.dom;
Pparams = Psop.params;

% Initialize an empty operator
Pop = opvar2d();
Pop.var1 = polynomial({'s1';'s2'});  % use default variables
Pop.var2 = polynomial({'s1_dum';'s2_dum'});

% Determine between what function spaces the operator maps
if numel(unique([vars.in,vars.out]))>2
    error("Operator maps between functions of more than two distinct variables.")
end
mapsx = numel(vars.in)>=1;
maps2x = numel(vars.out)>=1;
mapsy = numel(vars.in)>=2;
maps2y = numel(vars.out)>=2;
% NOTE: sopvar does not inherently distinguish between x and y variable
% --> we perform this distinction only in special case of default vars
if mapsx && ~mapsy && (strcmp(vars.in{1},'s2') || strcmp(vars.in{1},'y'))
    mapsx = false;
    mapsy = true;
end
if maps2x && ~maps2y && (strcmp(vars.out{1},'s2') || strcmp(vars.out{1},'y'))
    maps2x = false;
    maps2y = true;
end
if mapsx && maps2x && ~strcmp(vars.in{1},vars.out{1})
    % If input and output variables don't match, the space cannot be the
    % same
    maps2x = false;
    maps2y = true;
end
x_idx_in = mapsx;       
x_idx_out = maps2x;
y_idx_in = mapsx+mapsy;
y_idx_out = maps2x+maps2y;

% Set the domain of the variables
if (mapsx && maps2x && ~all(dom.in(x_idx_in,:)==dom.out(x_idx_out,:))) ||...
        (mapsy && maps2y && ~all(dom.in(y_idx_in,:)==dom.out(y_idx_out,:)))
    % Input and output domains should match for same variable
    error("Input and output domains of the same variable should match.")
end
if mapsx
    Pop.dom(1,:) = dom.in(x_idx_in,:);
elseif maps2x
    Pop.dom(1,:) = dom.out(x_idx_out,:);
end
if mapsy
    Pop.dom(2,:) = dom.in(y_idx_in,:);
elseif maps2y
    Pop.dom(2,:) = dom.out(y_idx_out,:);
end

% Determine what parameter in the opvar2d structure is non-empty
cidx = [mapsx,mapsy]*[1;2]+1;
ridx = [maps2x,maps2y]*[1;2]+1;
Pdim = zeros(4,2);
Pdim(ridx,1) = dims(1);     Pdim(cidx,2) = dims(2);
Pop.dim = Pdim;
Rparam_names = {'R00','R0x','R0y','R02';
                'Rx0','Rxx','Rxy','Rx2';
                'Ry0','Ryx','Ryy','Ry2';
                'R20','R2x','R2y','R22'};
Rname = Rparam_names{ridx,cidx};

% Determine the primary and dummy variable names used in the parameters
var1 = {};        var2 = {};
if maps2x
    var1 = [var1,{'s1'}];
end
if maps2y
    var1 = [var1,{'s2'}];
end
if mapsx
    if ~maps2x
        var2 = [var2,{'s1'}];
    else
        var2 = [var2,{'s1_dum'}];
    end
end
if mapsy
    if ~maps2y
        var2 = [var2,{'s2'}];
    else
        var2 = [var2,{'s2_dum'}];
    end
end

% Convert the parameters to 'polynomial' class objects
Rcell = cell(size(Pparams));
for i=1:numel(Pparams)
    [i1,i2] = ind2sub(size(Pparams),i);
    var2_i = var2;
    if i1==1
        % Multiplier parameter does not have dummy variables
        if mapsx && maps2x
            var2_i(x_idx_in) = var1(x_idx_out);
        elseif mapsy && maps2y
            var2_i(y_idx_in) = var1(y_idx_out);
        end
    end
    if i2==1 && (mapsx && maps2x)
        % In 2D case, i2==1 means multiplier along y-direction
        if mapsy && maps2y
            var2_i(y_idx_in) = var1(y_idx_out);
        end
    end
    Ri = quadPoly(Pparams{i}, Psop.ZL, Psop.ZR, dims, var1, var2_i, 0);
    %Ri = combine(Ri);
    Rcell{i} = combine(quadPoly.quadPoly2polynomial(Ri));
end

% Finally, set the parameters
if ~isa(Pop.(Rname),'cell')
    Pop.(Rname) = Rcell{1};
else
    Pop.(Rname) = reshape(Rcell,size(Pop.(Rname)));
end

end