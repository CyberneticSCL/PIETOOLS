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
% The primary variables are taken from Psop.vars. The dummies are still     % MMP, 08/30/2026
% <primary>_dum by convention: an 'sopvar' records only one name per        % MMP, 08/30/2026
% direction, so the input-side name is not stored and cannot be recovered.  % MMP, 08/30/2026
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
% MMP, 08/30/2026: Take the primary variables from Psop.vars rather than
%                  hardcoding s1,s2, which relabelled any operator not
%                  already in those. Dummy convention unchanged.

% Extract the dimension, variables and domain of the operator
dims = Psop.dims;
vars = Psop.vars;
dom = Psop.dom;
Pparams = Psop.params;

% Initialize an empty operator
Pop = opvar2d();
% Pop.var1/var2 are set below, once the variable names are read from Psop.  % MMP, 08/30/2026

% Determine between what function spaces the operator maps
if numel(unique([vars.in,vars.out]))>2
    error("Operator maps between functions of more than two distinct variables.")
end
% An 'sopvar' does not itself distinguish the x and y directions, so the    % MMP, 09/09/2026
% two distinct variables are assigned to them in sorted order, matching     % MMP, 09/09/2026
% the default {'s1','s2'}; a lone variable named 's2' or 'y' still          % MMP, 09/09/2026
% occupies the y direction. Each direction is then LOCATED IN EACH LIST     % MMP, 09/09/2026
% BY NAME rather than assumed to sit at position 1 and 2, because the       % MMP, 09/09/2026
% constructor stores vars.out as [S2,S3] and vars.in as [S3,S1]: for a      % MMP, 09/09/2026
% block L_2[x] -> L_2[x,y] that puts y first in vars.out, and reading       % MMP, 09/09/2026
% position 1 as x then took y for x and reported the two as unmatched.      % MMP, 09/09/2026
allv = unique([vars.in(:).',vars.out(:).']);                                % MMP, 09/09/2026
if isempty(allv)                                                            % MMP, 09/09/2026
    % R^n -> R^m: neither direction is involved.                            % MMP, 09/09/2026
    xy = {'',''};                                                           % MMP, 09/09/2026
elseif isscalar(allv) && (strcmp(allv{1},'s2') || strcmp(allv{1},'y'))      % MMP, 09/09/2026
    xy = {'',allv{1}};                                                      % MMP, 09/09/2026
elseif isscalar(allv)                                                       % MMP, 09/09/2026
    xy = {allv{1},''};                                                      % MMP, 09/09/2026
else                                                                        % MMP, 09/09/2026
    xy = {allv{1},allv{2}};                                                 % MMP, 09/09/2026
end                                                                         % MMP, 09/09/2026
[mapsx, x_idx_in ] = name_at(vars.in ,xy{1});                               % MMP, 09/09/2026
[maps2x,x_idx_out] = name_at(vars.out,xy{1});                               % MMP, 09/09/2026
[mapsy, y_idx_in ] = name_at(vars.in ,xy{2});                               % MMP, 09/09/2026
[maps2y,y_idx_out] = name_at(vars.out,xy{2});                               % MMP, 09/09/2026

% Name each direction from whichever role the operator involves; the        % MMP, 08/30/2026
% direction it does not involve keeps the default.                          % MMP, 08/30/2026
sv = {'s1','s2'};                                                           % MMP, 08/30/2026
if mapsx,   sv{1} = char(vars.in{x_idx_in});    end                         % MMP, 08/30/2026
if maps2x,  sv{1} = char(vars.out{x_idx_out});  end                         % MMP, 08/30/2026
if mapsy,   sv{2} = char(vars.in{y_idx_in});    end                         % MMP, 08/30/2026
if maps2y,  sv{2} = char(vars.out{y_idx_out});  end                         % MMP, 08/30/2026
Pop.var1 = polynomial(sv(:));                                               % MMP, 08/30/2026
Pop.var2 = polynomial(strcat(sv(:),'_dum'));                                % MMP, 08/30/2026

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
% These are handed to 'quadPoly' below alongside Psop.ZL and Psop.ZR, so     % MMP, 09/09/2026
% they must be in the STORED order of vars.out and vars.in, not in (x,y)     % MMP, 09/09/2026
% order: ZL{k} belongs to vars.out(k). A variable that appears only on the   % MMP, 09/09/2026
% input side is integrated out and keeps its own name; one that appears on   % MMP, 09/09/2026
% both sides is the dummy copy.                                             % MMP, 09/09/2026
var1 = vars.out(:).';                                                       % MMP, 09/09/2026
var2 = cell(1,numel(vars.in));                                              % MMP, 09/09/2026
for k = 1:numel(vars.in)                                                    % MMP, 09/09/2026
    if ismember(vars.in{k},vars.out)                                        % MMP, 09/09/2026
        var2{k} = [char(vars.in{k}),'_dum'];                                % MMP, 09/09/2026
    else                                                                    % MMP, 09/09/2026
        var2{k} = char(vars.in{k});                                         % MMP, 09/09/2026
    end                                                                     % MMP, 09/09/2026
end                                                                         % MMP, 09/09/2026

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

%%
function [tf,idx] = name_at(list,nm)                                        % MMP, 09/09/2026
% Position of variable 'nm' in 'list', or (false,0) when absent or unnamed. % MMP, 09/09/2026
if isempty(nm)                                                              % MMP, 09/09/2026
    tf = false;     idx = 0;                                                % MMP, 09/09/2026
else                                                                        % MMP, 09/09/2026
    [tf,idx] = ismember(nm,list);                                           % MMP, 09/09/2026
end                                                                         % MMP, 09/09/2026

end                                                                         % MMP, 09/09/2026
