function PIE_cat = vertcat(varargin)
% PIE_CAT = VERTCAT(PIE1,PIE2,...) takes a set of 'pie_struct' objects and
% returns a 'pie_struct' object representing the vertical concatenation
% of the input PIEs. In particular, if PIEi maps nwi exogenous and nui
% controlled inputs through nxi state variables to nzi regulated and nyi
% observed outputs, then [PIE1;PIE2] assumes the inputs of PIE1 and PIE2
% are the same, and maps those nw1=nw2 exogenous and nu1=nu2 controlled
% inputs through (nx1+nx2) state variables to nz1+nz2 regulated and
% ny1+ny2 observed outputs.
%
% INPUT
% - varargin:   list of 'pie_struct' objects which to concatenate.
%
% OUTPUT
% - PIE_cat:    'pie_struct' object representing the PIE defined by the
%               vertical concatenation of the input systems.
%
% NOTES
% For two PIE systems PIE_i defined as
%
% d/dt (Ti*vi(t)+Twi*wi(t)+Tui*ui(t)) = Ai*vi(t) +Bwi*wi(t) +Bui*ui(t)
%                                  zi(t) = Czi*vi(t)+Dzwi*wi(t)+Dzui*ui(t)
%                                  yi(t) = Cyi*vi(t)+Dywi*wi(t)+Dyui*ui(t)
%             Ki*vi(t)+Kwi*wi(t)+Kui*ui(t) = 0
%
% the vertical concatenation PIE = [PIE1;PIE2] takes the form
%
% d/dt (T*v(t)+Tw*w(t)+Tu*u(t)) = A*v(t) +Bw*w(t) +Bu*u(t)
%                          z(t) = Cz*v(t)+Dzw*w(t)+Dzu*u(t)
%                          y(t) = Cy*v(t)+Dyw*w(t)+Dyu*u(t)
%             K*v(t)+Kw*w(t)+Ku*u(t) = 0
%
% where v = [v1;v2], z = [z1;z2], y = [y1;y2], w = w1 = w2, u = u1 = u2,
% and
%
%   T  = [T1,0;0,T2];       Tw  = [Tw1;Tw2];        Tu  = [Tu1;Tu2];
%   A  = [A1,0;0,A2];       Bw  = [Bw1;Bw2];        Bu  = [Bu1;Bu2];
%   Cz = [Cz1,0;0,Cz2];     Dzw = [Dzw1;Dzw2];      Dzu = [Dzu1;Dzu2];
%   Cy = [Cy1,0;0,Cy2];     Dyw = [Dyw1;Dyw2];      Dyu = [Dyu1;Dyu2];
%   K  = [K1,0;0,K2];       Kw  = [Kw1;Kw2];        Ku  = [Ku1;Ku2];
%
% The components of the state and outputs of the PIEs are concatenated in
% the order of the PIEs, whereas the inputs are assumed to be the same. The
% PIEs are initialized if necessary, and the function spaces of the inputs
% of the two PIEs must match. A PIE without algebraic constraint has
% operators K, Kw, and Ku without any rows, so that PIEs with and without
% constraint can be concatenated.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - vertcat
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 12/21/2024: Initial coding;
% DJ, 10/29/2026: Update to new pie_struct structure (copvar operators,
%                   function spaces, and algebraic constraints K, Kw, and
%                   Ku);

% Extract the inputs
if nargin==1
    PIE_cat = varargin{1};
    return
end
PIE1 = varargin{1};
PIE2 = varargin{2};


% Make sure PIEs are properly specified.
try                                                                         % DJ, 10/29/2026
    if ~PIE1.is_initialized
        PIE1 = initialize(PIE1);    
    end
    if ~PIE2.is_initialized
        PIE2 = initialize(PIE2);
    end
catch
    error("One of the PIEs is not properly specified; concatenation not supported.")
end
% Make sure variables and domains match.
if PIE1.dim~=PIE2.dim
    error("PIEs are of different dimensionality; vertcat currently not supported.")
end
if ~isequal(PIE1.space.vars,PIE2.space.vars)
    error("Spatial variables of PIEs must match.")
end
if any(any(PIE1.dom~=PIE2.dom))
    error("Spatial domains of PIEs must match.")
end
% Make sure input spaces match.
if ~isequal(pie_space_canon(PIE1.space.u.depmat,PIE1.space.u.dims),...
            pie_space_canon(PIE2.space.u.depmat,PIE2.space.u.dims))
    error("Dimensions of actuator inputs are not consistent.")
end
if ~isequal(pie_space_canon(PIE1.space.w.depmat,PIE1.space.w.dims),...
            pie_space_canon(PIE2.space.w.depmat,PIE2.space.w.dims))
    error("Dimensions of exogenous inputs are not consistent.")
end


% Construct the concatenated PIE structure
PIE_cat = pie_struct();

% Dynamics
PIE_cat.T = blkdiag(PIE1.T,PIE2.T);
PIE_cat.Tw = [PIE1.Tw;PIE2.Tw];
PIE_cat.Tu = [PIE1.Tu;PIE2.Tu];
PIE_cat.A = blkdiag(PIE1.A,PIE2.A);
PIE_cat.Bw = [PIE1.Bw;PIE2.Bw];
PIE_cat.Bu = [PIE1.Bu;PIE2.Bu];
PIE_cat.Cz = blkdiag(PIE1.Cz,PIE2.Cz);
PIE_cat.Dzw = [PIE1.Dzw;PIE2.Dzw];
PIE_cat.Dzu = [PIE1.Dzu;PIE2.Dzu];
PIE_cat.Cy = blkdiag(PIE1.Cy,PIE2.Cy);
PIE_cat.Dyw = [PIE1.Dyw;PIE2.Dyw];
PIE_cat.Dyu = [PIE1.Dyu;PIE2.Dyu];

% Algebraic constraints on the separate states, for the same inputs.
PIE_cat.K = blkdiag(PIE1.K,PIE2.K);                                         % DJ, 10/29/2026
PIE_cat.Kw = [PIE1.Kw;PIE2.Kw];
PIE_cat.Ku = [PIE1.Ku;PIE2.Ku];

% Concatenate the components of the state and outputs; the inputs remain
% the same.
for type = {'x','z','y'}                                                    % DJ, 10/29/2026
    PIE_cat.space.(type{1}) = struct(...
        'dims',[PIE1.space.(type{1}).dims; PIE2.space.(type{1}).dims],...
        'depmat',[PIE1.space.(type{1}).depmat; PIE2.space.(type{1}).depmat]);
end
PIE_cat.space.w = PIE1.space.w;
PIE_cat.space.u = PIE1.space.u;

% Perform remaining concatenations.
if nargin>2
    PIE_cat = vertcat(PIE_cat,varargin{3:end});
end

end
