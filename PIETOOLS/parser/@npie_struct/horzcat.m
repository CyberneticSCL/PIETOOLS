function PIE_cat = horzcat(varargin)
% PIE_CAT = HORZCAT(PIE1,PIE2,...) takes a set of 'npie_struct' objects and
% returns an 'npie_struct' object representing the horizontal concatenation
% of the input PIEs. In particular, if PIEi maps nwi exogenous and nui
% controlled inputs through nxi state variables to nzi regulated and nyi
% observed outputs, then [PIE1,PIE2] assumes the outputs of PIE1 and PIE2
% are the same, and maps (nw1+nw2) exogenous and (nu1+nu2) controlled
% inputs through (nx1+nx2) state variables to those same nz1=nz2 regulated
% and ny1=ny2 observed outputs.
%
% INPUT
% - varargin:   list of 'npie_struct' objects which to concatenate.
%
% OUTPUT
% - PIE_cat:    'npie_struct' object representing the PIE defined by the
%               horizontal concatenation of the input systems.
%
% NOTES
% For two PIE systems PIE_i defined as
%
%  d/dt (Ti*vi(t)+Twi*wi(t)+Tui*ui(t)) = fi(vi(t))  +Bwi*wi(t) +Bui*ui(t)
%                                 z(t) = gzi(vi(t))+Dzwi*wi(t)+Dzui*ui(t)
%                                 y(t) = gyi(vi(t))+Dywi*wi(t)+Dyui*ui(t)
%         Ki*vi(t)+Kwi*wi(t)+Kui*ui(t) = 0
%
% the horizontal concatenation PIE = [PIE1,PIE2] takes the form
%
%        d/dt (T*v(t)+Tw*w(t)+Tu*u(t)) = f(v(t))  +Bw*w(t) +Bu*u(t)
%                                 z(t) = gz(v(t))+Dzw*w(t)+Dzu*u(t)
%                                 y(t) = gy(v(t))+Dyw*w(t)+Dyu*u(t)
%              K*v(t)+Kw*w(t)+Ku*u(t) = 0
%
% where v = [v1;v2], w = [w1;w2], u = [u1;u2], and
%
%   T  = [T1,0;0,T2];       Tw  = [Tw1,0;0,Tw2];    Tu  = [Tu1,0;0,Tu2];
%   f  = [f1;f2];           Bw  = [Bw1,0;0,Bw2];    Bu  = [Bu1,0;0,Bu2];
%   gz = gz1+gz2;           Dzw = [Dzw1,Dzw2];      Dzu = [Dzu1,Dzu2];
%   gy = gy1+gy2;           Dyw = [Dyw1,Dyw2];      Dyu = [Dyu1,Dyu2];
%   K  = [K1,0;0,K2];       Kw  = [Kw1,0;0,Kw2];    Ku  = [Ku1,0;0,Ku2];
%
% where f = [f1;f2] stacks the polynomials f1(v1) and f2(v2), and gz =
% gz1+gz2 adds the polynomials gz1(v1) and gz2(v2), in terms of the
% distinct state variables of v1 and v2.
%
% The components of the state and inputs of the PIEs are concatenated in
% the order of the PIEs, whereas the outputs are assumed to be the same.
% The PIEs are initialized if necessary, and the function spaces of the
% outputs of the two PIEs must match. A PIE without algebraic constraint
% has operators K, Kw, and Ku without any rows, so that PIEs with and
% without constraint can be concatenated. The names of the state components
% are kept if they are unique and sorted for the concatenated state (see
% `initialize`). Otherwise, they are replaced by the default names, and the
% state variables of the polynomials f, gz, and gy are renamed accordingly.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - horzcat
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
% DJ, 10/29/2026: Initial coding, based on the horzcat function of
%                   pie_struct;

% Extract the inputs
if nargin==1
    PIE_cat = varargin{1};
    return
end
PIE1 = varargin{1};
PIE2 = varargin{2};


% Check if concatenation is supported
% Make sure PIEs are properly specified.
try 
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
    error("PIEs are of different dimensionality; horzcat currently not supported.")
end
if ~isequal(PIE1.space.vars,PIE2.space.vars)
    error("Spatial variables of PIEs must match.")
end
if any(any(PIE1.dom~=PIE2.dom))
    error("Spatial domains of PIEs must match.")
end
% Make sure output spaces match.
if ~isequal(pie_space_canon(PIE1.space.y.depmat,PIE1.space.y.dims),...
            pie_space_canon(PIE2.space.y.depmat,PIE2.space.y.dims))
    error("Dimensions of observed outputs are not consistent.")
end
if ~isequal(pie_space_canon(PIE1.space.z.depmat,PIE1.space.z.dims),...
            pie_space_canon(PIE2.space.z.depmat,PIE2.space.z.dims))
    error("Dimensions of regulated outputs are not consistent.")
end


% Construct the concatenated PIE structure
PIE_cat = npie_struct();

% Keep the names of the state components if they are unique and sorted, as
% required by `initialize`. Otherwise, use the default names of `initialize`.
n1 = numel(PIE1.space.x.dims);      n2 = numel(PIE2.space.x.dims);
names_x = [PIE1.space.x.varname; PIE2.space.x.varname];
if numel(unique(names_x))<n1+n2 || ~issorted(names_x)
    fmt = ['%0',num2str(numel(num2str(n1+n2))),'d'];
    names_x = arrayfun(@(i) ['x',num2str(i,fmt),'_pie'],(1:n1+n2)','UniformOutput',false);
end
% Express the polynomials in terms of the names of the concatenated state.
poly = struct();
for name = {'f','gz','gy'}
    poly1 = subs(PIE1.(name{1}),PIE1.space.x.varname,names_x(1:n1));
    poly2 = subs(PIE2.(name{1}),PIE2.space.x.varname,names_x(n1+1:end));
    if strcmp(name{1},'f')
        poly.f = [poly1;poly2];         % states stacked, as the rows of T
    elseif size(poly1,1)>0
        poly.(name{1}) = poly1 + poly2;     % outputs of the PIEs are the same
    else
        poly.(name{1}) = [];            % no outputs, initialize adds zero polynomial
    end
end

% Concatenate the dynamics
PIE_cat.T = blkdiag(PIE1.T,PIE2.T);
PIE_cat.Tw = blkdiag(PIE1.Tw,PIE2.Tw);
PIE_cat.Tu = blkdiag(PIE1.Tu,PIE2.Tu);
PIE_cat.f = poly.f;
PIE_cat.Bw = blkdiag(PIE1.Bw,PIE2.Bw);
PIE_cat.Bu = blkdiag(PIE1.Bu,PIE2.Bu);
PIE_cat.gz = poly.gz;
PIE_cat.Dzw = [PIE1.Dzw,PIE2.Dzw];
PIE_cat.Dzu = [PIE1.Dzu,PIE2.Dzu];
PIE_cat.gy = poly.gy;
PIE_cat.Dyw = [PIE1.Dyw,PIE2.Dyw];
PIE_cat.Dyu = [PIE1.Dyu,PIE2.Dyu];
% Algebraic constraints on the separate states, inputs, and actuators.
PIE_cat.K = blkdiag(PIE1.K,PIE2.K);
PIE_cat.Kw = blkdiag(PIE1.Kw,PIE2.Kw);
PIE_cat.Ku = blkdiag(PIE1.Ku,PIE2.Ku);

% Concatenate the components of the state and inputs, the outputs remain
% the same. The names of the inputs are assigned by initialize.
for type = {'x','w','u'}
    PIE_cat.space.(type{1}) = struct(...
        'dims',[PIE1.space.(type{1}).dims; PIE2.space.(type{1}).dims],...
        'depmat',[PIE1.space.(type{1}).depmat; PIE2.space.(type{1}).depmat],...
        'varname',{cell(0,1)});
end
PIE_cat.space.x.varname = names_x;
PIE_cat.space.z = PIE1.space.z;
PIE_cat.space.y = PIE1.space.y;

% Perform remaining concatenations.
if nargin>2
    PIE_cat = horzcat(PIE_cat,varargin{3:end});
end

end
