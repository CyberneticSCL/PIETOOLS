function PIE_cat = vertcat(varargin)
% PIE_CAT = VERTCAT(PIE1,PIE2,...) takes a set of 'pie_struct' objects and
% returns a 'pie_struct' object representing the vertical concatentation
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
% - PIE_out:    'pie_struct' object representing the PIE defined by the
%               vertical concatenation of the input systems.
%
% NOTES
% For two PIE systems PIE_i defined as
% d/dt (Topi*v1(t)+Twopi*wi(t)+Tuopi*ui(t)) = Aopi*vi(t) +Bwopi*wi(t) +Buopi*ui(t)
%                                      z(t) = Czopi*vi(t)+Dzwopi*wi(t)+Dzuopi*ui(t)
%                                      y(t) = Cyopi*vi(t)+Dywopi*wi(t)+Dyuopi*ui(t)
% the horizontal concatenation PIE = [PIE1,PIE2] takes the form
% d/dt (Top*v(t)+Twop*w(t)+Tuop*u(t)) = Aop*v(t) +Bwop*w(t) +Buop*u(t)
%                                z(t) = Czop*v(t)+Dzwop*w(t)+Dzuop*u(t)
%                                y(t) = Cyop*v(t)+Dywop*w(t)+Dyuop*u(t)
% where v = [v1;v2], u = [u1;u2], w = [w1;w2], and
%   Top = [Top1,0;0,Top2];  Twop = [Twop1;Twop2];       Tuop = [Tuop1; Tuop2];
%   Aop = [Aop1,0;0,Aop2];  Bwop = [Bwop1;Bwop2];       Buop = [Buop1; Buop2];
%   Czop = [Czop1,0;0,Czop2];   Dzwop = [Dzwop1;Dzwop2];    Dzuop = [Dzuop1;Dzuop2];
%   Cyop = [Cyop1,0;0,Cyop2];   Dywop = [Dywop1;Dywop2];    Dyuop = [Dyuop1;Dyuop2];
% Keep in mind that, since 'opvar' objects can only map to and from 
% R^n x L2^m[a,b] -- not to or from R^n1 x L2^m1[a,b] x R^n2 x L2^m2[a,b]
% -- the order of state and output variables may be exchanged, as the
% finite-dimensional components are grouped together, and so are the
% infinite-dimensional components. That is, if e.g. v1=[v11;v12] and
% v2=[v21;v22] where v11 and v21 are finite-dimensional and v12 and v22 are
% infinite-dimensional, then v = [v11;v21;v12;v22];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - vertcat
%
% Copyright (C)2024  PIETOOLS Team
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

% % Extract the inputs
if nargin==1
    PIE_cat = varargin{1};
    return
end
PIE1 = varargin{1};
PIE2 = varargin{2};


% % Check if concatenation is supported
% Make sure PIEs are properly specified.
try PIE1 = initialize(PIE1);    PIE2 = initialize(PIE2);
catch
    error("One of the PIEs is not properly specified, concatenation not supported.")
end
% Make sure variables and domains match.
if any(any(~isequal(PIE1.vars,PIE2.vars)))
    error("Spatial variables of PIEs must match.")
end
if any(any(PIE1.dom~=PIE2.dom))
    error("Spatial domains of PIEs must match.")
end
% Make sure input dimensions match.
if any(PIE1.Bu.dim(:,2)~=PIE2.Bu.dim(:,2))
    error("Dimensions of observed outputs are not consistent.")
end
if any(PIE1.Bw.dim(:,2)~=PIE2.Bw.dim(:,2))
    error("Dimensions of regulated outputs are not consistent.")
end


% % Construct the concatenated PIE structure
PIE_cat = pie_struct();
PIE_cat.vars = PIE1.vars;
PIE_cat.dom = PIE1.dom;

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

% Perform remaining concatenations.
if nargin>2
    PIE_cat = vertcat(PIE_cat,varargin{3:end});
end

end