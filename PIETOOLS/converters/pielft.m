function PIE_lft = pielft(PIE1,PIE2)
% PIE_LFT = PIELFT(PIE1,PIE2) takes the linear fractional transformation of
% two PIE systems. Specifically, if PIE1 and PIE2 are defined by
%
% d/dt (T1*v1(t)+Tw1*w1(t)+Tu1*u1(t)) = A1*v1(t) +Bw1*w1(t) +Bu1*u1(t)
%                               z1(t) = Cz1*v1(t)+Dzw1*w1(t)+Dzu1*u1(t)
%                               y1(t) = Cy1*v1(t)+Dyw1*w1(t)+Dyu1*u1(t) 
%
% and
% 
% d/dt (T2*v2(t)+Tw2*w2(t)+Tu2*u2(t)) = A2*v2(t) +Bw2*w2(t) +Bu2*u2(t)
%                               z2(t) = Cz2*v2(t)+Dzw2*w2(t)+Dzu2*u2(t)
%                               y2(t) = Cy2*v2(t)+Dyw2*w2(t)+Dyu2*u2(t) 
% 
% respectively, then pielft(PIE1,PIE2) is obtained by enforcing u1 = y2 and
% u2 = y1, taking the form
%
%       d/dt (Top*v(t)+Twop*w(t)) = Aop*v(t) + Bop*w(t)
%                            z(t) = Cop*v(t) + Dop*w(t)
%
% where v = [v1;v2], w = [w1;w2], and z = [z1;z2], and
%   Top = [T1+Tu1*Dyu2*Cy1,  Tu1*Cy2        ]      
%         [Tu2*Cy1,          T2+Tu2*Dyu1*Cy2], 
%  Twop = [Tw1+Tu1*Dyu2*Dyw1,  Tu1*Dyw2       ]
%         [Tu2*Dyw1,         Tw2+Tu2*Dyu1*Dyw2],
%   Aop = [A1+Bu1*Dyu2*Cy1,  Bu1*Cy2        ]      
%         [Bu2*Cy1,          A2+Bu2*Dyu1*Cy2],
%   Bop = [Bw1+Bu1*Dyu2*Dyw1,  Bu1*Dyw2         ]
%         [Bu2*Dyw1,           Bw2+Bu2*Dyu1*Dyw2],
%   Cop = [Cz1+Dzu1*Dyu2*Cy1,  Dzu1*Cy2         ] 
%         [Dzu2*Cy1,           Cz2+Dzu2*Dyu1*Cy2],
%   Dop = [Dzw1+Dzu1*Dyu2*Dyw1,  Dzu1*Dyw2          ]
%         [Dzu2*Dyw1,            Dzw2+Dzu2*Dyu1*Dyw2],
%
% INPUT
% - PIE1, PIE2:     'pie_struct' objects representing linear PIEs of which
%                   to take the linear fractional transformation. The
%                   function space of the observed outputs of PIE1 must
%                   match that of the controlled inputs of PIE2, and vice
%                   versa, so that e.g.
%                       PIE1.space.u = PIE2.space.y;
%                       PIE1.space.y = PIE2.space.u;
%                   in terms of the spaces on which the components exist.
%
% OUTPUT
% - PIE_lft:        'pie_struct' object representing the linear
%                   fractional transformation of PIE1 and PIE2, setting the
%                   controlled inputs of PIE1 equal to the observed outputs
%                   of PIE2, and vice versa. The state variables, exogenous
%                   inputs, and regulated outputs of PIE1 and PIE2 are
%                   concatenated, and PIE_LFT has no controlled inputs or
%                   observed outputs.
%
% NOTES
% - Only linear PIEs, of type 'pie_struct', are supported, and the PIEs
%   cannot involve an algebraic constraint (K, Kw, Ku must be empty).
% - The operators of PIE_LFT are 'copvar' objects, with one space for each
%   of the spaces of PIE1, followed by those of PIE2, so that the state is
%   v = [v1;v2], and the components of PIE_LFT.space are one for each of
%   these spaces (see also 'initialize');
% - Currently, all controlled inputs and observed outputs are treated as
%   interconnection signals, and thus the dimensions of the controlled
%   inputs of PIE1 must match those of the observed outputs of PIE2, and
%   vice versa. 
% - If e.g. PIE1 does not have any controlled inputs, but does have 
%   exogenous inputs, and the dimensions of these inputs match those of the
%   observed outputs of PIE2, then the exogenous inputs of PIE1 are assumed
%   to actually be controlled inputs, and PIE1 is assumed to have no
%   exogenous inputs. Similarly, if PIE1 does not have any observed
%   outputs, but does have regulated outputs, and the dimensions of those
%   outputs match those of the controlled inputs of PIE2, then the
%   regulated outputs of PIE1 are treated as the interconnection signal.
%   Same reasoning applies for PIE2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pielft
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
% DJ, 12/22/2024: Initial coding;
% DJ, 12/30/2024: Add support for nonzero feedthrough;
% DJ, 10/29/2026: Update to new pie_struct structure, throw an error for
%                   PIEs with an algebraic constraint, and correct the
%                   operators of PIE2 that are swapped if the regulated
%                   output of PIE2 is the interconnection signal;

% % Right now, support only exactly two arguments.
if nargin~=2
    error("Two input arguments must be specified.")
end

% % Check if PIEs are suitable for LFT.
% Make sure the PIEs are linear PIEs.
if ~isa(PIE1,'pie_struct') || ~isa(PIE2,'pie_struct')
    error("LFT is only supported for linear PIEs, of type 'pie_struct'.")
end
% Make sure PIEs are properly specified.
try 
    if ~PIE1.is_initialized                                                 % DJ, 10/29/2026
        PIE1 = initialize(PIE1);
    end
    if ~PIE2.is_initialized
        PIE2 = initialize(PIE2);
    end
catch
    error("One of the PIEs is not properly specified; LFT not supported.")
end
% The inputs of the constraint Ku*u+Kw*w+K*x=0 would be substituted by the
% interconnection, which is not supported for now.
if size(PIE1.K,1)>0 || size(PIE2.K,1)>0                                     % DJ, 10/29/2026
    error("LFT is not supported for PIEs with an algebraic constraint.")
end              
% Make sure the variables and domains match. A finite-dimensional PIE does
% not involve any spatial variables, so these only have to match if both PIEs
% involve spatial variables.                                                 % DJ, 10/29/2026
uses_vars = @(P) any(any([P.space.x.depmat; P.space.w.depmat; P.space.u.depmat;...
                          P.space.z.depmat; P.space.y.depmat]));
PIE_vars = PIE1;
if ~uses_vars(PIE1) && uses_vars(PIE2)
    PIE_vars = PIE2;
elseif uses_vars(PIE1) && uses_vars(PIE2)
    if PIE1.dim~=PIE2.dim
        error("PIEs are of different dimensionality; LFT currently not supported.")
    end
    if ~isequal(PIE1.space.vars,PIE2.space.vars)
        error("Spatial variables of PIEs must match.")
    end
    if any(any(PIE1.dom~=PIE2.dom))
        error("Spatial domains of PIEs must match.")
    end
end


% % Check that the input spaces of PIE1 match the output spaces of PIE2, and
% % vice versa.
has_u1 = any(PIE1.space.u.dims);
has_y1 = any(PIE1.space.y.dims);
has_u2 = any(PIE2.space.u.dims);
has_y2 = any(PIE2.space.y.dims);
use_z2 = false;     use_z1 = false;
use_w1 = false;     use_w2 = false;
if has_u1 && ~has_y2 && same_space(PIE1.space.u,PIE2.space.z)
    % Treat regulated output as interconnection signal, u1 = z2;
    use_z2 = true;
elseif ~has_u1 && has_y2 && same_space(PIE2.space.y,PIE1.space.w)
    % Treat exogenous input as interconnection signal, w1 = y2;
    use_w1 = true;
elseif ~same_space(PIE1.space.u,PIE2.space.y)
    error("Controlled input (u) dimensions of first PIE should match observed output (y) dimensions of second PIE.")
end
if has_u2 && ~has_y1 && same_space(PIE2.space.u,PIE1.space.z)
    % Treat regulated output as interconnection signal, u2 = z1;
    use_z1 = true;
elseif ~has_u2 && has_y1 && same_space(PIE1.space.y,PIE2.space.w)
    % Treat exogenous input as interconnection signal, w2 = y1;
    use_w2 = true;
elseif ~same_space(PIE2.space.u,PIE1.space.y)
    error("Controlled input (u) dimensions of second PIE should match observed output (y) dimensions of first PIE.")
end


% % Extract PI operators defining each PIE
[Top1,Twop1,Tuop1,Aop1,Bwop1,Buop1,...
      Czop1,Dzwop1,Dzuop1,Cyop1,Dywop1,Dyuop1] = extract_ops(PIE1);
[Top2,Twop2,Tuop2,Aop2,Bwop2,Buop2,...
      Czop2,Dzwop2,Dzuop2,Cyop2,Dywop2,Dyuop2] = extract_ops(PIE2);
if use_w1
    % Suppose u1 = w1, and set w1 = [];
    Twop1_new = Tuop1;      Tuop1 = Twop1;      Twop1 = Twop1_new;
    Bwop1_new = Buop1;      Buop1 = Bwop1;      Bwop1 = Bwop1_new;
    Dzwop1_new = Dzuop1;    Dzuop1 = Dzwop1;    Dzwop1 = Dzwop1_new;
    Dywop1_new = Dyuop1;    Dyuop1 = Dywop1;    Dywop1 = Dywop1_new;
end
if use_w2
    % Suppose u2 = w2, and set w2 = [];
    Twop2_new = Tuop2;      Tuop2 = Twop2;      Twop2 = Twop2_new;
    Bwop2_new = Buop2;      Buop2 = Bwop2;      Bwop2 = Bwop2_new;
    Dzwop2_new = Dzuop2;    Dzuop2 = Dzwop2;    Dzwop2 = Dzwop2_new;
    Dywop2_new = Dyuop2;    Dyuop2 = Dywop2;    Dywop2 = Dywop2_new;
end
if use_z1
    % Suppose y1 = z1, and set z1 = [];
    Czop1_new = Cyop1;      Cyop1 = Czop1;      Czop1 = Czop1_new;
    Dzwop1_new = Dywop1;    Dywop1 = Dzwop1;    Dzwop1 = Dzwop1_new;
    Dzuop1_new = Dyuop1;    Dyuop1 = Dzuop1;    Dzuop1 = Dzuop1_new;
end
if use_z2
    % Suppose y2 = z2, and set z2 = [];
    Czop2_new = Cyop2;      Cyop2 = Czop2;      Czop2 = Czop2_new;         % DJ, 10/29/2026 (was Cyop1)
    Dzwop2_new = Dywop2;    Dywop2 = Dzwop2;    Dzwop2 = Dzwop2_new;       % DJ, 10/29/2026 (was Dywop1)
    Dzuop2_new = Dyuop2;    Dyuop2 = Dzuop2;    Dzuop2 = Dzuop2_new;       % DJ, 10/29/2026 (was Dyuop1)
end
if ~(Dyuop1==0) && ~(Dyuop2==0)
    error("LFT for PIEs with feedthrough from controlled input to observed output is not supported.")
end


% % Construct the actual LFT PIE
PIE_lft = pie_struct();
%   PIE_lft.dom = PIE_vars.dom;                                             % DJ, 10/29/2026 (was)
%   PIE_lft.vars = PIE_vars.vars;                                           % DJ, 10/29/2026 (was)
PIE_lft.space.vars = PIE_vars.space.vars;                                   % DJ, 10/29/2026
PIE_lft.space.dom = PIE_vars.space.dom;                                     % DJ, 10/29/2026
PIE_lft.T = [Top1+Tuop1*Dyuop2*Cyop1,Tuop1*Cyop2;     
             Tuop2*Cyop1,Top2+Tuop2*Dyuop1*Cyop2];            
PIE_lft.Tw = [Twop1+Tuop1*Dyuop2*Dywop1,Tuop1*Dywop2;
              Tuop2*Dywop1,Twop2+Tuop2*Dyuop1*Dywop2];
PIE_lft.A = [Aop1+Buop1*Dyuop2*Cyop1,Buop1*Cyop2;      
             Buop2*Cyop1,Aop2+Buop2*Dyuop1*Cyop2];           
PIE_lft.Bw = [Bwop1+Buop1*Dyuop2*Dywop1,Buop1*Dywop2;
              Buop2*Dywop1,Bwop2+Buop2*Dyuop1*Dywop2];
PIE_lft.Cz = [Czop1+Dzuop1*Dyuop2*Cyop1,Dzuop1*Cyop2;
              Dzuop2*Cyop1,Czop2+Dzuop2*Dyuop1*Cyop2];
PIE_lft.Dzw = [Dzwop1+Dzuop1*Dyuop2*Dywop1,Dzuop1*Dywop2;
               Dzuop2*Dywop1,Dzwop2+Dzuop2*Dyuop1*Dywop2];

PIE_lft = initialize(PIE_lft);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function logval = same_space(sp1,sp2)
% LOGVAL = SAME_SPACE(SP1,SP2) checks whether two function spaces,
% specified as e.g. PIE1.space.u and PIE2.space.y, are the same, in the
% sense that for each combination of spatial variables, the total number
% of components depending on exactly these variables is the same.
%
% DJ, 10/29/2026: Initial coding;
nvars = size(sp1.depmat,2);
logval = true;
for ii=1:2^nvars
    % Components depending on exactly the variables indicated by ii.
    depmat = logical(bitget(ii-1,1:nvars));
    in1 = all(sp1.depmat==depmat,2);      
    in2 = all(sp2.depmat==depmat,2);
    logval = logval && sum(sp1.dims(in1))==sum(sp2.dims(in2));
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [Top,Twop,Tuop,Aop,Bwop,Buop,...
            Czop,Dzwop,Dzuop,Cyop,Dywop,Dyuop] = extract_ops(PIE)
% Extracts the PI operators from a 'pie_struct' object.

Top = PIE.T;      Twop = PIE.Tw;      Tuop = PIE.Tu;
Aop = PIE.A;      Bwop = PIE.Bw;      Buop = PIE.Bu;
Czop = PIE.Cz;    Dzwop = PIE.Dzw;    Dzuop = PIE.Dzu;
Cyop = PIE.Cy;    Dywop = PIE.Dyw;    Dyuop = PIE.Dyu;

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %