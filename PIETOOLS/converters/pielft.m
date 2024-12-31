function PIE_lft = pielft(PIE1,PIE2)
% PIE_LFT = PIELFT(PIE1,PIE2) takes the linear fractional transformation of
% two PIE systems. Specifically, if PIE1 and PIE2 are defined by
% d/dt (Top1*v1(t)+Twop1*w1(t)+Tuop1*u1(t)) = Aop1*v1(t) +Bwop1*w1(t) +Buop1*u1(t)
%                                     z1(t) = Czop1*v1(t)+Dzwop1*w1(t)+Dzuop1*u1(t)
%                                     y1(t) = Cyop1*v1(t)+Dywop1*w1(t)+Dyuop1*u1(t) 
% and
% d/dt (Top2*v2(t)+Twop2*w2(t)+Tuop2*u2(t)) = Aop2*v2(t) +Bwop2*w2(t) +Buop2*u2(t)
%                                     z2(t) = Czop2*v2(t)+Dzwop2*w2(t)+Dzuop2*u2(t)
%                                     y2(t) = Cyop2*v2(t)+Dywop2*w2(t)+Dyuop2*u2(t) 
% respectively, then pielft(PIE1,PIE2) is obtained by enforcing u1 = y2 and
% u2 = y1, taking the form
%       d/dt (Top*v(t)+Twop*w(t)) = Aop*v(t) + Bop*w(t)
%                            z(t) = Cop*v(t) + Dop*w(t)
% where v = [v1;v2], w = [w1;w2], and z = [z1;z2], and
%   Top = [Top1+Tuop1*Dyuop2*Cyop1, Tuop1*Cyop2            ]      
%         [Tuop2*Cyop1,             Top2+Tuop2*Dyuop1*Cyop2], 
%  Twop = [Twop1+Tuop1*Dyuop2*Dywop1, Tuop1*Dywop2             ]
%         [Tuop2*Dywop1,              Twop2+Tuop2*Dyuop1*Dywop2],
%   Aop = [Aop1+Buop1*Dyuop2*Cyop1, Buop1*Cyop2            ]      
%         [Buop2*Cyop1,             Aop2+Buop2*Dyuop1*Cyop2],
%   Bop = [Bwop1+Buop1*Dyuop2*Dywop1, Buop1*Dywop2             ]
%         [Buop2*Dywop1,              Bwop2+Buop2*Dyuop1*Dywop2],
%   Cop = [Czop1+Dzuop1*Dyuop2*Cyop1, Dzuop1*Cyop2             ] 
%         [Dzuop2*Cyop1,              Czop2+Dzuop2*Dyuop1*Cyop2],
%   Dop = [Dzwop1+Dzuop1*Dyuop2*Dywop1, Dzuop1*Dywop2              ]
%         [Dzuop2*Dywop1,               Dzwop2+Dzuop2*Dyuop1*Dywop2],
%
% INPUT
% - PIE1, PIE2:     'pie_struct' objects representing PIEs of which to take
%                   the linear fractional transformation. The number of
%                   observed outputs of PIE1 must match the number of
%                   controlled inputs of PIE2, and vice versa, so that e.g.
%                       PIE1.Bu.dim(:,2)=PIE2.Cy.dim(:,1);
%                       PIE1.Cy.dim(:,1)=PIE2.Bu.dim(:,2);
%
% OUTPUT
% - PIE_lft:        'pie_struct' object representing to the linear
%                   fractional transformation of PIE1 and PIE2, setting the
%                   controlled inputs of PIE1 equal to the observed outputs
%                   of PIE2, and vice versa. The state variables, exogenous
%                   inputs, and regulated outputs of PIE1 and PIE2 are
%                   concatenated, and PIE_LFT has no controlled inputs or
%                   observed outputs.
%
% NOTES
% - Keep in mind that, since 'opvar' objects can only map to and from 
%   R^n x L2^m[a,b] -- not to or from R^n1 x L2^m1[a,b] x R^n2 x L2^m2[a,b]
%   -- the order of state, input, and output variables may be exchanged, as
%   the finite-dimensional components are grouped together, and so are the
%   infinite-dimensional components. That is, if e.g. v1=[v11;v12] and
%   v2=[v21;v22] where v11 and v21 are finite-dimensional and v12 and v22
%   are infinite-dimensional, then the LFT has state v = [v11;v21;v12;v22];
% - Currently, all controlled inputs and observed outputs are treated as
%   interconnection signals, and thus the dimensions of the controlled
%   inputs of PIE1 must match those of the observed outputs of PIE2, and
%   vice versa. 
% - If e.g. PIE1 does not have any controlled inputs, but does have 
%   exogenous inputs, and the dimensions of these inputs match those of the
%   observed outputs of PIE2, then the exogenous inputs of PIE1 are assumed
%   to actually be controlled inputs, and PIE1 is assumed to have no
%   exogenous inputs. Similarly, if PIE1 does not have any observed outputs,
%   but does have regulated outputs, and the dimensions of those outputs
%   match those of the controlled inputs of PIE2, then the regulated
%   outputs of PIE1 are treated as the interconnection signal. Same
%   reasoning applies for PIE2. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pielft
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
% DJ, 12/22/2024: Initial coding;
% DJ, 12/30/2024: Add support for nonzero feedthrough;

% % Right now, support only exactly two arguments.
if nargin~=2
    error("Two input arguments must be specified.")
end

% % Check if PIEs are suitable for LFT.
% Make sure PIEs are properly specified.
try PIE1 = initialize(PIE1);    PIE2 = initialize(PIE2);
catch
    error("One of the PIEs is not properly specified; LFT not supported.")
end
% Make sure the variables and domains match.
if PIE1.dim~=PIE2.dim
    error("PIEs are of different dimensionality; LFT currently not supported.")
end
if any(any(~isequal(PIE1.vars,PIE2.vars)))
    error("Spatial variables of PIEs must match.")
end
if any(any(PIE1.dom~=PIE2.dom))
    error("Spatial domains of PIEs must match.")
end


% % Check that the input dimensions of PIE1 match the output dimensions of
% % PIE2, and vice versa.
nu1_op = PIE1.Bu.dim(:,2);        nw1_op = PIE1.Bw.dim(:,2);
ny1_op = PIE1.Cy.dim(:,1);        nz1_op = PIE1.Cz.dim(:,1);
nu2_op = PIE2.Bu.dim(:,2);        nw2_op = PIE2.Bw.dim(:,2);
ny2_op = PIE2.Cy.dim(:,1);        nz2_op = PIE2.Cz.dim(:,1);
use_z2 = false;     use_z1 = false;
use_w1 = false;     use_w2 = false;
if any(nu1_op) && ~any(ny2_op) && all(nu1_op==nz2_op)
    % Treat regulated output as interconnection signal, u1 = z2;
    use_z2 = true;
elseif ~any(nu1_op) && any(ny2_op) && all(ny2_op==nw1_op)
    % Treat exogenous input as interconnection signal, w1 = y2;
    use_w1 = true;
elseif any(nu1_op~=ny2_op)
    error("Controlled input (u) dimensions of first PIE should match observed output (y) dimensions of second PIE.")
end
if any(nu2_op) && ~any(ny1_op) && all(nu2_op==nz1_op)
    % Treat regulated output as interconnection signal, u2 = z1;
    use_z1 = true;
elseif ~any(nu2_op) && any(ny1_op) && all(ny1_op==nw2_op)
    % Treat exogenous input as interconnection signal, w2 = y1;
    use_w2 = true;
elseif any(nu2_op~=ny1_op)
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
    Czop2_new = Cyop1;      Cyop2 = Czop2;      Czop2 = Czop2_new;
    Dzwop2_new = Dywop1;    Dywop2 = Dzwop2;    Dzwop2 = Dzwop2_new;
    Dzuop2_new = Dyuop1;    Dyuop2 = Dzuop2;    Dzuop2 = Dzuop2_new;
end
if ~(Dyuop1==0) && ~(Dyuop2==0)
    error("LFT for PIEs with feedthrough from controlled input to observed output is not supported.")
end


% % Construct the actual LFT PIE
PIE_lft = pie_struct;
PIE_lft.dom = PIE1.dom;     
PIE_lft.vars = PIE1.vars;
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



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [Top,Twop,Tuop,Aop,Bwop,Buop,...
            Czop,Dzwop,Dzuop,Cyop,Dywop,Dyuop] = extract_ops(PIE)
% Extracts the PI operators from a 'pie_struct' object.

Top = PIE.T;      Twop = PIE.Tw;      Tuop = PIE.Tu;
Aop = PIE.A;      Bwop = PIE.Bw;      Buop = PIE.Bu;
Czop = PIE.Cz;    Dzwop = PIE.Dzw;    Dzuop = PIE.Dzu;
Cyop = PIE.Cy;    Dywop = PIE.Dyw;    Dyuop = PIE.Dyu;

end