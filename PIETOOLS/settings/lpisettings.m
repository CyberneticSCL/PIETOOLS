function sttngs = lpisettings(type,derivative_strictness,simplify,solver)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sttngs = lpisettings(type,dim,derivative_strictness,simplify,solver)
% defines a "settings" structure for the PIETOOLS executive functions,
% calling one of the predefined settings scripts.
%
% INPUTS
%   type:   One of {'light','heavy','veryheavy','stripped','extreme','custom'}.
%           Different types will invoke different degrees for the monomials
%           in the opvar objects, as well as whether to add a Psatz term,
%           and whether to use "lpi_ineq" in solving LPIs. As such, the
%           type will determine the complexity/accuracy of the LPI, with
%           complexity increasing as:
%               extreme < stripped < light < heavy < veryheavy
%           Defaults to 'light'.
%   derivative_strictness:  Nonnegative scalar quantifying how strictly
%                           negative the derivative of the LF in the
%                           different executive functions is supposed to
%                           be. Set to 0 to enforce only nonngetavity of
%                           the derivative.
%                           Defaults to 0.
%   simplify:   '' or 'psimplify'. Set to 'psimplify' to perform
%               simplification of the SOS program produced by the
%               executives, before it is solved.
%               Defaults to ''.
%   solver:     One of {'sedumi','mosek','sdpnalplus','sdpt3'}, specifying
%               which SOS solver to use to solve the SOS program produced
%               by the executives. 
%               Defaults to 'sedumi'.
%
% OUTPUTS
%   sttngs: A struct specifying settings for the PIETOOLS executive
%           functions. See e.g. the "settings_light" function for more
%           details.
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - lpisettings
%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding SS - 15_07_2022

arguments
    type {mustBeMember(type,{'light','heavy','veryheavy','stripped','extreme','custom'})}='light';
    derivative_strictness {mustBeNonnegative}= 0;
    simplify {mustBeMember(simplify,{'','psimplify'})}= '';
    solver {mustBeMember(solver,{'sedumi','mosek','sdpnalplus','sdpt3'})} = 'sedumi';
end

switch type
    % Set the general LPI settings for 1D case, adding 2D settings as a
    % separate field.
    case 'light'
        sttngs = settings_PIETOOLS_light;
        sttngs.settings_2d = settings_PIETOOLS_light_2D;
    case 'custom'
        sttngs = settings_PIETOOLS_custom;
        sttngs.settings_2d = settings_PIETOOLS_custom_2D;
    case 'extreme'
        sttngs = settings_PIETOOLS_extreme;
        sttngs.settings_2d = settings_PIETOOLS_extreme_2D;
    case 'heavy'
        sttngs = settings_PIETOOLS_heavy;
        sttngs.settings_2d = settings_PIETOOLS_heavy_2D;
    case 'stripped'
        sttngs = settings_PIETOOLS_stripped;
        sttngs.settings_2d = settings_PIETOOLS_stripped_2D;
    case 'veryheavy'
        sttngs = settings_PIETOOLS_veryheavy;
        sttngs.settings_2d = settings_PIETOOLS_veryheavy_2D;
    otherwise
        warning('Unknown settings parameter requested. Defaulting to light settings');
        sttngs = settings_PIETOOLS_light;
        sttngs.settings_2d = settings_PIETOOLS_light_2D;
end

% Set strictness of positivity of LF and negativity of its derivative.
sttngs.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
sttngs.eppos2 = 1*1e-6;   % Positivity of Lyapunov Function with respect to spatially distributed states
sttngs.epneg = derivative_strictness;    % Negativity of Derivative of Lyapunov Function in both ODE and PDE state

% Set same settings for 2D case.
sttngs.settings_2d.eppos = [1e-4; 1e-6; 1e-6; 1e-6];    % Positivity of LF wrt R x L2[s1] x L2[s2] x L2[s1,s2]
sttngs.settings_2d.epneg = derivative_strictness;       % Negativity of Derivative of Lyapunov Function in both ODE and PDE state


% Set the SOS solve settings (independent of dimension).
sttngs.sos_opts.simplify = strcmp(simplify,'psimplify');    % Use psimplify in solving the SOS?
sttngs.sos_opts.solver = solver;    % Use psimplify in solving the SOS?

end