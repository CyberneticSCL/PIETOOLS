function [prog_sol, varargout] = lpiscript_sop(PIE,lpi,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG_SOL,VARARGOUT] = LPISCRIPT_SOP(PIE,LPI,OPTS) solves the predefined
% LPI named LPI for the system PIE with the container executives
% (executives_sopvar), and returns the solved program and the outputs of
% the executive. It is LPISCRIPT on the container path:
%   LPI                 executive, form
%   'stability'         PIETOOLS_stability_sop, 'Q'
%   'stability-dual'    PIETOOLS_stability_sop, 'Qd'
%   'well-posedness'    PIETOOLS_well_posedness_sop
%   'l2gain'            PIETOOLS_Hinf_gain_sop, 'Q'
%   'l2gain-dual'       PIETOOLS_Hinf_gain_sop, 'Qd'
%   'h2norm'            PIETOOLS_H2_norm_sop, 'o'
%   'h2norm-dual'       PIETOOLS_H2_norm_sop, 'c'
%   'hinf-observer'     PIETOOLS_Hinf_estimator_sop
%   'hinf-controller'   PIETOOLS_Hinf_control_sop
%   'h2-observer'       PIETOOLS_H2_estimator_sop
%   'h2-controller'     PIETOOLS_H2_control_sop
% the executive and form of each stock executive that LPISCRIPT calls.
%
% INPUT
% - PIE:   'pie_struct' object of a 1-D or 2-D PIE; or a 'pde_struct'
%          object, or a 'sys' object of type 'pde' or 'pie', which is
%          converted to a PIE first;
% - lpi:   'char' or 'string', one of the names above; or a function handle
%          f, called as [prog_sol, vout] = f(PIE,OPTS), with vout a 'cell'
%          of further outputs;
% - opts:  'struct' of LPI settings (lpisettings_sop or lpisettings); or a
%          'char' or 'string' naming a level: 'stripped', 'light', 'heavy',
%          'veryheavy' or 'custom' (lpisettings_sop), or 'extreme'
%          (lpisettings); or [] for lpisettings_sop('light').
%
% OUTPUT
% - prog_sol:   'struct', the solved LPI program;
% - varargout:  the outputs of the executive after the program, in its
%               order: those LPISCRIPT returns (P for 'stability'; P, R,
%               omega for 'well-posedness'; R, gam for 'l2gain', R the
%               solved positive operator; W, gam, R, Q for 'h2norm'; L or
%               K, gam, P, Z for the H-infinity observer and controller; L
%               or K, gam, P, Z, W for the H2 observer and controller),
%               then the 'struct' info; for a function handle, the entries
%               of vout.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpiscript_sop
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
% Initial coding MMP, 10/28/2026: the container counterpart of LPISCRIPT
%                (DJ, 11/01/2024), which had none (maintainer: create
%                one); the same names, inputs and error messages, the
%                executives and forms of PIETOOLS_auto_execute_sop.
%                'custom', which LPISCRIPT accepts but has no case for
%                (it then indexes the char), is rejected; an unknown
%                settings name raises an error listing the accepted ones
%                (LPISCRIPT passes it to lpisettings, whose list lacks the
%                container tiers).

% Three inputs, as LPISCRIPT
if nargin<3
    error('Insufficient number of arguments. Correct syntax: lpiscript_sop(PIE,lpi,settings)');
end

% The LPI: a function handle, or one of the names
names = {'stability','stability-dual','well-posedness','l2gain','l2gain-dual', ...
         'h2norm','h2norm-dual','hinf-observer','hinf-controller','h2-observer','h2-controller'};
if ~isa(lpi,'function_handle')
    if ~(isa(lpi,'string') || isa(lpi,'char'))
        error("The input 'lpi' must be a string value or a function handle.");
    end
    lpi = char(lpi);
    if ~ismember(lpi,names)
        error("Specified LPI type must be one of 'stability', 'stability-dual', 'well-posedness'," + ...
              " 'l2gain', 'l2gain-dual', 'h2norm', 'h2norm-dual'," + ...
              " 'hinf-observer', 'hinf-controller', 'h2-observer', or 'h2-controller'.")
    end
end

% The settings: a struct as given, a level by name (the container tiers of
% lpisettings_sop; 'extreme', which they lack, from lpisettings), or the
% light tier
if isempty(opts)
    opts = lpisettings_sop('light');
elseif isa(opts,'string') || isa(opts,'char')
    opts = char(opts);
    if ismember(opts,{'stripped','light','heavy','veryheavy','custom'})
        opts = lpisettings_sop(opts);
    elseif strcmp(opts,'extreme')
        opts = lpisettings(opts);
    else
        error("Settings should be one of 'stripped', 'light', 'heavy', 'veryheavy', 'custom' or 'extreme'.")
    end
elseif ~isa(opts,'struct')
    error("Settings must either be a string value or a settings structure similar to the output of lpisettings_sop().")
end

% The system: a PIE, or a PDE converted to one
if isa(PIE,'sys') && strcmpi(PIE.type,'pie')
    PIE = PIE.params;
elseif (isa(PIE,'sys') && strcmpi(PIE.type,'pde')) || isa(PIE,'pde_struct')
    PIE = convert(PIE,'pie');
elseif ~isa(PIE,'pie_struct')
    error('PIE system for which to run the executive should be specified as object of type ''pie_struct''.')
end

% The executive and form of the LPI; the outputs after the program are
% passed on in the executive's order
if isa(lpi,'function_handle')
    [prog_sol, vout] = lpi(PIE,opts);
    varargout = reshape(vout,1,[]);
    return
end
switch lpi
    case 'stability',       f = @() PIETOOLS_stability_sop(PIE,opts,'Q');           nout = 3;
    case 'stability-dual',  f = @() PIETOOLS_stability_sop(PIE,opts,'Qd');          nout = 3;
    case 'well-posedness',  f = @() PIETOOLS_well_posedness_sop(PIE,opts);          nout = 5;
    case 'l2gain',          f = @() PIETOOLS_Hinf_gain_sop(PIE,opts,'Q');           nout = 4;
    case 'l2gain-dual',     f = @() PIETOOLS_Hinf_gain_sop(PIE,opts,'Qd');          nout = 4;
    case 'h2norm',          f = @() PIETOOLS_H2_norm_sop(PIE,opts,'o');             nout = 6;
    case 'h2norm-dual',     f = @() PIETOOLS_H2_norm_sop(PIE,opts,'c');             nout = 6;
    case 'hinf-observer',   f = @() PIETOOLS_Hinf_estimator_sop(PIE,opts);          nout = 6;
    case 'hinf-controller', f = @() PIETOOLS_Hinf_control_sop(PIE,opts);            nout = 6;
    case 'h2-observer',     f = @() PIETOOLS_H2_estimator_sop(PIE,opts);            nout = 7;
    case 'h2-controller',   f = @() PIETOOLS_H2_control_sop(PIE,opts);              nout = 7;
end
out = cell(1,nout);
[out{:}] = f();
prog_sol = out{1};
varargout = out(2:end);
end
