function Pop = sdopvar2dopvar(Psop)
% POP = SDOPVAR2DOPVAR(PSOP) takes an sdopvar object representing a
% 4-PI operator component with decision variables and returns a dopvar
% object representing the same operator.
%
% INPUTS
% - Psop:       'sdopvar' object representing a 1D PI decision operator. It
%               can map between different function spaces, but the input
%               and output domains cannot be both finite- and infinite-
%               dimensional (i.e. L2^n and R^m are supported, but R^m x
%               L2^n is not);
%
% OUTPUTS
% - Pop:        'dopvar' object representing the same operator as the
%               input;
%
% NOTES
% Modeled after sopvar2opvar. The primary variable is taken from
% objSdopvar.vars and the dummy is that name with '_dum' appended.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - sdopvar2dopvar
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
% DJ, 10/05/2026: Initial coding;

if ~isa(Psop,'sdopvar')
    error("Input must be of type 'sdopvar'.")
end

% Extract the dimension, variables, domain and decision variables
dims = Psop.dims;
vars = Psop.vars;
dom = Psop.dom;
Zd = Psop.Zd;

% Make sure the operator is 1D; no 2D decision-operator counterpart exists
if numel(unique([vars.in,vars.out]))>1
    error("Operator maps between functions of more than one variable, which is not supported.")
end

% The primary variable is recorded on the object, so use it. NB: 'dopvar
% Pop' below does evalin('caller','pvar s1 s1_dum'), so do not name this s1.
sv = unique([vars.in,vars.out]);
if isempty(sv),  vname = 's1';  else,  vname = char(sv{1});  end
dname = [vname,'_dum'];

% Initialize an empty operator
dopvar Pop;
Pop.var1 = polynomial({vname});
Pop.var2 = polynomial({dname});

% Distinguish cases of different input/output domains
if isempty(vars.in) && isempty(vars.out)
    % Operator maps R to R
    P = struct('A',Psop.params.A{1},'B',Psop.params.B{1});
    Pop.P = sdvar2dpvar(P,dims,struct('out',{cell(1,0)},'in',{cell(1,0)}),cell(1,0),cell(1,0),Zd);
elseif isempty(vars.out)
    % Operator maps L2 to R
    Pop.I = dom.in;
    P = struct('A',Psop.params.A{1},'B',Psop.params.B{1});
    Pop.Q1 = sdvar2dpvar(P,dims,struct('out',{cell(1,0)},'in',{{vname}}),cell(1,0),Psop.ZR,Zd);
elseif isempty(vars.in)
    % Operator maps R to L2
    Pop.I = dom.out;
    P = struct('A',Psop.params.A{1},'B',Psop.params.B{1});
    Pop.Q2 = sdvar2dpvar(P,dims,struct('out',{{vname}},'in',{cell(1,0)}),Psop.ZL,cell(1,0),Zd);
else
    % Operator maps L2 to L2. 'sdvar2dpvar' builds an explicit two-variable
    % dpvar, so the output and (dummy) input roles need DISTINCT names here
    % -- unlike 'sdopvar.vars', which uses the same physical name for both
    % (the dummy is implicit there; see dopvar2sdopvar.m's R-block comment).
    Pop.I = dom.in;
    vars_blk = struct('out',{{vname}},'in',{{dname}});
    P0 = struct('A',Psop.params.A{1},'B',Psop.params.B{1});
    P1 = struct('A',Psop.params.A{2},'B',Psop.params.B{2});
    P2 = struct('A',Psop.params.A{3},'B',Psop.params.B{3});
    Pop.R.R0 = sdvar2dpvar(P0,dims,vars_blk,Psop.ZL,Psop.ZR,Zd);
    Pop.R.R1 = sdvar2dpvar(P1,dims,vars_blk,Psop.ZL,Psop.ZR,Zd);
    Pop.R.R2 = sdvar2dpvar(P2,dims,vars_blk,Psop.ZL,Psop.ZR,Zd);
end

end
