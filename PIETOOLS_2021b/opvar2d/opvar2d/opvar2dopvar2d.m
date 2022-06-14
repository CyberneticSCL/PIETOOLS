function D = opvar2dopvar2d(P,dvarname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D = opvar2dopvar2d(P) converts an opvar2d object to a corresponding
% dopvar2d object with decision variable names dvarname
% 
% INPUT
% P: opvar2d object (or dopvar2d object)
% dvarname: nx1 cell of chars (optional)
%       - Each element should correspond to a decision variable name for
%         the polynomial objects that appear in the opvar2d
%       - If no dvarname is specified, polynomial varibales with "coeff" in
%         the name will be converted to decision variables. If no such
%         variables appear, dpvars will be constructed with no decision
%         variables.
% 
% OUTPUT
% D: dopvar2d object equivalent of P
%       - Object has same spatial domain and dimensions as input P
%       - All components of D are dpvar class objects
%  
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, or Declan Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ, MP, SS - 08/05/2021

% Error checking
if ~isa(P,'opvar2d') && ~isa(P,'dopvar2d')
    error('Input must be a opvar2d (or dopvar2d) type object')
end

% Initialize the dopvar2d object with the same spatial domain and
% dimensions as the opvar2d object
D = dopvar2d([],P.dim,P.I,P.var1,P.var2);

if nargin==1
    % If no decision variable names are specified just use "dpvar" to convert
    % each component to a dpvar
    fset = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
    for f=fset
        if ~isempty(P.(f{:}))
            D.(f{:}) = dpvar(P.(f{:}));
        else
            D.(f{:}) = [];
        end
    end
    fset = {'Rxx','Rx2','R2x'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(3,1);
        for i=1:3
            if ~isempty(tmp_inR{i,1})
                tmp_outR{i,1} = dpvar(tmp_inR{i,1});
            else
                tmp_outR{i,1} = [];
            end
        end
        D.(f{:}) = tmp_outR;
    end
    fset = {'Ryy','Ry2','R2y'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(1,3);
        for j=1:3
            if ~isempty(tmp_inR{1,j})
                tmp_outR{1,j} = dpvar(tmp_inR{1,j});
            else
                tmp_outR{1,j} = [];
            end
        end
        D.(f{:}) = tmp_outR;
    end
    fset = {'R22'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(3,3);
        for i=1:3
            for j=1:3
                if ~isempty(tmp_inR{i,j})
                    tmp_outR{i,j} = dpvar(tmp_inR{i,j});
                else
                    tmp_outR{i,j} = [];
                end
            end
        end
        D.(f{:}) = tmp_outR;
    end
    
elseif nargin==2
    % If decision variable names are specified, use poly2dpvar to convert
    % each element to a dpvar with specified dvarnames
    if ~iscellstr(dvarname) && ~ischar(dvarname)
        error('Decision variable names must be specified as a cell of ''char'' objects')
    end    
    fset = {'R00', 'R0x', 'R0y', 'R02', 'Rx0', 'Rxy', 'Ry0', 'Ryx', 'R20'};
    for f=fset
        if ~isempty(P.(f{:}))
            D.(f{:}) = poly2dpvar(polynomial(P.(f{:})),dvarname);
        else
            D.(f{:}) = [];
        end
    end
    fset = {'Rxx','Rx2','R2x'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(3,1);
        for i=1:3
            if ~isempty(tmp_inR{i,1})
                tmp_outR{i,1} = poly2dpvar(polynomial(tmp_inR{i,1}),dvarname);
            else
                tmp_outR{i,1} = [];
            end
        end
        D.(f{:}) = tmp_outR;
    end
    fset = {'Ryy','Ry2','R2y'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(1,3);
        for j=1:3
            if ~isempty(tmp_inR{1,j})
                tmp_outR{1,j} = poly2dpvar(polynomial(tmp_inR{1,j}),dvarname);
            else
                tmp_outR{1,j} = [];
            end
        end
        D.(f{:}) = tmp_outR;
    end
    fset = {'R22'};
    for f=fset
        tmp_inR = P.(f{:});
        tmp_outR = cell(3,3);
        for i=1:3
            for j=1:3
                if ~isempty(tmp_inR{i,j})
                    tmp_outR{i,j} = poly2dpvar(polynomial(tmp_inR{i,j}),dvarname);
                else
                    tmp_outR{i,j} = [];
                end
            end
        end
        D.(f{:}) = tmp_outR;
    end
    
elseif nargin>=3
    error('At most two inputs are allowed: an opvar2d to convert, and a cell of decision variable names.') 
end

end