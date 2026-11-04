function PIE_out = npie2pie_struct(PIE_in)
% PIE_OUT = NPIE2PIE_STRUCT(PIE_IN) takes a nonlinear PIE in 'npie_struct'
% format, and returns a 'pie_struct' object representing the same PIE, if
% the PIE is linear.
%
% INPUTS
% - PIE_in:     'npie_struct' object representing a PIE in any number of
%               spatial variables (it is initialized if it is not yet),
%                   d/dt (T*x+Tw*w+Tu*u) = f(x) + Bw*w + Bu*u,
%                                      z = gz(x) + Dzw*w + Dzu*u,
%                                      y = gy(x) + Dyw*w + Dyu*u,
%                   Ku*u + Kw*w + K*x = 0,
%               with the polynomials f, gz, and gy linear in the state x,
%               and the operators T, ..., Ku 'tensopmat' objects involving
%               only 'sopvar' objects;
%
% OUTPUTS
% - PIE_out:    'pie_struct' representing the same PIE as the input,
%                   d/dt (T*x+Tw*w+Tu*u) = A*x + Bw*w + Bu*u,
%                                      z = Cz*x + Dzw*w + Dzu*u,
%                                      y = Cy*x + Dyw*w + Dyu*u,
%                   Ku*u + Kw*w + K*x = 0,
%               with all operators as 'copvar' objects, mapping between the
%               components of the state, inputs and outputs of PIE_in.space,
%               each a separate space of the 'copvar'. The function spaces
%               and maps between the PDE and the PIE (space_pde, pie2pde,
%               pde2pie) and PIE_in.misc are copied;
%
% NOTES
% - An error is thrown if the PIE is nonlinear, i.e. if f, gz, or gy has a
%   monomial in the state of degree different from 1 with a coefficient
%   above 1e-12, or if an operator involves tensor products of operators,
%   see `tensopmat2copvar`.
% - This is the inverse of `pie2npie_struct`.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - npie2pie_struct
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
% DJ, 11/03/2026: Initial coding;


if ~PIE_in.is_initialized
    PIE_in = initialize(PIE_in);
end
obj_types = {'x','w','u','y','z'};
space = PIE_in.space;
space_pde = PIE_in.space_pde;

% Copy the function spaces, without the names of the components, which are
% not part of the 'pie_struct'.
PIE_out = pie_struct();
PIE_out.space.vars = space.vars;
PIE_out.space.dom = space.dom;
PIE_out.space_pde.vars = space_pde.vars;
PIE_out.space_pde.dom = space_pde.dom;
for obj = obj_types
    for f = {'dims','depmat'}
        PIE_out.space.(obj{1}).(f{1}) = space.(obj{1}).(f{1});
    end
    for f = {'dims','depmat','ID'}
        PIE_out.space_pde.(obj{1}).(f{1}) = space_pde.(obj{1}).(f{1});
    end
end

% Convert the operators to 'copvar' objects, with a space for each component
ops = {'T','Tw','Tu','Bw','Bu','Dzw','Dzu','Dyw','Dyu','K','Kw','Ku'};
for k=1:numel(ops)
    P = PIE_in.(ops{k});
    if size(P,1)>0 && size(P,2)>0
        PIE_out.(ops{k}) = tensopmat2copvar(P);
    end
end

% Convert the polynomials f, gz, and gy to the operators A, Cz, and Cy,
% acting on the state, so that f(x) = A*x. If the polynomial is linear, each
% monomial is a single component x_k of the state, and the coefficient of
% this monomial is the block of A acting on x_k.
polys = {'f','A';   'gz','Cz';  'gy','Cy'};
[vsort,idx] = sort(space.vars);
for k=1:size(polys,1)
    f = PIE_in.(polys{k,1});
    degs = zeros(0,1);
    if size(f,1)>0
        degs = sum(f.Z>0,2);
    end
    % Monomials of degree other than 1 may only have zero coefficients (e.g.
    % kept by `clean_poly` to retain the shape of a row), so that f is linear.
    for j = find(degs'~=1)
        for i = 1:size(f.C.ops,1)
            op_ij = f.C.ops{i,j};
            if isempty(op_ij)
                continue
            elseif isa(op_ij,'tensopvar') || isa(op_ij,'intop')
                [~,is_zero_ij] = clean_op(op_ij,1e-12);
            else
                [~,is_zero_ij] = clean_dpvar(op_ij,1e-12);
            end
            if ~is_zero_ij
                error(['The PIE is nonlinear: the polynomial ',polys{k,1},' has a monomial of degree ',num2str(degs(j)),' in the state x.'])
            end
        end
    end
    if ~any(degs==1)
        % A zero polynomial is a zero operator, which is only declared for A
        % as the other operators are initialized as zero if not specified.
        if k==1
            PIE_out.A = 0*PIE_out.T;
        end
        continue
    end
    % Place the coefficient of each monomial in the column of the component.
    is_lin = degs==1;
    [~,kcol] = ismember(f.varname(sum(f.Z(is_lin,:),2)),space.x.varname);
    if numel(unique(kcol))<numel(kcol)
        error(['The polynomial ',polys{k,1},' has multiple monomials in the same component of the state, merge these first, see `combine_terms`.'])
    end
    Pt = f.C;
    Pt.ops = Pt.ops(:,is_lin);
    ops_cell = cell(size(Pt.ops,1),numel(space.x.dims));
    ops_cell(:,kcol) = Pt.ops;
    Pt.ops = ops_cell;
    % The input space consists of all components of the state, also those
    % on which the polynomial does not depend.
    Pt.dims.in = space.x.dims;
    Pt.vars.in = vsort;
    Pt.dom.in = space.dom(idx,:);
    Pt.depmat.in = double(space.x.depmat(:,idx));
    PIE_out.(polys{k,2}) = tensopmat2copvar(Pt);
end

% Copy the maps between the PDE and the PIE, converting the operators
for obj = obj_types
    p2p = PIE_in.pde2pie.(obj{1});
    for f = setdiff(fieldnames(p2p)','C','stable')
        PIE_out.pde2pie.(obj{1}).(f{1}) = p2p.(f{1});
    end
    if isfield(p2p,'C') && isa(p2p.C,'tensopmat') && size(p2p.C,1)>0
        PIE_out.pde2pie.(obj{1}).C = tensopmat2copvar(p2p.C);
    end
    if isa(PIE_in.pie2pde.(obj{1}),'tensopmat') && size(PIE_in.pie2pde.(obj{1}),1)>0
        PIE_out.pie2pde.(obj{1}) = tensopmat2copvar(PIE_in.pie2pde.(obj{1}));
    end
end
PIE_out.misc = PIE_in.misc;
PIE_out = initialize(PIE_out);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

