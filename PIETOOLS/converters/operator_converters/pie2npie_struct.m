function PIE_out = pie2npie_struct(PIE_in)
% PIE_OUT = PIE2NPIE_STRUCT(PIE_IN) takes a linear PIE in 'pie_struct'
% format and returns an 'npie_struct' object representing the same PIE,
% with the right-hand sides of the PIE expressed as (linear) distributed
% polynomials.
%
% INPUTS
% - PIE_in:     'pie_struct' object representing a linear PIE in any
%               number of spatial variables, with the operators as
%               'copvar' objects (it is initialized if it is not yet);
%
% OUTPUTS
% - PIE_out:    'npie_struct' representing the same PIE as the input,
%                   d/dt (T*x+Tw*w+Tu*u) = f(x) + Bw*w + Bu*u,
%                                     z = gz(x) + Dzw*w + Dzu*u,
%                                     y = gy(x) + Dyw*w + Dyu*u,
%                   Ku*u + Kw*w + K*x = 0,
%               with the polynomials f = A*x, gz = Cz*x, gy = Cy*x, as
%               'polyopvar' objects, and the operators T, ..., Ku as
%               'tensopmat' objects. The components of the state x, the
%               inputs, and the outputs are those of PIE_in.space, with the
%               state components named 'x1_pie', ..., as returned by
%               `convert_PIETOOLS_PDE`, needed for e.g. `piesos_LF`. The
%               function spaces and maps between the PDE and the PIE
%               (space_pde, pie2pde, pde2pie) and PIE_in.misc are copied;
%
% NOTES
% - The result is the same PIE as returned by `convert_PIETOOLS_PDE` with
%   the flag 'nonlinear', for the same PDE.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie2npie_struct
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
% DJ, 10/15/2026: Initial coding;
% DJ, 10/29/2026: Update to new PIE structures: accept a PIE with operators
%                   in 'copvar' format (converting using `pie2opvar`), and
%                   return an 'npie_struct' object, instead of a struct;
% DJ, 10/30/2026: Convert all operators, the polynomials gz, gy, and
%                   the spaces and maps between the PDE and the PIE. Rename
%                   from PIE2polyPIE;


if ~PIE_in.is_initialized
    PIE_in = initialize(PIE_in);
end
obj_types = {'x','w','u','y','z'};
vnames = PIE_in.space.vars;     
space_pde = PIE_in.space_pde;

% Declare the function spaces, naming the state components using zero-
% padded numbers, so that the names are sorted in the order of the
% components.
PIE_out = npie_struct();
PIE_out.space.vars = vnames;        
PIE_out.space.dom = PIE_in.space.dom;
for obj = obj_types
    PIE_out.space.(obj{1}).dims = PIE_in.space.(obj{1}).dims;
    PIE_out.space.(obj{1}).depmat = PIE_in.space.(obj{1}).depmat;
end
ncomps = numel(PIE_in.space.x.dims);
PIE_out.space.x.varname = arrayfun(@(i) sprintf(['x%0',num2str(numel(num2str(ncomps))),'d_pie'],i),(1:ncomps)','UniformOutput',false);
PIE_out.space_pde.vars = space_pde.vars;      
PIE_out.space_pde.dom = space_pde.dom;
for obj = obj_types
    for f = {'dims','depmat','ID'}
        PIE_out.space_pde.(obj{1}).(f{1}) = space_pde.(obj{1}).(f{1});
    end
end

% Describe the components of each type of space by their sizes (dims) and
% the variables they depend on (depmat), to match them with the spaces of
% the operators.
space_pie = PIE_in.space;      
space_pie.k = struct('dims',[],'depmat',[]);
cop_convert = @(P,so,vo,si,vi) copvar2tensopmat(P,struct('out',so.dims,'in',si.dims),...
        struct('out',{vo},'in',{vi}),struct('out',so.depmat,'in',si.depmat));

% Convert the operators to 'tensopmat' objects
ops = {'T','T','x','x';     'Tw','Tw','x','w';      'Tu','Tu','x','u';
       'Bw','Bw','x','w';   'Bu','Bu','x','u';      'Dzw','Dzw','z','w';
       'Dzu','Dzu','z','u'; 'Dyw','Dyw','y','w';    'Dyu','Dyu','y','u';
       'K','K','k','x';      'Kw','Kw','k','w';      'Ku','Ku','k','u'};
for k=1:size(ops,1)
    P = PIE_in.(ops{k,2});
    if size(P,1)>0 && size(P,2)>0
        PIE_out.(ops{k,1}) = cop_convert(P,space_pie.(ops{k,3}),vnames,space_pie.(ops{k,4}),vnames);
    end
end
% Convert A, Cz, and Cy to polynomials in the state
x = PIE_out.x;
polys = {'f','A','x';   'gz','Cz','z';  'gy','Cy','y'};
for k=1:size(polys,1)
    P = PIE_in.(polys{k,2});
    if size(P,1)>0 && size(P,2)>0
        PIE_out.(polys{k,1}) = cop_convert(P,space_pie.(polys{k,3}),vnames,space_pie.x,vnames)*x;
    end
end

% Copy the maps between the PDE and the PIE, converting the operators
for obj = obj_types
    p2p = PIE_in.pde2pie.(obj{1});
    for f = setdiff(fieldnames(p2p)','C','stable')
        PIE_out.pde2pie.(obj{1}).(f{1}) = p2p.(f{1});
    end
    if isfield(p2p,'C') && isa(p2p.C,'copvar')
        PIE_out.pde2pie.(obj{1}).C = cop_convert(p2p.C,space_pie.(obj{1}),vnames,space_pie.k,[]);
    end
    if isfield(PIE_in.pie2pde,obj{1}) && isa(PIE_in.pie2pde.(obj{1}),'copvar')
        PIE_out.pie2pde.(obj{1}) = cop_convert(PIE_in.pie2pde.(obj{1}),space_pde.(obj{1}),space_pde.vars,space_pie.(obj{1}),vnames);
    end
end
PIE_out.misc = PIE_in.misc;
PIE_out = initialize(PIE_out);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
