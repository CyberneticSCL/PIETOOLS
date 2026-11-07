function PDE = static_constraints(PDE,mode)
% PDE = STATIC_CONSTRAINTS(PDE,MODE) returns a copy of the initialized PDE
% for displaying the static constraints 0=E(s) of the PDE.
%
% A static constraint is an equation that holds on the whole domain, such as
% the continuity equation 0 = d/ds1 v2 + d/ds2 v3, and is not a boundary
% condition. 'initialize' stores each one as an output z, flagged by its ID
% in PDE.constraint_ID, so that the converter can derive its operator. The
% display routines, however, print outputs as z = ..., so they get a copy of
% the PDE from this function, in which the constraints are no longer outputs.
%
% INPUT
% - PDE:    initialized 'pde_struct' object;
% - mode:   (optional) char:
%           'list' (default), the constraints are no outputs, and are
%                   stored in the field PDE.BC of the copy, only because
%                   this field is printed as equations 0 = ...;
%           'only', as 'list', but PDE.BC of the copy contains only the
%                   constraints, so that they can be summarized apart
%                   from the boundary conditions;
%           'drop', the constraints are removed from the outputs, and not
%                   listed anywhere;
%
% OUTPUT
% - PDE:    copy of the PDE, to be used for display only. If the PDE has no
%           static constraints, it is unchanged (for mode 'only', without
%           boundary conditions).
%
% NOTES
% - The static constraints of the actual PDE remain outputs flagged in
%   PDE.constraint_ID: they are not boundary conditions, and are not
%   handled as such by the converter.
%
% DB, 11/07/2026: Initial coding;

if nargin<2
    mode = 'list';
end
if isempty(PDE.constraint_ID) || isempty(PDE.z)
    if strcmp(mode,'only')
        PDE.BC = cell(0,1);     PDE.BC_tab = zeros(0,size(PDE.BC_tab,2));
    end
    return
end
z_IDs = cellfun(@(zi) zi.ID,PDE.z);
is_con = ismember(z_IDs(:),PDE.constraint_ID(:));
if ~any(is_con)
    if strcmp(mode,'only')
        PDE.BC = cell(0,1);     PDE.BC_tab = zeros(0,size(PDE.BC_tab,2));
    end
    return
end
if ~strcmp(mode,'drop')
    % Store the constraints in the field that the display routines print
    % as equations 0 = ...
    shown = PDE.BC(:);
    shown_tab = PDE.BC_tab;
    if strcmp(mode,'only')
        shown = cell(0,1);      shown_tab = zeros(0,size(PDE.BC_tab,2));
    end
    for ii=find(is_con)'
        con_ii = rmfield(PDE.z{ii},'ID');
        shown{end+1,1} = con_ii;
        row = zeros(1,size(shown_tab,2));
        row(1:size(PDE.z_tab,2)) = PDE.z_tab(ii,:);
        row(1) = size(shown_tab,1)+1;
        shown_tab(end+1,:) = row;
    end
    PDE.BC = shown;     PDE.BC_tab = shown_tab;
end
PDE.z(is_con) = [];
PDE.z_tab(is_con,:) = [];
PDE.z_tab(:,1) = (1:size(PDE.z_tab,1))';
PDE.constraint_ID = [];
end
