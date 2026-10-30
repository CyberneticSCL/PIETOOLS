function PIE = pie_pde_maps(PIE,PDE0,del_tab,Wop,Uop,Fop)
% PIE = PIE_PDE_MAPS(PIE,PDE0,DEL_TAB,WOP,UOP,FOP) declares the function
% spaces of a PDE, and the maps between this PDE and the PIE constructed
% from it, in the fields PIE.space_pde, PIE.pie2pde and PIE.pde2pie.
%
% INPUT
% - PIE:    'pie_struct' or 'npie_struct' object, with the operators and
%           the components of PIE.space specified, and with the tables
%           x_tab, w_tab, u_tab, y_tab, and z_tab set, as obtained from a
%           conversion of PDE0, so that the first column of each table
%           refers to the components of the PDE after the delays have been
%           expanded. For the state, PIE.pde2pie.x.tdiff must have been set
%           if the PIE involves higher-order temporal derivatives;
% - PDE0:   Initialized 'pde_struct' object, the PDE to which the PIE is
%           related, before expanding any delays;
% - del_tab:    (optional) table returned by `expand_delays`, relating the
%           components of the PDE added to represent delays to the
%           components of PDE0. Empty if there are no delays;
% - Wop, Uop:   (optional) operators mapping the PIE inputs w and u to the
%           PDE inputs, as 'opvar' or 'opvar2d' objects. The default 1
%           indicates that the PDE inputs are the same as the PIE inputs;
% - Fop:    (optional) 'opvar' or 'opvar2d' object mapping the PDE state to
%           the first PIE.misc.periodic.n_aux (moment) components of the
%           PIE state, for PIEs for PDEs with generalized periodic BCs.
%           Linear PIEs only;
%
% OUTPUT
% - PIE:    The same PIE, with the fields
%       space_pde:  function spaces of the PDE0, with the components of
%                   each of x, w, u, y, and z identified by
%                   PDE0.(.){i}.ID (and, for an 'npie_struct', unnamed:
%                   the names are assigned by initialize);
%       pie2pde:    operators X_pde, W_pde, U_pde, Y_pde, and Z_pde,
%                   'copvar' objects for a 'pie_struct', and 'tensopmat'
%                   objects for an 'npie_struct';
%       pde2pie:    for each of the fields x, w, u, y, and z, the PDE
%                   components idx corresponding to each of the PIE
%                   components, and for x, w, and u also the fields diff,
%                   tdiff, loc, and C, and for x the field delay, as
%                   described in `pie_struct`;
%
% NOTES
% - The component i of the PDE state is represented in the PIE by the
%   state component j for which T maps to this component, with tdiff(j)=0,
%   delay(j)=0, and idx(j,1)=i. The operator X_pde has an identity block
%   at position (i,j) for each such pair, and is zero otherwise.
% - Unless the PIE is constructed for a PDE with periodic BCs, the
%   operator C is the identity operator on the PIE components. For the
%   moment components of a PIE for periodic BCs, C is instead the
%   associated row of Fop, mapping the PDE component to a real number.
% - If Wop or Uop is not 1 (i.e. if the PIE inputs include boundary values
%   or derivatives of the PDE inputs), the associated operator W_pde or
%   U_pde has the rows and columns of the PDE inputs ordered and grouped by
%   the type of space they exist in, as for an 'opvar2d' object.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_pde_maps
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
% DJ (CC), 10/29/2026: Initial coding;

if nargin<3,    del_tab = [];   end
if nargin<4,    Wop = 1;        end
if nargin<5,    Uop = 1;        end
if nargin<6,    Fop = [];       end
is_lin = isa(PIE,'pie_struct');
vnames = PIE.space.vars;        dom = PIE.space.dom;
N = numel(vnames);
types = {'x','w','u','y','z'};
n_aux = 0;
if ~isempty(Fop),   n_aux = PIE.misc.periodic.n_aux;    end


% % % The function spaces of the PDE, and their dependence on the variables
% % % of the PIE. If the variables were merged and renamed (by `combine_vars`,
% % % in expanding delays), the k-th PDE variable is the k-th PIE variable.
N0 = size(PDE0.vars,1);
vnames0 = cell(1,N0);
for ii=1:N0
    vnames0{ii} = PDE0.vars(ii,1).varname{1};
end
[~,loc0] = ismember(vnames0,vnames);
if any(loc0==0),    loc0 = 1:N0;    end
sp_pde = struct('vars',{vnames0},'dom',reshape(double(PDE0.dom),N0,2));
dep_pde = struct();
for t = types
    tab = PDE0.([t{1},'_tab']);
    sp_pde.(t{1}).dims = tab(:,2);
    sp_pde.(t{1}).depmat = logical(tab(:,3:2+N0));
    sp_pde.(t{1}).ID = zeros(size(tab,1),1);
    if ~is_lin,     sp_pde.(t{1}).varname = cell(0,1);      end     % names, assigned by initialize
    for ii=1:size(tab,1)
        sp_pde.(t{1}).ID(ii) = PDE0.(t{1}){ii}.ID;
    end
    dep_pde.(t{1}) = false(size(tab,1),N);
    dep_pde.(t{1})(:,loc0) = sp_pde.(t{1}).depmat;
end
PIE.space_pde = sp_pde;


% % % The maps from the PDE to the PIE, other than C
p2p = PIE.pde2pie;
nx = size(p2p.x.idx,1);
idx = zeros(nx,3);          loc = zeros(nx,N);          delay = zeros(nx,1);
diff = p2p.x.diff;          tdiff = p2p.x.tdiff;
if isempty(tdiff),          tdiff = zeros(nx,1);        end
nx0 = numel(PDE0.x);        nw0 = numel(PDE0.w);
for ii=1:nx
    p = p2p.x.idx(ii,1);
    if p<=nx0
        idx(ii,1) = p;
        continue
    end
    % Components added to represent delays refer to the PDE component of
    % which they are the history, see the table returned by expand_delays.
    if isempty(del_tab)
        error("Components of the PDE not in the original PDE should represent delays.")
    end
    r = p-nx0;
    k = del_tab(r,1);
    if k<=nx0,          idx(ii,1) = k;
    elseif k<=nx0+nw0,  idx(ii,2) = k-nx0;
    else,               idx(ii,3) = k-nx0-nw0;
    end
    % The delay variable is the variable of the history that the PDE
    % component does not depend on (or the merged variable, if merged).
    col = find(idx(ii,:));
    dv = find(PIE.space.x.depmat(ii,:));
    rest = dv(~dep_pde.(types{col})(idx(ii,col),dv));
    if isempty(rest),   rest = dv;  end
    if isempty(rest)
        error("A component representing a delay should depend on a spatial variable.")
    end
    delay(ii) = rest(1);
    % The history x_i(t,s_l) of the PDE component is differentiated wrt the
    % spatial variable s_l, which is the same as a temporal derivative.
    tdiff(ii) = tdiff(ii) + diff(ii,delay(ii));
    diff(ii,:) = 0;
    loc(ii,loc0) = del_tab(r,2+(1:N0));
    diff(ii,loc0) = del_tab(r,2+size(del_tab,2)/2-1+(1:N0));
end
p2p.x.idx = idx;            p2p.x.loc = loc;        p2p.x.delay = delay;
p2p.x.diff = diff;          p2p.x.tdiff = tdiff;
% Inputs: boundary values of an input are evaluated at the lower boundary
% of the variables on which the PIE input no longer depends.
for t = {'w','u'}
    n = size(p2p.(t{1}).idx,1);
    dep_pie = double(PIE.space.(t{1}).depmat);
    if isempty(dep_pie),    dep_pie = zeros(n,N);   end
    p2p.(t{1}).loc = -double(dep_pde.(t{1})(p2p.(t{1}).idx(:,1),:) & ~dep_pie);
    p2p.(t{1}).tdiff = zeros(n,1);
    if isempty(p2p.(t{1}).diff),    p2p.(t{1}).diff = zeros(n,N);   end
end


% % % The operators pie2pde (k=1), mapping the PIE components to the PDE
% % % components, and C (k=2), mapping the PDE components to the PIE components,
% % % are identity operators between the associated components, so they are
% % % defined by a matrix with identity blocks at (rows(j),cols(j)).
for t = types
    n = size(p2p.(t{1}).idx,1);
    if n==0,    continue,   end
    for k=1:2
        if k==2 && any(strcmp(t{1},{'y','z'})),     continue,   end
        dims_pie = PIE.space.(t{1}).dims(:);
        dep_pie = logical(PIE.space.(t{1}).depmat);
        if k==1
            % PDE components i=idx(j) from the PIE components j. For the
            % state, only the components that are neither derivatives in
            % time, nor histories, nor moments, map to the PDE state.
            sel = true(n,1);
            if strcmp(t{1},'x'),    sel = tdiff==0 & delay==0 & idx(:,1)>0;  sel(1:n_aux) = false;  end
            rows = p2p.(t{1}).idx(sel,1);       cols = find(sel);
            dims_out = sp_pde.(t{1}).dims;      dep_out = dep_pde.(t{1});
            dims_in = dims_pie;                 dep_in = dep_pie;
            Wmap = 1;
            if strcmp(t{1},'w'),    Wmap = Wop;     end
            if strcmp(t{1},'u'),    Wmap = Uop;     end
            if ~isequal(Wmap,1)
                % Inputs include derivatives or boundary values, so the map
                % is not an identity, but an operator from W/U.
                if is_lin,  PIE.pie2pde.(t{1}) = op2copvar_sop(Wmap);
                else,       PIE.pie2pde.(t{1}) = tensopmat(Wmap);
                end
                continue
            end
        else
            % PIE components j from the PDE components, after taking the
            % derivatives and evaluating. The moments (F*x) map a PDE
            % component, so these have the size and dependence of the PDE.
            rows = (1:n)';      cols = rows;
            dims_out = dims_pie;                dep_out = dep_pie;
            dims_in = dims_pie;                 dep_in = dep_pie;
            if strcmp(t{1},'x') && n_aux>0
                rows = (n_aux+1:n)';        cols = rows;
                dims_in(1:n_aux) = sp_pde.x.dims(idx(1:n_aux,1));
                dep_in(1:n_aux,:) = dep_pde.x(idx(1:n_aux,1),:);
            end
        end
        % Sparse matrix with identity blocks.
        off_out = [0;cumsum(dims_out(:))];      off_in = [0;cumsum(dims_in(:))];
        II = cell(numel(rows),1);               JJ = cell(numel(rows),1);
        for jj=1:numel(rows)
            m = dims_in(cols(jj));
            if dims_out(rows(jj))~=m
                error("Components of the PIE and the PDE are of different sizes.")
            end
            II{jj} = off_out(rows(jj)) + (1:m)';
            JJ{jj} = off_in(cols(jj)) + (1:m)';
        end
        Mat = sparse(cell2mat(II),cell2mat(JJ),1,off_out(end),off_in(end));
        [vs,perm] = sort(vnames);
        dims = struct('out',dims_out(:),'in',dims_in(:));
        if is_lin
            spaces = struct('out',{cell(1,numel(dims_out))},'in',{cell(1,numel(dims_in))});
            for ii=1:numel(dims_out),   spaces.out{ii} = vnames(dep_out(ii,:));     end
            for ii=1:numel(dims_in),    spaces.in{ii} = vnames(dep_in(ii,:));       end
            Pop = mat2copvar_sop(Mat,dims,spaces,struct('vars',{vs},'dom',dom(perm,:)));
        else
            Pop = mat2tensopmat(Mat,dims,struct('in',{vs},'out',{vs}),...
                struct('in',dom(perm,:),'out',dom(perm,:)),...
                struct('in',double(dep_in(:,perm)),'out',double(dep_out(:,perm))));
        end
        if k==1
            PIE.pie2pde.(t{1}) = Pop;
        else
            % Replace the identity blocks of the moments by the associated
            % rows of Fop, mapping the PDE component to a real number.
            if strcmp(t{1},'x') && n_aux>0
                % Fop has rows for the moment components, and columns for the
                % PDE state components, both ordered by type of space (R,
                % L2[s1], L2[s2], L2[s1,s2]), as in the PIE.
                Cblk = Pop.C;
                kinds = 1 + double(dep_pie(n_aux+1:end,:))*(1:N)';
                for ii=1:n_aux
                    q = find(idx(n_aux+1:end,1)==idx(ii,1),1);
                    rows = sum(dims_pie(1:ii-1)) + (1:dims_pie(ii));
                    cols = sum(Fop.dim(1:kinds(q)-1,2)) + sum(dims_pie(n_aux+find(kinds(1:q-1)==kinds(q)))) + (1:dims_in(ii));
                    Fblk = op2copvar_sop(Fop(rows,cols));
                    Cblk{ii,ii} = Fblk.C{1,1};
                end
                Pop = copvar(Cblk,metadata(Pop));
            end
            p2p.(t{1}).C = Pop;
        end
    end
end
PIE.pde2pie = p2p;

end
