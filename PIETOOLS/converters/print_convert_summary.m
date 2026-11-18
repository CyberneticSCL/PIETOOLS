function print_convert_summary(PIE,obj,aux_idcs,F_aux)
% PRINT_CONVERT_SUMMARY(PIE,obj,aux_idcs,F_aux) prints how the components
% of the state, inputs or outputs of the PIE (PIE.space) relate to those
% of the PDE (PIE.space_pde), as described by PIE.pde2pie.
%
% INPUT
% - PIE:    'pie_struct' or 'npie_struct' object;
% - obj:    'x', 'u', 'w', 'y' or 'z', the object to display;
% - aux_idcs, F_aux:  (optional) the components of the PIE state that are
%           auxiliary states xf=F*x, displayed as F*x, and their F;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - print_convert_summary
%
% Copyright (C)2026 PIETOOLS Team
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
% DB, 11/18/2026: Move from local function of convert_PIETOOLS_PDE to its own
%                   file, to also display the conversion of a BVP;

names = struct('x','fundamental state component','y','observed output','z','regulated output',...
               'u','actuator input','w','exogenous input');
object_name = names.(obj);
sp = PIE.space.(obj);       p2p = PIE.pde2pie.(obj);
ncomps = numel(sp.dims);    nvars = numel(PIE.space.vars);
% Derivatives in space (diff) and time (tdiff), position at a boundary
% (loc: -1 lower, 1 upper) and delay variable (delay) of the PDE component
% that each component is; zero if not specified
maps = struct('diff',zeros(ncomps,nvars),'tdiff',zeros(ncomps,1),'loc',zeros(ncomps,nvars),'delay',zeros(ncomps,1));
for fn=fieldnames(maps)'
    if isfield(p2p,fn{1}) && size(p2p.(fn{1}),1)==ncomps
        maps.(fn{1}) = p2p.(fn{1});
    end
end
% The PDE component of each component: of type src_type and index src (a
% component of the PIE state may be the history of a delayed input)
src_type = repmat({obj},ncomps,1);
[src,col] = max(p2p.idx,[],2);
if strcmp(obj,'x')
    xwu = {'x','w','u'};
    src_type = xwu(col)';
end
% The variables of the PDE in those of the PIE (the same, if renamed)
[~,vmap] = ismember(PIE.space_pde.vars,PIE.space.vars);
if any(vmap==0)
    vmap = 1:numel(PIE.space_pde.vars);
end
has_sdiff = any(maps.diff(:));
has_diffs = has_sdiff || any(maps.tdiff) || any(maps.delay) || ~all(strcmp(src_type,obj));
if ~has_diffs && all(src==(1:ncomps)')
    % The components in the PIE correspond exactly to those in the PDE
    return
elseif ~has_diffs
    % Otherwise, we list the new order of the components...
    fprintf(['\n','The ',object_name,'s have been reindexed as:\n']);
else
    % ...or how the new components relate to those in the PDE
    fprintf(['\n','The following ',object_name,'s have been introduced:\n']);
end

% Use UNICODE for the partial derivatives, and for the superscript orders
% of the derivatives (sup_num{k+1} for order k)
thin_space = char(8201);
partial = '\x2202';
sup_num = {'\x2070';'\xB9';'\xB2';'\xB3';'\x2074';'\x2075';'\x2076';'\x2077';'\x2078';'\x2079'};
% The PIE state and inputs are indicated as such in their names (a hat),
% unlike those of the PDE; the names are given by 'varname2display'
pie_tag = '';
if any(strcmp(obj,{'x','u','w'}))
    pie_tag = '_pie';
end
% Determine how many subscripts will be needed.
n_digits = 1+floor(log(ncomps)/log(10));


% % For a clearer display, estimate how long the display will be: the
% % names of the variables, and the length of the component names
varname_cell = {'s'};
if nvars>1
    varname_cell = arrayfun(@(kk) varname2display(sprintf('s%d',kk)),(1:nvars)','UniformOutput',false);
end
lngth_varnames = 3;      % just three characters: (t)
if any(sp.depmat(:)) || has_sdiff
    lngth_varnames = lngth_varnames+2*nvars;    % add ,sj for each variable
end
nvars_max = max(sum(sp.depmat,2));
lngth_varnames_mean = ceil(lngth_varnames*nvars_max/nvars);
LHS_length_max = 1+1 + n_digits + lngth_varnames_mean+3; % e.g. ' x13(t,s1,s2,s3)',


% % For each of the components, display its size, and which variables it
% % depends on.
for ii=1:ncomps
    % Establish the names of the variables on which the component depends.
    varnames_LHS_t = 't';
    LHS_length = 1;
    for kk=find(sp.depmat(ii,:))
        varnames_LHS_t = [varnames_LHS_t,',',varname_cell{kk}];
        LHS_length = LHS_length+3;  % add 3 characters: ,s1
    end
    % The PDE component, on its own variables, at time t+s_k for a delay
    % in s_k, and at the boundary loc (lower by default) in the variables
    % that the component of the PIE does not depend on
    varnames_RHS_t = 't';
    if maps.delay(ii)>0
        varnames_RHS_t = ['t+',varname_cell{maps.delay(ii)}];
    end
    dep_src = false(1,nvars);
    dep_src(vmap) = PIE.space_pde.(src_type{ii}).depmat(src(ii),:);
    for kk=find(dep_src)
        if sp.depmat(ii,kk) && maps.loc(ii,kk)==0
            varnames_RHS_t = [varnames_RHS_t,',',varname_cell{kk}];
        else
            varnames_RHS_t = [varnames_RHS_t,',',num2str(PIE.space.dom(kk,1+(maps.loc(ii,kk)==1)))];
        end
    end

    % Establish the index for the new component, displayed as a subscript.
    Lcomp_idx = '';
    if ncomps>1
        Lcomp_idx = num2str(ii);
        LHS_length = LHS_length + length(Lcomp_idx);
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    LHS_name = [' ',varname2display([obj,pie_tag,Lcomp_idx]),'(',varnames_LHS_t,')'];
    LHS_length = 1 + LHS_length + 3;

    % Establish the index for the old component
    Rcomp_idx = '';
    if ncomps>1
        Rcomp_idx = num2str(src(ii));
    end
    % Set the name of the component, including its dependence on spatial
    % variables.
    RHS_name = [thin_space,varname2display([src_type{ii},Rcomp_idx]),'(',varnames_RHS_t,')'];
    % For state components, also indicate what spatial/temporal derivative
    % of the PDE state they correspond to.
    if strcmp(obj,'x') && any(maps.tdiff)
        tdiff = maps.tdiff(ii);
        if tdiff==0
            diff_str = ['       ',thin_space];
        elseif tdiff==1
            diff_str = [' (',partial,'/',partial,'t)',thin_space];
        else
            t_sup = cell2mat(sup_num(str2num(num2str(tdiff)')+1)');
            diff_str = ['(',partial,'/',partial,'t)',t_sup,thin_space];
        end
    else
        diff_str = '';
    end
    if has_sdiff
        for jj=1:nvars
            sdiff = maps.diff(ii,jj);
            if sdiff==0
                diff_str = ['        ',diff_str,thin_space];
            elseif sdiff==1
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',thin_space];
            else
                s_sup = cell2mat(sup_num(str2num(num2str(sdiff)')+1)');
                diff_str = [diff_str,'(',partial,'/',partial,varname_cell{jj},')',s_sup,thin_space];
            end
        end
    end
    if has_diffs
        RHS_name = [diff_str,RHS_name];
    end
    if nargin>=3 && ismember(ii,aux_idcs)
        % DB 11/07/2026: auxiliary state xf=F*x, displayed as the integral
        % form of F*x (not as a boundary value of x)
        F_str = 'F*x';
        try
            F_blk = F_aux{aux_idcs==ii};
            F_d = get_display(F_blk,varname2display(['x',Rcomp_idx]));
            if F_d.ok
                F_str = '';
                for kk=1:numel(F_d.terms)
                    if kk==1
                        sgn = repmat('-',1,F_d.signs(kk)<0);
                    else
                        sgn = [' ',char(43+2*(F_d.signs(kk)<0)),' '];
                    end
                    F_str = [F_str,sgn,F_d.terms{kk}];
                end
            end
        catch
        end
        RHS_name = [diff_str,thin_space,F_str];
    end
    % % % Finally, display:
    MT_space = max(LHS_length_max-LHS_length,1);
    fprintf(['  ',LHS_name,repmat(' ',[1,MT_space]),' <--  ',RHS_name,'\n']);
end
end
