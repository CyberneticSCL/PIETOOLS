function Pm = op2copvar_sop(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM = OP2COPVAR_SOP(POP) converts a legacy PI operator to a container:
%   opvar      -> opvar2copvar (Testfolder/converters);                     % MMP, 10/16/2026 (was)
%   opvar      -> opvar2copvar (converters/operator_converters);            % MMP, 10/16/2026
%   opvar2d    -> a copvar over the present spaces of POP (R^n, L2[x],
%                 L2[y], L2[x,y]), one block per nonzero component, with a
%                 zero block kept in every all-zero row and column so that
%                 no space is dropped (opvar2d2copvar errors there), as the
%                 2-D transcriptions do; with no space at all on one side,  % MMP, 10/28/2026
%                 a container with no space on that side;                   % MMP, 10/28/2026
%   copvar/cdopvar -> returned as is.
%
% Initial coding MMP, 10/08/2026 (the opvar2d branch is the logic of the
% frozen test helper cx_hinf_op2d, copied so that the executives depend on
% the library only).
% MMP, 10/16/2026: moved from lpi_programming_sopvar/ to converters/
%                operator_converters/, with the operator converters of
%                sopvar/Testfolder/converters: a conversion between operator
%                classes, not LPI programming, and the executives must not
%                depend on a test folder. Comment only otherwise.
% MMP, 10/28/2026: An opvar2d with no space on one side (an empty channel of
%                a PIE, such as C1 with no regulated output) gives a
%                container with no space on that side, its spaces stated
%                from POP.dim; before, the side with no rows or columns made
%                the zero-block loops fail (maintainer: empty channels
%                should not raise an error; opvar2copvar does the same in
%                1-D).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(Pop,'copvar') || isa(Pop,'cdopvar'),     Pm = Pop;   return,     end
if isa(Pop,'opvar') || isa(Pop,'dopvar')
    Pm = opvar2copvar(Pop);     return
end
if ~isa(Pop,'opvar2d')
    error('op2copvar_sop:class','POP should be an opvar, opvar2d or a container.');
end
nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
d = Pop.dim;
rows = find(d(:,1)>0);  cols = find(d(:,2)>0);
keep = false(4,4);
for i = rows',  for j = cols',  keep(i,j) = ~is_zero(Pop.(nm{i,j}));  end,  end
% for i = rows',  if ~any(keep(i,cols)),  keep(i,cols(1)) = true;   end,    end % MMP, 10/28/2026 (was)
% for j = cols',  if ~any(keep(rows,j)),  keep(rows(1),j) = true;   end,    end % MMP, 10/28/2026 (was)
% A side with no space has no block to place (see the end)                  % MMP, 10/28/2026
for i = rows',  if ~isempty(cols) && ~any(keep(i,cols)),  keep(i,cols(1)) = true;   end,    end % MMP, 10/28/2026
for j = cols',  if ~isempty(rows) && ~any(keep(rows,j)),  keep(rows(1),j) = true;   end,    end % MMP, 10/28/2026
C = cell(numel(rows),numel(cols));
for a = 1:numel(rows)
    for b = 1:numel(cols)
        i = rows(a);    j = cols(b);
        if ~keep(i,j),  continue,   end
        Pblk = opvar2d();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        sel = zeros(4,2);   sel(i,1) = d(i,1);  sel(j,2) = d(j,2);
        Pblk.dim = sel;
        Pblk.(nm{i,j}) = Pop.(nm{i,j});
        C{a,b} = opvar2d2sopvar(Pblk);
    end
end
% Pm = copvar(C);                                                           % MMP, 10/28/2026 (was)
if ~isempty(rows) && ~isempty(cols)                                         % MMP, 10/28/2026
    Pm = copvar(C);                                                         % MMP, 10/28/2026
else                                                                        % MMP, 10/28/2026
    % No space on one side (an empty channel of a PIE): no block to read    % MMP, 10/28/2026
    % the spaces off, so they are stated from POP.dim. Space k of R, L2[x], % MMP, 10/28/2026
    % L2[y], L2[x,y] holds the variables of row k of MASK; the registry is  % MMP, 10/28/2026
    % those the present spaces hold, sorted, with their domains             % MMP, 10/28/2026
    nmv = reshape(Pop.var1.varname,1,[]);                                   % MMP, 10/28/2026
    mask = logical([0 0; 1 0; 0 1; 1 1]);                                   % MMP, 10/28/2026
    used = any(mask(unique([rows;cols]),:),1);                              % MMP, 10/28/2026
    [vs,ix] = sort(nmv(used));                                              % MMP, 10/28/2026
    dI = Pop.I(used,:);         mk = mask(:,used);      mk = mk(:,ix);      % MMP, 10/28/2026
    meta = struct('vars',{vs},'dom',dI(ix,:),'space_out',mk(rows,:),...
                  'space_in',mk(cols,:),'dim_out',d(rows,1),'dim_in',d(cols,2)); % MMP, 10/28/2026
    Pm = copvar(C,meta);                                                    % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
end


function tf = is_zero(comp)
if iscell(comp)
    tf = true;
    for k = 1:numel(comp),  tf = tf && is_zero(comp{k});  end
    return
end
if isempty(comp),   tf = true;  return,     end
comp = polynomial(comp);
tf = isempty(comp.coefficient) || ~any(comp.coefficient(:));
end
