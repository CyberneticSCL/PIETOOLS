function Pm = cx_hinf_op2d(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM = CX_HINF_OP2D(POP) is opvar2d2copvar(POP), except that a row or column
% left with no block by zero components gets one explicit zero block.
%
% opvar2d2copvar drops every zero component and then derives the metadata
% from the blocks, so it errors (copvar:emptyRow / emptyColumn) on e.g.
% D11 = 0 - a 1 x 1 grid whose only block is zero. opvar2copvar (1-D)
% received the explicit-zero-block fix on 09/25/2026; opvar2d2copvar did
% not. This is that fix, done locally rather than in the library: each
% component is still converted by opvar2d2sopvar on a single-component
% opvar2d, as opvar2d2copvar does, and the zero block is simply the
% converted zero component (degree-0 basis, no decision variables).
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nm = {'R00','R0x','R0y','R02';
      'Rx0','Rxx','Rxy','Rx2';
      'Ry0','Ryx','Ryy','Ry2';
      'R20','R2x','R2y','R22'};
d = Pop.dim;
rows = find(d(:,1)>0);  cols = find(d(:,2)>0);
keep = false(4,4);
for i = rows',  for j = cols',  keep(i,j) = ~is_zero(Pop.(nm{i,j}));  end,  end
for i = rows',  if ~any(keep(i,cols)),  keep(i,cols(1)) = true;   end,    end
for j = cols',  if ~any(keep(rows,j)),  keep(rows(1),j) = true;   end,    end
C = cell(numel(rows),numel(cols));
for a = 1:numel(rows)
    for b = 1:numel(cols)
        i = rows(a);    j = cols(b);
        if ~keep(i,j),  continue,   end
        Pblk = opvar2d();
        Pblk.I = Pop.I;     Pblk.var1 = Pop.var1;   Pblk.var2 = Pop.var2;
        sel = zeros(4,2);   sel(i,1) = d(i,1);  sel(j,2) = d(j,2);
        Pblk.dim = sel;
        Pblk.(nm{i,j}) = Pop.(nm{i,j});         % may be the zero component
        C{a,b} = opvar2d2sopvar(Pblk);
    end
end
Pm = copvar(C);
end

function tf = is_zero(comp)
% Same test as opvar2d2copvar's is_zero_component: cells are the alpha split.
if iscell(comp)
    tf = true;
    for k = 1:numel(comp),  tf = tf && is_zero(comp{k});  end
    return
end
if isempty(comp),   tf = true;  return,     end
comp = polynomial(comp);
tf = isempty(comp.coefficient) || ~any(comp.coefficient(:));
end
