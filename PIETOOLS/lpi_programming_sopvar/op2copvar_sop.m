function Pm = op2copvar_sop(Pop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PM = OP2COPVAR_SOP(POP) converts a legacy PI operator to a container:
%   opvar      -> opvar2copvar (Testfolder/converters);
%   opvar2d    -> a copvar over the present spaces of POP (R^n, L2[x],
%                 L2[y], L2[x,y]), one block per nonzero component, with a
%                 zero block kept in every all-zero row and column so that
%                 no space is dropped (opvar2d2copvar errors there), as the
%                 2-D transcriptions do;
%   copvar/cdopvar -> returned as is.
%
% Initial coding MMP, 10/08/2026 (the opvar2d branch is the logic of the
% frozen test helper cx_hinf_op2d, copied so that the executives depend on
% the library only).
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
        Pblk.(nm{i,j}) = Pop.(nm{i,j});
        C{a,b} = opvar2d2sopvar(Pblk);
    end
end
Pm = copvar(C);
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
