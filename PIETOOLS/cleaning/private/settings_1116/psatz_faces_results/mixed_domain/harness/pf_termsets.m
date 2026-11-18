function [optl,names] = pf_termsets(which)
% [OPTL,NAMES] = PF_TERMSETS(WHICH) the term-set options of the term run as
% LPI_INEQ_SOP option structs (2-D registry): 'none' (plain term only),
% 'auto' (the reader's default), 'facesS' (faces, scalar offset), 'faces0'
% (faces at offset 0), 'faces' (faces, per-direction offset 1), 'box'
% (Putinar box factors, per direction), 'corners', 'faces_s1', 'faces_s2'
% (the faces of one direction), 'box_s1', 'box_s2', 'box_s1+faces' (the 1-D
% space's box factor with the faces). WHICH a cellstr of names (default all).
% SCRATCH, MMP, 11/16/2026.
all_names = {'none','auto','facesS','faces0','faces','box','corners','faces_s1','faces_s2','box_s1','box_s2','box_s1+faces'};
if nargin<1 || isempty(which),  which = all_names;   end
optl = cell(1,numel(which));    names = which;
for i = 1:numel(which)
    nm = which{i};
    switch nm
        case 'none',    optl{i} = struct('psatz','none');
        case 'auto',    optl{i} = struct();
        case 'faces_s1',        optl{i} = rowopt([1 0;2 0],[1 1],true);
        case 'faces_s2',        optl{i} = rowopt([0 1;0 2],[1 1],true);
        case 'box_s1',          optl{i} = rowopt([3 0],1,true);
        case 'box_s2',          optl{i} = rowopt([0 3],1,true);
        case 'box_s1+faces',    optl{i} = rowopt([3 0;1 0;2 0;0 1;0 2],[1 1 1 1 1],true);
        otherwise
            [rows,offs,pd] = pf_terms(nm,2,[]);
            optl{i} = rowopt(rows,offs,pd);
    end
end
end

function o = rowopt(rows,offs,pd)
o = struct('psatz',0,'psatz_offset',0,'psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd);
end
