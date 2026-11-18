function [optl,names] = pf_termsets3(which)
% [OPTL,NAMES] = PF_TERMSETS3(WHICH) the product term sets of the follow-up
% run T3 (maintainer, 11/16 11:30: the theory supports products of faces
% and products of box quadratics, and mixed-direction products are what a
% mixed domain may need), as LPI_INEQ_SOP option structs on the 2-D
% registry (rows of per-direction codes: 1 lower face, 2 upper face, 3 the
% box quadratic g_k = (s_k - a_k)(b_k - s_k)):
%   'product2'  {g1, g2, g1 g2}: the preordering of the box quadratics
%   'mixedfb'   the four faces and the face x box-quadratic products of
%               different directions, f_{1,i} g2 and g1 f_{2,j}
%   'full16'    every product of distinct faces (both faces of one
%               direction = its box quadratic): the 15 nonzero rows of
%               {0,1,2,3}^2, the Schmuedgen preordering of the four faces
% each at the per-direction offset 1 (the default rule of the proposal,
% eq. d4) and, with suffix 0, at offset 0. SCRATCH, MMP, 11/16/2026.
all_names = {'product2','product20','mixedfb','mixedfb0','full16','full160'};
if nargin<1 || isempty(which),  which = all_names;   end
optl = cell(1,numel(which));    names = which;
for i = 1:numel(which)
    nm = which{i};  o = 1;
    if endsWith(nm,'0'),    o = 0;  nm = nm(1:end-1);   end
    switch nm
        case 'product2',    rows = [3 0;0 3;3 3];
        case 'mixedfb',     rows = [1 0;2 0;0 1;0 2;1 3;2 3;3 1;3 2];
        case 'full16',      rows = [1 0;2 0;3 0;0 1;0 2;0 3;1 1;1 2;2 1;2 2;1 3;2 3;3 1;3 2;3 3];
        otherwise,          error('pf_termsets3:name','unknown term set %s',which{i});
    end
    optl{i} = struct('psatz',0,'psatz_offset',0,'psatz_rows',rows,'psatz_rows_offset',o*ones(1,size(rows,1)),'psatz_rows_perdir',true);
end
end
