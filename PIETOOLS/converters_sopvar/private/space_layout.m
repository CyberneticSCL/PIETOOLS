function [masks,dims,comp_idx] = space_layout(sizes,dep,nvars)
% [MASKS,DIMS,COMP_IDX] = SPACE_LAYOUT(SIZES,DEP,NVARS) groups components by
% the space they exist in.
%
% INPUTS
% - sizes:  ncomps x 1 array with the size of each component;
% - dep:    ncomps x nvars logical array, indicating on which variables each
%           component depends;
% - nvars:  number of spatial variables;
%
% OUTPUTS
% - masks:  1 x nspaces array of the binary indices of the spaces that
%           contain a component, in increasing order, where the index of a
%           space is sum_j 2^(j-1)*[space contains variable j];
% - dims:   1 x nspaces array with the total size of the components in each space;
% - comp_idx:   ncomps x 1 cell, with for each component the indices in the
%           space it lives in that are associated to this component.
ncomps = numel(sizes);
if nvars>0 && ncomps>0
    cmask = double(dep)*(2.^(0:nvars-1)');
else
    cmask = zeros(ncomps,1);
end
masks = unique(cmask(:))';
dims = zeros(size(masks));
comp_idx = cell(ncomps,1);
for ii=1:numel(masks)
    for comp=find(cmask==masks(ii))'
        comp_idx{comp} = dims(ii) + (1:sizes(comp));
        dims(ii) = dims(ii) + sizes(comp);
    end
end
end
