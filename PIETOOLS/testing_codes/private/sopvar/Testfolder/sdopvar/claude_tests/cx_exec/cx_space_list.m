function [sp,dm] = cx_space_list(P,side)
% [SP,DM] = CX_SPACE_LIST(P,SIDE) lists a container's output ('out') or
% input ('in') spaces as poscopvar / lpivar_cdopvar take them: one cellstr of
% variable names per space ({} for R^n), and the component counts. Same as
% the helper of the same job in test_copvar_kyp.
%
% Initial coding MMP, 09/25/2026
if strcmp(side,'out'),  S = P.space_out;    dm = P.dim_out(:);
else,                   S = P.space_in;     dm = P.dim_in(:);
end
sp = cell(1,size(S,1));
for k = 1:size(S,1),    sp{k} = reshape(P.vars(S(k,:)),1,[]);    end
end
