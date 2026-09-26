function [prog,P] = cx_h2_pos2d(prog,sp,dm,dom,blk,terms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,P] = CX_H2_POS2D(PROG,SP,DM,DOM,BLK,TERMS) is the container call
% standing in for a 2-D executive's
%
%   [prog,P] = poslpivar_2d(prog,dim,deg,opts);
%   for j: if use_psatz(j), [prog,P2] = poslpivar_2d(prog,dim,deg_j,opts_j);
%          P = P + P2; end
%
% on a state held in L2[x,y] ONLY. BLK and TERMS are one block and the
% psatz list of settings2possopvar's output (out.LF, out.LF_psatz or out.eq,
% out.eq_psatz): degrees with the full subset-cap array, include (from
% exclude) and per-direction sep. That translation covers the L2[x,y]
% block alone - possopvar is L2 -> L2 - so any other space is an error
% here rather than a guessed degree mapping.
%
% poscopvar takes the same vocabulary for one space holding every registry
% variable (copquadvar.m:441-447), so BLK.deg is passed as that space's
% per-basis-operator cell and BLK.opts.include as its include entry.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if numel(sp)~=1 || numel(sp{1})~=2
    error('cx_h2_pos2d:spaces',['Only an L2[x,y]-only operator is transcribed: '...
          'settings2possopvar translates no R, L2[x] or L2[y] degrees.'])
end
[prog,P] = decl(prog,sp,dm,dom,blk);
for j = 1:numel(terms)
    [prog,P2] = decl(prog,sp,dm,dom,terms{j});
    P = P + P2;
end
end

function [prog,P] = decl(prog,sp,dm,dom,b)
opt = struct('include',{{b.opts.include}},'sep',b.opts.sep,'psatz',b.opts.psatz);
[prog,P] = poscopvar(prog,dm,sp,dom,{b.deg},opt);
end
