function Test_mdopvar()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_MDOPVAR exercises the decision-variable container: the single-Zd
% invariant, promotion of a 'mopvar' operand, mixed block classes, and the
% refusal to compose two decision containers.
%
% SCOPE. The block algebra is not retested here - that is @sdopvar's, and it
% is covered by Testfolder/sdopvar. What is tested is what the CONTAINER is
% responsible for: that every block ends on one decision variable list and
% stays there under each operation, that dispatch against a 'mopvar' works
% in both argument positions, and that the block grid and spaces come out
% right. The operator semantics of the container - plus, adjoint and
% composition against the independently evaluated action - are tested in
% Test_mopvar, whose blocks are 'sopvar' and can therefore be applied to a
% polynomial by 'apply_sopvar'; the container logic is the same code path.
%
% SPACE STRUCTURE
% Spaces are R^n, L_2[s1] and L_2[s1,s2], with distinct component counts in
% every position so that a swapped block index cannot pass.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_mdopvar
%
% Copyright (C) 2026 PIETOOLS Team
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
% Initial coding MMP, 09/17/2026

warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');
np = 0;     nf = 0;
rng(0);

% Spaces: S{1} = R, S{2} = L2[s1], S{3} = L2[s1,s2]; all on [0,1].
S = {cell(1,0),{'s1'},{'s1','s2'}};
dm = [0 1];
out_D = [1 2];      din_D = [2 3];      % D: (L2[s1],L2[s1,s2]) -> (R,L2[s1])
qD = [2 3];         pD = [2 3];

fprintf('\n===== Test_mdopvar =====\n');

%% ---------------------------------------------------------------- 1
fprintf('\n-- 1. the single decision variable list --\n');
% Blocks built independently, so their Zd lists differ. The container must
% put them all on one list, and share the one array rather than copies.
D = cell(2,2);
for i = 1:2
    for j = 1:2
        vs = struct('out',{S{out_D(i)}},'in',{S{din_D(j)}});
        ds = struct('out',repmat(dm,numel(vs.out),1),'in',repmat(dm,numel(vs.in),1));
        gs = struct('out',ones(1,numel(vs.out)),'in',ones(1,numel(vs.in)));
        D{i,j} = rand_sdopvar([qD(i),pD(j)],vs,ds,gs,20+5*(i+j),0.4);
    end
end
[np,nf] = chk(np,nf,'blocks were built on different lists',...
    ~isequal(D{1,1}.Zd(:),D{2,2}.Zd(:)));
Dm = mdopvar(D);
[np,nf] = chk(np,nf,'constructor puts every block on one list',...
    all(cellfun(@(b) isequal(b.Zd,Dm.Zd),Dm.C(:))));
[np,nf] = chk(np,nf,'merged list contains every block''s variables',...
    all(cellfun(@(b) all(ismember(b.Zd,Dm.Zd)),D(:))));
info = verify(Dm);
[np,nf] = chk(np,nf,'verify(Dm) passes',info.true==1,strjoin(info.flags,' | '));
% A block swapped in with a foreign list must be caught.
Dbad = Dm;      Dbad.C{1,1} = D{1,1};
[np,nf] = chk(np,nf,'verify catches a block on a foreign list',...
    verify(Dbad).true==0);

%% ---------------------------------------------------------------- 2
fprintf('\n-- 2. operations preserve the invariant --\n');
Ds = Dm + mdopvar(D);
[np,nf] = chk(np,nf,'plus preserves it',...
    all(cellfun(@(b) isequal(b.Zd,Ds.Zd),Ds.C(:))) && verify(Ds).true==1);
Dt = Dm';
[np,nf] = chk(np,nf,'ctranspose preserves it',verify(Dt).true==1);
[np,nf] = chk(np,nf,'ctranspose transposes grid and metadata',...
    isequal(size(Dt),[2,2]) && isequal(Dt.space_out,Dm.space_in) ...
    && isequal(Dt.dim_out(:),Dm.dim_in(:)));
% Composed against a FIXED container, which is the case that arises:
% P*A, T'*P*T. G: (R,L2[s1]) -> (L2[s1],L2[s1,s2]).
G = mopvar(mkgrid(S,dm,din_D,out_D,pD,[2 3],9,true(2,2)));
Dc = Dm*G;
[np,nf] = chk(np,nf,'mtimes preserves it',...
    all(cellfun(@(b) isempty(b)||~isa(b,'sdopvar')||isequal(b.Zd,Dc.Zd),Dc.C(:))) ...
    && verify(Dc).true==1);
[np,nf] = chk(np,nf,'D*G takes D''s output and G''s input structure',...
    isequal(size(Dc),[2,2]) && isequal(Dc.space_out,Dm.space_out) ...
    && isequal(Dc.space_in,G.space_in) ...
    && isequal(Dc.dim_out(:),Dm.dim_out(:)) && isequal(Dc.dim_in(:),G.dim_in(:)));
[np,nf] = chk(np,nf,'scalar factor keeps the metadata',...
    verify(2*Dm).true==1 && isequal(size(2*Dm),[2,2]));

%% ---------------------------------------------------------------- 3
fprintf('\n-- 3. two decision factors are refused up front --\n');
% Decided on the argument types, before promotion and before the space
% checks, so it fires for both a space-compatible and an incompatible pair.
[np,nf] = chkerr(np,nf,'P*P'' refused',...
    'mdopvar:decisionTimesDecision',@() Dm*Dm');
[np,nf] = chkerr(np,nf,'refused before the space check too',...
    'mdopvar:decisionTimesDecision',@() Dm*Dm);

%% ---------------------------------------------------------------- 4
fprintf('\n-- 4. dispatch and promotion against mopvar --\n');
% Gsq is FIXED and square in spaces, so it can sit on either side of Dm.
Gsq = mopvar(mkgrid(S,dm,out_D,din_D,qD,pD,11,true(2,2)));
Sp1 = Dm + Gsq;
Sp2 = Gsq + Dm;
[np,nf] = chk(np,nf,'mdopvar + mopvar dispatches to mdopvar',...
    isa(Sp1,'mdopvar') && verify(Sp1).true==1);
[np,nf] = chk(np,nf,'mopvar + mdopvar dispatches to mdopvar',...
    isa(Sp2,'mdopvar') && verify(Sp2).true==1);
% Gsq' has input spaces out_D, which are Dm's output spaces, so this is the
% T'*P shape with the FIXED factor on the left.
Mp = Gsq'*Dm;
[np,nf] = chk(np,nf,'mopvar * mdopvar dispatches to mdopvar',...
    isa(Mp,'mdopvar') && verify(Mp).true==1);
[np,nf] = chk(np,nf,'explicit promotion keeps blocks and empties Zd',...
    isa(mdopvar(Gsq),'mdopvar') && isempty(mdopvar(Gsq).Zd) ...
    && verify(mdopvar(Gsq)).true==1);

%% ---------------------------------------------------------------- 5
fprintf('\n-- 5. mixed and all-fixed block classes --\n');
% Reached by ordinary arithmetic: adding a fixed operator to a decision
% operator that has a zero block leaves a 'sopvar' block among 'sdopvar'
% ones, because the zero block passes through addition untouched.
Dz = D;     Dz{1,2} = [];
Dmix = mdopvar(Dz) + Gsq;
[np,nf] = chk(np,nf,'mixed-class container is legal and verifies',...
    numel(unique(cellfun(@class,Dmix.C,'uni',0)))>1 && verify(Dmix).true==1);
[np,nf] = chk(np,nf,'a mixed container still composes',...
    verify(Dmix*G).true==1);
% An all-'sopvar' mdopvar is legal too, with an empty Zd. It is still
% refused as the second factor of a product with another mdopvar, by type.
Dfix = mdopvar(mkgrid(S,dm,out_D,din_D,qD,pD,13,true(2,2)));
[np,nf] = chk(np,nf,'all-fixed mdopvar is legal with empty Zd',...
    verify(Dfix).true==1 && isempty(Dfix.Zd));

%% ---------------------------------------------------------------- 6
fprintf('\n-- 6. the two containers reject each other''s blocks --\n');
[np,nf] = chkerr(np,nf,'mopvar rejects an sdopvar block',...
    'mopvar:badBlockClass',@() mopvar(D));
[np,nf] = chkerr(np,nf,'mdopvar rejects a non-operator block',...
    'mdopvar:badBlockClass',@() mdopvar({D{1,1},42}));
[np,nf] = chk(np,nf,'class states whether decision variables are present',...
    isa(Dm,'mdopvar') && ~isa(Gsq,'mdopvar') && isa(Gsq,'mopvar'));

%% ----------------------------------------------------------------
fprintf('\n===== %d passed, %d FAILED =====\n\n',np,nf);
if nf>0
    error('Test_mdopvar:failures','%d check(s) failed.',nf);
end
end

%% ======================= helpers =======================
function [np,nf] = chk(np,nf,name,tf,extra)
if nargin<5,    extra = '';     end
if tf
    np = np+1;      fprintf('   ok    %s\n',name);
else
    nf = nf+1;      fprintf('   FAIL  %s  %s\n',name,extra);
end
end

function [np,nf] = chkerr(np,nf,name,id,fh)
% Require the specific identifier, not merely that something threw: a
% different error would mean the guard under test never ran.
try
    fh();
    nf = nf+1;      fprintf('   FAIL  %s  (no error raised)\n',name);
catch ME
    if strcmp(ME.identifier,id)
        np = np+1;  fprintf('   ok    %s\n',name);
    else
        nf = nf+1;
        fprintf('   FAIL  %s  (got %s: %s)\n',name,ME.identifier,ME.message);
    end
end
end

function C = mkgrid(S,dm,out_idx,din_idx,q,p,seed,occ)
% M x N cell of random 'sopvar' blocks over the given output/input spaces.
rng(seed);
M = numel(out_idx);     N = numel(din_idx);
C = cell(M,N);
for i = 1:M
    for j = 1:N
        if ~occ(i,j),   continue,   end
        vs = struct('out',{S{out_idx(i)}},'in',{S{din_idx(j)}});
        ds = struct('out',repmat(dm,numel(vs.out),1),'in',repmat(dm,numel(vs.in),1));
        gs = struct('out',ones(1,numel(vs.out)),'in',ones(1,numel(vs.in)));
        C{i,j} = rand_sopvar([q(i),p(j)],vs,ds,gs,0.6);
    end
end
end
