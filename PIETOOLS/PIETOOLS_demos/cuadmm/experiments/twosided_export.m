% twosided_export.m -- write a dump with a DELIBERATELY WRONG svec factor.
%
% The one-sided control in cuimport_test proved nothing about the trap that
% matters.  Corrupting only the import is not an inverse pair, and rel_b caught
% it at once (6.9e-01 and 9.7e-01).  The real hazard is the CLAUDE.md item 4
% pattern: the SAME error in export and import, which cancels in the row
% residual.  That requires re-exporting, re-solving, and re-importing with the
% matched wrong factor.
%
% Factor 0.5 is chosen, not 1.0, because it is the ADVERSARIAL direction: the
% one-sided run showed 0.5 shrinks off-diagonals toward a diagonal matrix and
% left psd_min POSITIVE (+5.4e-09), while 1.0 broke positivity outright
% (-1.2e-01).  If the PSD check has power it must catch the case that already
% survived it once.
%
% PREDICTION, recorded before running: rel_b will be SMALL (the linear system is
% consistently transformed on both sides and the error cancels) but psd_min will
% be BAD, because cuADMM projects onto the PSD cone using ITS OWN fixed sqrt(2)
% svec convention internally -- not the factor I exported with -- so the cone it
% enforces is not the cone I intended.  If psd_min comes back clean too, the
% verification has no power against this trap and needs a different check.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
MF = fullfile(cuadmm_outdir(),'cu_fish','fish_nobnd.mat');
S = load(MF);
K = S.K; Ks = double(K.s(:)'); Kf = K.f;
nvar = size(S.At,1);  vec_len = Kf + sum(Ks.*(Ks+1)/2);

FAC = 0.5;
Tw = build_T(Kf,Ks,vec_len,nvar,FAC);
out = fullfile(cuadmm_outdir(),'cu_fish','fish_w05');
if ~exist(out,'dir'), mkdir(out); end
store_coo(Tw*S.At, fullfile(out,'At.txt'));
store_coo(full(S.b(:)), fullfile(out,'b.txt'));
store_coo(Tw*S.c(:),  fullfile(out,'C.txt'));
fid = fopen(fullfile(out,'blk.txt'),'w');
if Kf>0, fprintf(fid,'u %d\n',Kf); end
for i=1:numel(Ks), fprintf(fid,'s %d\n',Ks(i)); end
fclose(fid);
fid = fopen(fullfile(out,'con_num.txt'),'w'); fprintf(fid,'%d\n',numel(S.b)); fclose(fid);
fprintf('TW wrote|%s|factor=%.4f|vec_len=%d|nnz=%d\n',out,FAC,vec_len,nnz(Tw*S.At));
fprintf('TWDONE\n');

function T = build_T(Kf,Ks,vec_len,nvar,r2)
I={};J={};V={};
if Kf>0, p=(1:Kf)'; I{end+1}=p; J{end+1}=p; V{end+1}=ones(Kf,1); end
bs=Kf; bv=Kf;
for k=1:numel(Ks)
    N=Ks(k); idx=find(triu(true(N))); [ii,jj]=ind2sub([N N],idx);
    p=(1:numel(idx))'; dg=(ii==jj); od=~dg;
    I{end+1}=[bs+p; bs+p(od)];                                       %#ok<AGROW>
    J{end+1}=[bv+idx; bv+sub2ind([N N],jj(od),ii(od))];              %#ok<AGROW>
    V{end+1}=[dg + od*r2; repmat(r2,nnz(od),1)];                     %#ok<AGROW>
    bs=bs+N*(N+1)/2; bv=bv+N^2;
end
T=sparse(cat(1,I{:}),cat(1,J{:}),cat(1,V{:}),vec_len,nvar);
end

function store_coo(A,fn)
% COO, 0-based, "row col val", rows ascending -- the format dump2cuadmm writes.
A = sparse(A);
[i,j,v] = find(A);
[~,p] = sortrows([i j]);  i=i(p); j=j(p); v=v(p);
fid = fopen(fn,'w');
fprintf(fid,'%d %d %.17g\n',[i-1 j-1 v]');
fclose(fid);
end
