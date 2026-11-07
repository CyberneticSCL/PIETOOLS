% regress_lpieq.m -- does the lpi_eq nnz change alter any existing result?
%
% lpi_eq is on the path of EVERY case in the suite, and the baseline table was
% measured BEFORE the patch (only the n=24/32 rungs were built after it). The
% change is meant to be exactly equivalent -- nnz(C.C)~=0 and ~all(all(C.C==0))
% both ask whether every entry is zero -- but "meant to be" is not a measurement.
% Re-run a spread of cases and compare m, K.f, blocks and rel_b against the
% recorded values: Kf=43 and Kf=0 variants, an objective case that carries an
% 'ineq' expression, the only case declaring two positive operators, and a 2-D
% case. A single differing digit means the guard was not equivalent in practice.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
REF = containers.Map();
% The reference MUST be the frozen pre-patch table shipped in results/, not the
% live one under cuadmm_outdir(): that one holds whatever bl_run last wrote, so
% after a run on current code this check would compare the code against itself.
REFTSV = fullfile(fileparts(HERE),'results','2026-09-24','bl_mosek.tsv');
fid = fopen(REFTSV,'r');
if fid < 0, error('regress_lpieq:ref','reference table not found: %s',REFTSV); end
hdr = fgetl(fid);
H = strsplit(hdr,sprintf('\t'));
ix = @(n) find(strcmp(H,n));
while true
    l = fgetl(fid); if ~ischar(l), break; end
    p = strsplit(l,sprintf('\t'),'CollapseDelimiters',false);
    REF(p{1}) = p;
end
fclose(fid);

CASES = {'stab_rd1','stabpde_rd1','hinf_rd1','wellposed_rd1','stab2_rd'};
C = bl_cases();
ids = {C.id};
fprintf('RG id|m_ref|m_new|Ks_ref|Ks_new|relb_ref|relb_new|rel_reldiff|MATCH\n');
allok = true;
for k = 1:numel(CASES)
    id = CASES{k};
    c = C(strcmp(ids,id));
    evalc('[sol,M] = cuadmm_private(c.builder,c.args{:});');
    S = cuadmm_private('sdpshape',sol);
    Atf=[];bf=[];
    for q=1:sol.expr.num, Atf=[Atf,sol.expr.At{q}]; bf=[bf;sol.expr.b{q}]; end
    x = sol.solinfo.RRx(:);
    rb = norm(full(Atf'*x-bf))/max(norm(full(bf)),eps);
    p = REF(id);
    m_ref = str2double(p{ix('m')});  Ks_ref = p{ix('Ks')};
    rb_ref = str2double(p{ix('rel_b')});
    d = abs(rb-rb_ref)/max(rb_ref,eps);
    ok = (S.m==m_ref) && strcmp(mat2str(S.Ks),Ks_ref) && d < 1e-6;
    allok = allok && ok;
    fprintf('RG %s|%d|%d|%s|%s|%.6e|%.6e|%.2e|%d\n', ...
        id,m_ref,S.m,Ks_ref,mat2str(S.Ks),rb_ref,rb,d,ok);
end
fprintf('RGVERDICT %s\n', string(allok));
fprintf('RGDONE\n');
