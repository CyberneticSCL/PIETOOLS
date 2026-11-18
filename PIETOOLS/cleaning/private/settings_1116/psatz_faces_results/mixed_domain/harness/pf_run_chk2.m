% PF_RUN_CHK2  The family check with both families (Step 2 oscillatory,
% Step 1 boundary layer): S17 under SeDuMi (N-D default storage, the M2-type
% 'inaccurate' certificate); D4 under MOSEK with the stock storage and the
% corner rows at offset 0 and at offset 1 (the 11/12 record certified D4 at
% 664k variables with the corner rows then on by default); S3 under MOSEK
% with the N-D default storage (control). SCRATCH, MMP, 11/16/2026.
pf_check_fam('def','sedumi',3);
C0 = pf_cfg({'corners0','corners'},2);
pf_check_fam('stock','mosek',1,C0(1));
pf_check_fam('stock','mosek',1,C0(2));
pf_check_fam('def','mosek',2);
disp('## CHK2 ALL DONE');
