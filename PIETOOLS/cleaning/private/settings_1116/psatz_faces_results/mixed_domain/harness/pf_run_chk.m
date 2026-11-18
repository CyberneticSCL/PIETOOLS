% PF_RUN_CHK  The family check of the certificates: S17 (SeDuMi, the M2
% 'inaccurate' certificate), D4 and S3 under SeDuMi and under MOSEK, N-D
% default and stock storages. SCRATCH, MMP, 11/16/2026.
pf_check_fam('def','sedumi',[3 1 2]);
pf_check_fam('stock','mosek',[1 2]);
pf_check_fam('def','mosek',[1 2]);
disp('## CHK ALL DONE');
