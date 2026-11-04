function lb_diag()
% LB_DIAG, 10/28/2026: the reruns behind the sopvar failure diagnosis.
%   1. DDE19, DDF03: the output and input dimensions of the PIE (nz = 0 is
%      the hypothesis), and DDE19 with its legacy output fields renamed to
%      the current DDE names (C0 -> C1, C{1} -> C1i{1}, D1 -> D11, D2 -> D12)
%      on both paths at light (documented 1.964);
%   2. PDE 42, 43, 50: the state spaces of the 2-D PIE (T.dim rows R,
%      L2[s1], L2[s2], L2[s1,s2]);
%   3. PDE 45 (exact 0.11101): the sopvar H-infinity LPI with the stock
%      equality slack ('stock', settings_PIETOOLS_stripped_2D, which declares
%      the eq fields) against stock stripped (0.25218, unknown), and the
%      default slack ('like') at light with dw_P raised from 2 to 4. That
%      second test is void: form 'Q' reads dw_Q, not dw_P
%      (PIETOOLS_Hinf_gain_sop, dwbase), and the program was unchanged.
here = fileparts(mfilename('fullpath'));
lb_setup();
S = load(fullfile(here,'catalog.mat'));  C = S.C;
cat_of = @(k) C(find(strcmp({C.key},k),1));
pr = @(lbl,r) fprintf('%-34s %-10s gam %-12.8g %7.1f s  ndec %-8d m %-6d fr %.3g numerr %g\n', ...
    lbl, r.status, r.gam, r.sec, r.ndec, r.m, r.fr, r.numerr);

% 1. Output and input dimensions
for k = {'DDE:DDE19','NDSDDF:DDF03'}
    try
        PIE = lb_build_system(cat_of(k{1}), [], true);
        fprintf('%s: T %s | B1 %s | B2 %s | C1 %s | D11 %s | D12 %s\n', k{1}, ...
            mat2str(PIE.T.dim), mat2str(PIE.B1.dim), mat2str(PIE.B2.dim), ...
            mat2str(PIE.C1.dim), mat2str(PIE.D11.dim), mat2str(PIE.D12.dim));
    catch ME
        fprintf('%s: %s\n', k{1}, ME.message);
    end
end
% DDE19 with the current field names
try
    code = fileread(fullfile(here,'blocks','DDE19.m'));
    code = strrep(code,'DDE.C0=','DDE.C1=');     code = strrep(code,'DDE.C{1}=','DDE.C1i{1}=');
    code = strrep(code,'DDE.D1=','DDE.D11=');    code = strrep(code,'DDE.D2=','DDE.D12=');
    sys = lb_eval_block(code,'');
    sys = initialize_PIETOOLS_DDE(sys);  PIE = convert_PIETOOLS_DDE(sys,'pie');
    fprintf('DDE19 renamed: C1 %s | D11 %s | D12 %s\n', mat2str(PIE.C1.dim), mat2str(PIE.D11.dim), mat2str(PIE.D12.dim));
    pr('DDE19 renamed, stock light', lb_exec('PIETOOLS_Hinf_control', PIE, lpisettings('light',0,'','mosek')));
    pr('DDE19 renamed, sop light', lb_exec('sop:PIETOOLS_Hinf_control', PIE, lpisettings_sop('light',0,'','mosek')));
catch ME
    fprintf('DDE19 renamed: %s (%s)\n', ME.message, ME.stack(1).name);
end

% 2. State spaces of the 2-D PIEs
for k = {'PDE:42','PDE:43','PDE:50'}
    try
        PIE = lb_build_system(cat_of(k{1}), [], true);
        fprintf('%s: T.dim rows R, L2[s1], L2[s2], L2[s1,s2] = %s\n', k{1}, mat2str(PIE.T.dim));
    catch ME
        fprintf('%s: %s\n', k{1}, ME.message);
    end
end

% 3. PDE 45: stock slack at the stripped 2-D degrees, then 'like' at dw_P 4
PIE = lb_build_system(cat_of('PDE:45'), [], true);
try
    St = lpisettings_sop('light',0,'','mosek');
    St.settings_2d = settings_PIETOOLS_stripped_2D();
    St.sop.slack = 'stock';     St.sop.max_raise = 0;
    pr('PDE45 sop, stock slack, stripped 2D', lb_exec('sop:PIETOOLS_Hinf_gain', PIE, St));
catch ME
    fprintf('PDE45 stock slack: %s (%s)\n', ME.message, ME.stack(1).name);
end
try
    St = lpisettings_sop('light',0,'','mosek');
    St.sop.dw_P = 4;            St.sop.max_raise = 0;
    pr('PDE45 sop, like slack, dw_P 4', lb_exec('sop:PIETOOLS_Hinf_gain', PIE, St));
catch ME
    fprintf('PDE45 dw_P 4: %s (%s)\n', ME.message, ME.stack(1).name);
end
end
