% examples_DDE_library_PIETOOLS.m wind tunnel block, one comment level removed (11/19/2026: C1, C1i, D11, D12; Ai)
 % % % % Wind Tunnel Problem from [4], adapted from [5]
 Hinf_control=1
 aM=1/1.964; kM=-.0117;zM=.8;wM=6;
 DDE.A0=[-aM 0 0; 0 0 1; 0 -wM^2 -2*zM*wM];
 %  DDE.Adi{1}=[0 kM*aM 0; 0 0 0; 0 0 0];                                  % MMP, 11/19/2026 (was)
  DDE.Ai{1}=[0 kM*aM 0; 0 0 0; 0 0 0];                                     % MMP, 11/19/2026
 % [5] has the discrete delay x2(t-tau): kappa m' + m = k theta(t-tau0),   % MMP, 11/19/2026
 % eq. (6.7) of Mazanti, Boussaada, Niculescu, arXiv 2002.06128: Ai{1},    % MMP, 11/19/2026
 % the discrete delay, in place of Adi{1}, a distributed delay (since      % MMP, 11/19/2026
 % 63831988, 12/2022; Ad{1} before, outside the struct).                   % MMP, 11/19/2026
 % z = [x1; x2; .1u]: the output of this block since PIETOOLS 2021b, in    % MMP, 11/19/2026
 % field names initialize_PIETOOLS_DDE no longer reads (C0, C{1}, D1,      % MMP, 11/19/2026
 % D2: z was 0). Neither [5] nor arXiv 1806.08071 ([4]) defines z for      % MMP, 11/19/2026
 % this model; the choice follows the examples of [4], z = [x; .1u].       % MMP, 11/19/2026
 %  DDE.C0=[1 0 0; 0 1 0;0 0 0]; DDE.C{1}=[0 0 0;0 0 0;0 0 0];             % MMP, 11/19/2026 (was)
  DDE.C1=[1 0 0; 0 1 0;0 0 0]; DDE.C1i{1}=[0 0 0;0 0 0;0 0 0];             % MMP, 11/19/2026
  DDE.B2=[0;0;wM^2]; %B in paper
  DDE.B1=[1 0;0 0;0 10]; % D in paper
  DDE.tau(1)=.33;
 %  DDE.D1=[0 0;0 0;0 0];                                                  % MMP, 11/19/2026 (was)
  DDE.D11=[0 0;0 0;0 0];                                                   % MMP, 11/19/2026
 %  DDE.D2=[0;0;.1];                                                       % MMP, 11/19/2026 (was)
  DDE.D12=[0;0;.1];                                                        % MMP, 11/19/2026
 % st: example-specific settings (11/19/2026); closed-loop gain            % MMP, 11/19/2026
 %     1.963986, 265 decision variables (stripped as it is: 1,252)         % MMP, 11/19/2026
 st = lpisettings_sop('stripped');                                         % MMP, 11/19/2026
 st.lpi.pos.len_inc{1} = 0;                                                % MMP, 11/19/2026
 st.lpi.pos.len_com{1} = 0;                                                % MMP, 11/19/2026
 st.lpi.pos.len_sos{1} = [1;2;2];                                          % MMP, 11/19/2026
 st.lpi.ineq.dD = 0;                                                       % MMP, 11/19/2026
 st.lpi.var.deg{1} = [0 2 2];                                              % MMP, 11/19/2026
 % %  gam=1.9369   @6,12, eps1=eps2=.000001, tau=.33 IPM=6.61
 % %  gam=1.964   @4,8, eps1=eps2=.000001, tau=.33 IPM=1.98
 % %  gam=1.964   @2,4, eps1=eps2=.000001, tau=.33 IPM=.743
 % %  % C=[0 10]
