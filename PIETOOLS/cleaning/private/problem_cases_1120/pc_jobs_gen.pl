# pc_jobs_gen.pl OUT - the job grid of the problem-case run of 11/20/2026
# (cleaning/private/problem_cases_1120/PLAN.md). Columns: tag, key, form,
# tier, mods, variant, param, ndvcap. mods carry no commas.
use strict; use warnings;
my @rows;
sub job { my ($tag,$key,$form,$tier,$mods,$cap) = @_; push @rows, [$tag,$key,$form,$tier,$mods,'A','',$cap//'']; }
my @marg = ('st.eppos = 1e-3','st.eppos = 1e-2','st.eppos = 1e-1',
            'st.epneg = 1e-4','st.epneg = 1e-3','st.epneg = 1e-2',
            'st.eppos = 1e-2; st.epneg = 1e-3');
my $SD = "st.sos_opts.solver = 'sdpt3'";
my $SE = "st.sos_opts.solver = 'sedumi'";
job('SANITY','PDE:41','P','stripped','st.eppos = 1e-3');
# --- PDE 41: heat with delayed boundary input (1-D; every start inconclusive)
for my $ft (['Q','light'],['P','light'],['P','stripped'],['Q','stripped']) {
    job('W41','PDE:41',$ft->[0],$ft->[1],$_) for @marg;
}
job('W41','PDE:41','Q','light',$SD);  job('W41','PDE:41','P','light',$SD);
job('W41','PDE:41','P','stripped',$SD);
job('W41','PDE:41','P','stripped',"st.eppos = 1e-2; $SE");
job('W41','PDE:41','Q','light',"st.eppos = 1e-2; $SE");
job('W41','PDE:41','Q','light',$_) for ('st.sop.deg_Q{1} = [6 2 2]','st.sop.deg_Q{1} = [8 4 4]',
    'st.sop.deg_Q{1} = [8 8 8]','st.sop.deg_Q{1} = [8 4 4]; st.lpi.ineq.dD = 0',
    'st.sop.deg_Q{1} = [2 2 2]; st.lpi.ineq.dD = 2');
job('W41','PDE:41','P','light','st.lpi.pos.len_sos{1} = [3;3;3]; st.lpi.ineq.dD = 0');
job('W41','PDE:41','P','veryheavy','st.lpi.ineq.dD = 0');
job('W41','PDE:41','P','light','st.lpi.pos.len_com{1} = 0; st.lpi.ineq.dD = 2');
job('W41','PDE:41','P','heavy','st.sop.dw_feas = [2 2 1]; st.lpi.ineq.dD = 0');
# --- PDE 44: wave with partly delayed boundary damping (1-D; infeasible)
job('W44','PDE:44','Q','light',$_) for ('st.sop.deg_Q{1} = [8 8 8]','st.sop.deg_Q{1} = [8 4 4]',
    'st.sop.deg_Q{1} = [10 6 6]','st.sop.deg_Q{1} = [8 8 8]; st.lpi.ineq.dD = 0',
    'st.sop.deg_Q{1} = [2 6 6]','st.sop.deg_Q{1} = [6 6 6]; st.lpi.ineq.dD = 2',
    'st.eppos = 1e-2','st.eppos = 1e-1','st.epneg = 1e-3',$SD);
job('W44','PDE:44','P','light',$_) for ('st.eppos = 1e-2','st.eppos = 1e-1','st.epneg = 1e-3',
    'st.lpi.ineq.dD = 2; st.lpi.pos.len_com{1} = 0');
job('W44','PDE:44','P','veryheavy',$_) for ('st.lpi.pos.len_sos{1} = [4;4;4]','st.lpi.ineq.dD = 0');
job('W44','PDE:44','P','heavy','st.lpi.pos.len_com{1} = 3; st.lpi.pos.len_inc{1} = 2');
# --- X3: the system of 44 with the delay on a second variable (2-D)
job('WX3','PDE:X3','Q','light',$_) for ('st.sop.deg_Q{1} = [8 8 8]','st.sop.deg_Q{2} = 2',
    'st.sop.deg_Q{1} = [8 4 4]; st.sop.deg_Q{2} = 2','st.eppos = 1e-2','st.epneg = 1e-3',$SD);
job('WX3','PDE:X3','P','light',$_) for ("st.settings_2d = settings_2d_stock('light')",'st.eppos = 1e-2');
job('WX3','PDE:X3','P','heavy','');
# --- PDE 36: 2-D telegraph (inconclusive); base B36 = the 52,880-variable start
my $B36 = 'st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1';
job('W36','PDE:36','Q','light',"$B36; $_") for (@marg[0..4],'st.sop.deg_Q{2} = [1 1 1]',
    'st.sop.deg_Q{2} = [2 2 2]','st.sop.deg_Q{2} = [0 1 1]',"st.lpi.ineq.psatz{2} = 'faces1'",
    "st.lpi.ineq.psatz{2} = 'none'","st.lpi.ineq.cap_inc = 'auto'",$SD,"st.eppos = 1e-2; $SE");
job('W36','PDE:36','Q','stripped','st.lpi.ineq.dD = 0');
job('W36','PDE:36','Q','heavy',"st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'");
# --- PDE 39: 2-D Navier-Stokes (feasibility ratio 2.64 at 316,737)
my $B39 = 'st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 1 1]; st.lpi.ineq.wcross = 1';
# 06:32: SeDuMi dropped: it raised to 1,828,161 variables (nnz(L) 2.3e8, about 20 min per iteration)
job('W39','PDE:39','Q','light',"$B39; $_") for (@marg[0..4],"st.lpi.ineq.psatz{2} = 'faces1'");
job('W39','PDE:39','Q','light','st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [1 1 1]; st.lpi.ineq.wcross = 1');
job('W39','PDE:39','Q','light',"st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 2 1]; st.lpi.ineq.wcross = 1; st.lpi.ineq.cap_inc = 'auto'");
# --- W36F (05:20): PDE 36 certified at epneg 1e-3 (W36); robustness and cheaper neighbours
job(q{W36F},q{PDE:36},q{Q},q{light},"$B36; $_") for (q{st.epneg = 2e-3},q{st.epneg = 5e-3},q{st.epneg = 1e-2},
    q{st.epneg = 3e-2},q{st.epneg = 1e-1},q{st.epneg = 1e-3; st.eppos = 1e-2},q{st.epneg = 1e-3; st.sop.deg_Q{2} = [1 1 1]},
    "st.epneg = 1e-3; st.lpi.ineq.psatz{2} = 'none'","st.epneg = 1e-3; st.lpi.ineq.cap_inc = 'auto'");
job(q{W36F},q{PDE:36},q{Q},q{stripped},q{st.lpi.ineq.dD = 0; st.epneg = 1e-3});
job(q{W36F},q{PDE:36},q{Q},q{light},q{st.epneg = 1e-3});
# --- wave 2, exclusive: stock 2-D storage in 'P' (certified PDE 65 on 11/19)
# 06:33: 36 moved to W2S36 (MATLAB access violation in its second pass, 751,936 variables; 36 certifies in Q)
job(q{W2S36},q{PDE:36},'P','light',"st.settings_2d = settings_2d_stock('light')");
job('W2S','PDE:69','P','light',"st.settings_2d = settings_2d_stock('light')");
job('W2S','PDE:70','P','light',"st.settings_2d = settings_2d_stock('light')");
job('W2S','PDE:39','P','light',"st.settings_2d = settings_2d_stock('heavy')");
# --- wave 2, exclusive: degree of Q against slack on 69 and 70
for my $k ('PDE:69','PDE:70') {
    job('W2Q',$k,'Q','light',$_) for ("st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2; st.lpi.ineq.cap_inc = 'auto'",
        "st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [2 1 1]; st.lpi.ineq.cap_inc = 'auto'",
        "st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = 'faces1'; st.lpi.ineq.wcross = 1",
        'st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.eppos = 1e-2');
    job('W2Q',$k,'P','stripped','st.lpi.ineq.dD = 0; st.eppos = 1e-2');
}
# --- wave 3, exclusive: PDE 56 probes within the workstation (gate 2e6, m 45,000)
job('W56','PDE:56','P','stripped','');
job('W56','PDE:56','Q','stripped',"st.lpi.ineq.cap_inc = 'auto'");
job('W56','PDE:56','Q','light','st.lpi.ineq.dD = 0; st.eppos = 1e-2');
job('W56','PDE:56','Q','light','st.lpi.ineq.dD = 0; st.epneg = 1e-3');
# --- WXC (05:30, cross terms where not tried with the new levers; small): X3, 39 in 'Q'
job('WXC','PDE:X3','Q','light',$_) for ('st.lpi.ineq.wcross = 1; st.epneg = 1e-3',
    'st.lpi.ineq.wcross = 1; st.lpi.ineq.dD = 0; st.epneg = 1e-3','st.lpi.ineq.wcross = 1; st.sop.deg_Q{1} = [8 8 8]',
    'st.lpi.ineq.wcross = 1; st.sop.deg_Q{2} = 2');
job('WXC','PDE:X3','P','light',$_) for ('st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1; st.epneg = 1e-3',
    'st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1; st.eppos = 1e-2');
job('WXC','PDE:39','Q','light',$_) for ('st.lpi.ineq.wcross = 1; st.lpi.ineq.dD = 0; st.epneg = 1e-3',
    'st.lpi.ineq.wcross = 1; st.lpi.ineq.dD = 0; st.epneg = 1e-2');
# --- W2X (05:30, exclusive): cross terms in the storage form and on 56, never solved before
job('W2X','PDE:39','P','light',$_) for ('st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1',
    'st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1; st.epneg = 1e-3');
job('W2X',$_,'P','light','st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1; st.eppos = 1e-2') for ('PDE:69','PDE:70');
job('W2X',$_,'Q','light',"st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2; st.lpi.ineq.cap_inc = 'auto'; st.lpi.ineq.wcross = 1") for ('PDE:69','PDE:70');
job('W2X','PDE:56','Q','light',$_) for ("st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.lpi.ineq.cap_inc = 'auto'",
    'st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1; st.epneg = 1e-3');
job('W2X','PDE:56','P','light','st.lpi.ineq.dD = 0; st.lpi.pos.len_cross{2} = 2; st.lpi.ineq.wcross = 1');
# --- W2P (05:35, exclusive): higher-degree Psatz terms on 69 and 70. Codes in 2-D: 0 plain, 1 box product,
#     3-6 faces (s_d - a_d), (b_d - s_d), 7-8 box quadratics (s_d - a_d)(b_d - s_d); 'corners' = faces plus the
#     products of two faces of different directions; wmin, dw raise the weight degree of every term.
my $CQ = "st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'";
for my $k ('PDE:69','PDE:70') {
    job('W2P',$k,'Q','light',"$CQ; $_") for ('st.lpi.ineq.psatz{2} = [0 3 4 5 6 7 8]',
        "st.lpi.ineq.psatz{2} = 'corners'",'st.lpi.ineq.psatz{2} = [0 1 3 4 5 6]',
        'st.lpi.ineq.wmin = [3 3 0]','st.lpi.ineq.wmin = [4 4 0]','st.lpi.ineq.dw = 2');
    job('W2P',$k,'P','stripped',$_) for ('st.lpi.ineq.dD = 0; st.lpi.ineq.psatz{2} = [0 3 4 5 6 7 8]',
        'st.lpi.ineq.dD = 0; st.lpi.ineq.wmin = [4 4 0]');
}
# --- W36V (05:55): the PDE 36 settings again with the executive status and eta recorded (diag_amt 11/20)
job('W36V','PDE:36','Q','light',"$B36; $_") for ('st.epneg = 1e-3','st.sop.deg_Q{2} = [0 1 1]','',
    'st.epneg = 1e-4','st.epneg = 2e-3','st.epneg = 5e-3','st.epneg = 1e-3; st.sop.deg_Q{2} = [0 1 1]');
# --- W39S (07:30, exclusive): PDE 39 certified in 'P' light with settings_2d_stock('heavy') at 1,692,801 (pass 1);
#     the slack holds the variables: slack cuts on that start, smallest expected first
my $SH = "st.settings_2d = settings_2d_stock('heavy')";
job('W39S','PDE:39','P','light',"$SH; $_") for ("st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'",'st.lpi.ineq.dD = 0');
job('W39S','PDE:39','P','stripped',$SH);
job('W39S','PDE:39','P','light',"$SH; $_") for ("st.lpi.ineq.cap_inc = 'auto'","st.lpi.ineq.psatz{2} = 'none'");
# --- W2SH (08:50; veryheavy rows 09:00; exclusive): the stock 2-D storage at heavy (the storage of the 39 certificate) on 69 and 70,
#     never tried (W2S ran them at the stock light level); 39 at the stock light level with dD 0: is heavy needed
job('W2SH',$_,'P','light',"$SH; st.lpi.ineq.dD = 0") for ('PDE:69','PDE:70');
job('W2SH','PDE:39','P','light',"st.settings_2d = settings_2d_stock('light'); st.lpi.ineq.dD = 0");
job('W2SH',$_,'P','light',$SH) for ('PDE:69','PDE:70');
job(q{W2SH},$_,q{P},q{light},"st.settings_2d = settings_2d_stock('veryheavy'); st.lpi.ineq.dD = 0") for (q{PDE:69},q{PDE:70});
# --- W2E (11:20, exclusive): deg_Q{2} = [0 1 1], the storage degree that certifies PDE 38 (2-D undamped wave, position
#     state u, 'Q' light dD 0, 380,804) and PDE 36; never run on 69 or 70; cap first (smaller), then uncapped, then [0 2 2]
my $Q0 = "st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 1 1]";
job('W2E',$_,'Q','light',"$Q0; st.lpi.ineq.cap_inc = 'auto'") for ('PDE:70','PDE:69');
job('W2E',$_,'Q','light',$Q0) for ('PDE:70','PDE:69');
job('W2E',$_,'Q','light',"st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = [0 2 2]") for ('PDE:70','PDE:69');
# --- W2D (11:20, exclusive): PDE 70 at deg_Q{2} = 2 (2,015,785 variables, m 35,390; W2Q stopped it at the 2e6 gate)
job('W2D','PDE:70','Q','light',"st.lpi.ineq.dD = 0; st.sop.deg_Q{2} = 2; st.lpi.ineq.cap_inc = 'auto'",2.1e6);
# --- W2R (11:23, exclusive): a second slack-weight raise in form 'Q' on 70 (st.sop.max_raise_Q; st.lpi.ineq.dw is not
#     read by the stability executive); 70 pass 2 has 457,393 variables, so pass 3 should stay under the gate (69 not)
job('W2R','PDE:70','Q','light',"st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'; st.sop.max_raise_Q = [2 2 0]");
# --- W1V, W2V (11:26): the dual forms, never run on a problem case. PIETOOLS_stability_sop: 'Pd' (P >= eppos I,
#     T*P*A' + A*P*T' + 2*epneg*T*P*T' <= 0) in every dimension, 'Qd' in 1-D only (posed as 'Pd' in two or more
#     variables). PDE 18 (1-D damped wave in the position state u) certifies only in 'Pd' (336 variables, 11/19).
#     W1V: the 1-D cases 44, 41 and X3 (two variables); W2V: 70, 69, 56, smallest programs first.
for my $k ('PDE:44','PDE:41') {
    job('W1V',$k,$_->[0],$_->[1],'') for (['Pd','light'],['Qd','light'],['Pd','stripped'],['Qd','stripped']);
}
job('W1V','PDE:X3','Pd',$_,'') for ('light','stripped');
# 11:31: split by size after DIAG69 (the 1-D position forms of 69 and 70 certify in 'Pd' only): W2VS the five
#     programs with dD 0 on a second runner now, W2VL the two at dD 1 after W2E on runner 8 (memory)
job('W2VS',$_,'Pd','stripped','st.lpi.ineq.dD = 0') for ('PDE:70','PDE:69');
job('W2VS',$_,'Pd','light','st.lpi.ineq.dD = 0') for ('PDE:70','PDE:69','PDE:56');
job('W2VL',$_,'Pd','light','') for ('PDE:70','PDE:69');
# --- W56D (11:34): PDE 56 with a raised storage degree. Its 1-D twin, PDE 43 (delay as a 2-D transport state),
#     certifies only at deg_Q{2} = [1 2 2] ('Q' light dD 0 with the cap, 132,609; 11/19 EXS2); 56 is a 3-variable PIE
#     without 1-D blocks, so every block of its Q takes deg_Q{3} (PIETOOLS_stability_sop, form 'Q'); never varied
job('W56D','PDE:56','Q','light',"st.lpi.ineq.dD = 0; st.lpi.ineq.cap_inc = 'auto'; $_") for ('st.sop.deg_Q{3} = [1 2 2]','st.sop.deg_Q{3} = 2');
# --- W2V38 (11:33): calibration – PDE 38 (2-D undamped wave in the position state, certified in 'Q' light dD 0
#     deg_Q{2} [0 1 1]) in the dual form; 70 'Pd' stripped dD 0 infeasible (362,596)
job('W2V38','PDE:38','Pd',$_,'st.lpi.ineq.dD = 0') for ('stripped','light');
# --- W55 (12:41): PDE 55 (3-D coupled heat, |lam| = 1 against 16.2) – P = I is an energy certificate (C2*3*pi^2 = 8.9
#     > lam), yet 'P' stripped dD 0 (114,704) is infeasible and dD 1 (7.34M) ran out the Sol limit. The slack weight in
#     three variables was never raised: the 11/19 row with st.lpi.ineq.wmin [0 0 1] changed nothing (not read by the
#     executive); the executive takes dwbase = st.sop.dw_feas(min(nv,end)) (stripped [0 1 1]: 1 in 3-D) and no raise
#     (max_raise_feas [1 1 0]). Rows: the 3-D entry of dw_feas at 2 and 3.
job('W55','PDE:55','P','stripped',"st.lpi.ineq.dD = 0; $_") for ('st.sop.dw_feas = [0 1 2]','st.sop.dw_feas = [0 1 3]');
# --- W55C (12:55): PDE 55 at dD 1 between the 'auto' cap (1,808, infeasible) and no cap (7.34M): cap_inc takes one
#     cap per direction ('auto' = Dmax - 1; GET_LIFT_DEGS); W55: dD 0 with dw_feas 2, 3 infeasible (284,944; 684,304)
job('W55C','PDE:55','P','stripped',"st.lpi.ineq.cap_inc = $_") for ('[2 2 2]','[3 3 3]','[4 4 4]');
# --- W50S (13:13): the PDE 50 H2 dual 'cco' (W >= eppos I, (A*W)*T' + T*(W*A') + B1*B1' <= 0) is infeasible at every
#     start; a coercive W there is a 'Pd' stability certificate (P = W), so 'cco' needs 'Pd' feasible. The stability
#     test of the PDE 50 system in the dual and primal forms:
job('W50S','PDE:50',$_->[0],$_->[1],$_->[2]) for (['Pd','light',''],['Pd','stripped','st.lpi.ineq.dD = 0'],['Qd','light',''],['P','light','st.lpi.ineq.dD = 0']);
# --- W43P (13:43): PDE 43 (1-D twin of 56) has been tested in 'Q' only (certified at deg_Q{2} [1 2 2], 132,609);
#     in 'P' a cheaper route for 56 would show here first
job('W43P','PDE:43','P',$_->[0],$_->[1]) for (['light','st.lpi.ineq.dD = 0'],['light',''],['heavy','st.lpi.ineq.dD = 0'],['stripped','']);
open(my $fh,'>',$ARGV[0]) or die;
print $fh "tag,key,form,tier,mods,variant,param,ndvcap\n";
for my $r (@rows) { die "comma in mods: $r->[4]\n" if $r->[4] =~ /,/;
    print $fh join(',',@$r[0..3],"\"$r->[4]\"",@$r[5..7]),"\n"; }
close $fh;
my %n; $n{$_->[0]}++ for @rows; print join(' ', map {"$_=$n{$_}"} sort keys %n), "\n";
