#!/usr/bin/perl
# make_starts.pl - the starts files of the settings rerun (PLAN.md §3):
#   starts_P1.csv   two and three variables (1a, 1b);
#   starts_P1d.csv  their dual tests ('Qd' for a 'Q' example, 'Pd' for 'P');
#   starts_P2.csv   one variable (2a, 2b).
# Columns key,kind,form,tier,mods,best (ex_search). kind and best from
# cleaning/private/example_settings_1119/phaseA_choices.csv; a PDE row starts
# from the st its example file returns (mods 'examplest'); DDE, NDS and DDF
# rows from the tier and changes of PHASE_C_TABLE.md. MMP, 11/21/2026.
use strict; use warnings;
my $root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS/cleaning/private';
my $out  = "$root/settings_rerun_1121";
my (%kind, %best, %formA);
open(my $A, '<', "$root/example_settings_1119/phaseA_choices.csv") or die "phaseA: $!";
my $h = <$A>;
while (my $l = <$A>) {
    $l =~ s/\r?\n$//;
    my @v = ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g);
    my $key = "$v[0]:$v[1]";
    $kind{$key} = $v[2]; $formA{$key} = $v[3]; my $b = $v[7] // ''; $b =~ s/"//g; $best{$key} = $b;
}
close $A;
# Phase C table: case -> (tier, changes, form)
my (%tierC, %modsC, %formC);
open(my $C, '<', "$root/example_settings_1119/PHASE_C_TABLE.md") or die "table: $!";
while (my $l = <$C>) {
    next unless $l =~ /^\| /;
    my @c = map { s/^\s+|\s+$//gr } split /\|/, $l;
    next if @c < 5 || $c[2] eq 'Case' || $c[2] =~ /^-+$/;
    my ($case, $form, $tc) = ($c[2], $c[3], $c[4]);
    next unless $tc =~ /^(stripped|light|heavy|veryheavy|extreme)\s*:\s*(.*)$/ || $tc =~ /^(stripped|light|heavy|veryheavy|extreme)\s*$/;
    my ($tier, $ch) = ($1, $2 // '');
    $ch =~ s/^\s*–\s*$//;
    my @parts = grep { /\S/ } split /;\s*(?![^\[\{\(]*[\]\}\)])/, $ch;   # a ';' inside brackets separates entries
    my $mods = join('; ', map { "st.$_" } @parts);
    $tierC{$case} = $tier; $modsC{$case} = $mods; $formC{$case} = $form;
}
close $C;
# The form each test requests in the current example files (survey of 11/21)
my %formS;
open(my $S, "<", "$root/stability_form_swap_1121/results_SWAPSURVEY.csv") or die "survey: $!";
<$S>;
while (my $l = <$S>) { my @v = split /,/, $l; $formS{$v[0]} = $v[4]; }
close $S;
sub start_row {
    my ($key, $form, $tier, $mods) = @_;
    my $k = $kind{$key} || 'stab'; my $b = $best{$key}; $b = 1 if !defined $b || $b eq '';
    $mods =~ s/"/""/g;
    return sprintf("%s,%s,%s,%s,\"%s\",%s\n", $key, $k, $form, $tier, $mods, $b);
}
sub form_of {
    my ($key) = @_;
    # PDE rows: the form the example file names (survey); the harness copies of
    # the DDE, NDS and DDF blocks predate the 11/19 forms, so those rows take
    # the Phase A form
    return $formS{$key} if $key =~ /^PDE:/ && defined $formS{$key} && $formS{$key} ne "";
    my $f = $formA{$key} // '';
    return $f if $f ne '';
    my ($case) = $key =~ /^[A-Z]+:([^\/]+)/;
    my $fc = $formC{$case} // ''; $fc =~ s/^(stab|gain|value)\s+//;
    return $fc ne '' ? $fc : 'Q';
}
sub case_start {             # the example st: PDE rows examplest, else the table
    my ($key) = @_;
    my ($lib, $case) = $key =~ /^([A-Z]+):([^\/]+)/;
    return ('light', 'examplest') if $lib eq 'PDE';
    return ($tierC{$case} // 'light', $modsC{$case} // '');
}
# Priority 1 (two and three variables)
my %alts = (
  'PDE:36' => [['Q','light','st.lpi.ineq.dD = 0; st.lpi.ineq.wcross = 1']],
  'PDE:39' => [['Q','light',''], ['P','light',"st.settings_2d = settings_2d_stock('light'); st.lpi.ineq.dD = 0"]],
  'PDE:53' => [['P','light',''], ['P','heavy','']],
  'PDE:54' => [['P','light',''], ['P','heavy','']],
  'PDE:55' => [['P','light',''], ['P','heavy','']],
  'PDE:56' => [['Q','light','st.lpi.ineq.dD = 0']],
);
my @P1 = qw(PDE:36 PDE:39 PDE:50 PDE:50/H2_norm_dual PDE:53 PDE:54 PDE:55 PDE:56 PDE:57 PDE:62 PDE:63 PDE:65 PDE:72);
open(my $O, '>', "$out/starts_P1.csv") or die; print $O "key,kind,form,tier,mods,best\n";
for my $key (@P1) {
    for my $a (@{$alts{$key} // []}) { print $O start_row($key, @$a); }
    my ($t, $m) = case_start($key);
    print $O start_row($key, form_of($key), $t, $m);
}
close $O;
# Priority 1d: the dual test of each certified two- and three-variable stability row
my @P1d = qw(33 34 35 37 38 40 42 43 57 58 59 60 61 62 63 64 66);
open($O, '>', "$out/starts_P1d.csv") or die; print $O "key,kind,form,tier,mods,best\n";
for my $id (@P1d) {
    my $key = "PDE:$id"; my $f = form_of($key);
    my $d = ($f eq 'P') ? 'Q>Pd' : 'Q>Qd';
    print $O start_row($key, $d, 'light', 'examplest');
    print $O start_row($key, $d, 'light', '');
    print $O start_row($key, $d, 'heavy', '');
}
close $O;
# Priority 2 (one variable)
my @P2 = qw(DDE:DDE18 PDE:31 DDE:DDE20 PDE:32 PDE:48 DDE:DDE14 DDE:DDE24 DDE:DDE11/Hinf_gain_dual PDE:49
            NDSDDF:NDS11 DDE:DDE10/Hinf_gain_dual PDE:9/stability_dual DDE:DDE26 DDE:DDE25/H2_norm_dual DDE:DDE22
            PDE:8/stability_dual PDE:52 PDE:30 NDSDDF:NDS11/stability_dual DDE:DDE25 DDE:DDE13 DDE:DDE10 DDE:DDE09
            PDE:9 PDE:49/H2_norm_dual PDE:28 PDE:22 DDE:DDE23 PDE:15 NDSDDF:DDF03 DDE:DDE15 DDE:DDE11 PDE:41
            PDE:27 PDE:21/Hinf_gain_dual PDE:13.2 PDE:13.2/stability_dual);
open($O, '>', "$out/starts_P2.csv") or die; print $O "key,kind,form,tier,mods,best\n";
for my $key (@P2) {
    my ($t, $m) = case_start($key);
    print $O start_row($key, form_of($key), $t, $m);
}
close $O;
print "P1 ", scalar(@P1), " keys; P1d ", scalar(@P1d), "; P2 ", scalar(@P2), "\n";
