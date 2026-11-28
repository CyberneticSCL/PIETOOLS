# rename_table.tsv from cands.tsv (name, path, new, 2025 file) and counts.txt (files/lines by area).
use strict; use warnings;
my ($cands, $counts) = @ARGV;
my %keep = ('exec_sop_settings' => 'keep: names the settings field st.sop',
            'lb_full_sop' => 'keep: dated baseline helper; lb_full exists beside it',
            'lb_sop_plan' => 'keep: dated baseline helper (records)',
            'bl_sop' => 'keep: test helper, the sop variant of a suite');
my %new = ('eq_opts_sopvar' => 'get_eq_opts');
my %note = ('eq_opts_sopvar' => '2025 had get_eq_opts_2D (executives/2D/utility); serves dopvar and sdopvar',
            'PIETOOLS_stability_sop' => '2025: four executives PDEstability(_dual), PIE2PDEstability(_dual) = forms P, Pd, Q, Qd',
            'PIETOOLS_H2_norm_sop' => '2025: H2_norm_o, _c, _o_coercive, _c_coercive = forms o, c, oco, cco',
            'PIETOOLS_Hinf_gain_sop' => '2025 same name = form Q; also _coercive, _dual, _dual_coercive = P, Qd, Pd',
            'PIETOOLS_Hinf_OFcontrol_sop' => 'new in this release',
            'getController_sop' => '2025 third input tol (numeric); now opts (struct)',
            'getObserver_sop' => '2025 third input tol (numeric); now opts (struct)',
            'lpi_ineq_sop' => '2025 opts.psatz 0/1, opts.pure, opts.sep; now term sets, no pure/sep',
            'poslpivar_settings_sop' => 'C7 of the container rename: named in this pass');
my %cnt; open my $c, '<', $counts or die; <$c>;
while (<$c>) { chomp; s/\r//; my @f = split ' '; next unless @f >= 6 && $f[0] !~ /^LINES/; $cnt{$f[0]} = [@f[1..5]] } close $c;
open my $in, '<', $cands or die; <$in>;
print join("\t", qw(current kind location new_name 2025_counterpart toolbox_f/l doc_f/l tests_f/l note)), "\n";
my %done;
while (<$in>) { chomp; s/\r//; my ($n, $p, $nw, $y) = split /\t/;
  next if $done{"$n\t$p"}++;
  my $kind = $p =~ m{^testing_codes/} ? 'test' : 'function';
  my $new = $keep{$n} ? $n : ($new{$n} // $nw);
  my $k = $cnt{$n} // ['0/0','0/0','0/0'];
  print join("\t", $n, $kind, $p, $new, $y, $k->[0], $k->[1], $k->[2], $keep{$n} // $note{$n} // ''), "\n" }
close $in;
for my $d (['executives_sopvar','executives'], ['lpi_programming_sopvar','lpi_programming'], ['settings_sopvar','settings']) {
  my $k = $cnt{$d->[0]} // ['0/0','0/0','0/0'];
  print join("\t", $d->[0], 'folder', "$d->[0]/ (and testing_codes/private/ where present)", $d->[1], "$d->[1]/", @$k[0..2], 'retired 2025 folder of the same name in cleaning/private/ (off the path)'), "\n" }
