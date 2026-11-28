# Candidates (basename contains _sop as a suffix token, or settings_sopvar_*, eq_opts_sopvar, exec_sop_settings):
# proposed new name, 2025 file of that name, current files of that name (collision).
use strict; use warnings;
my ($allf, $f25) = @ARGV;
my (%cur, %y25);
open my $a, '<', $allf or die; while (<$a>) { chomp; s/\r//; my $p=$_; (my $b=$p)=~s#.*/##; $b=~s/\.m$//; push @{$cur{$b}}, $p } close $a;
open my $y, '<', $f25 or die; while (<$y>) { chomp; s/\r//; my $p=$_; (my $b=$p)=~s#.*/##; $b=~s/\.m$//; push @{$y25{$b}}, $p } close $y;
print join("\t", qw(name path new 2025_file current_same_name)), "\n";
for my $b (sort keys %cur) {
  my $new;
  if    ($b =~ /^settings_sopvar_(\w+)$/) { $new = "settings_PIETOOLS_$1" }
  elsif ($b eq 'eq_opts_sopvar')          { $new = 'eq_opts' }
  elsif ($b eq 'exec_sop_settings')       { $new = 'exec_settings' }
  elsif ($b =~ /^(\w+)_sop$/)             { $new = $1 }
  elsif ($b =~ /^(\w*)_sop_(\w+)$/)       { $new = "$1_$2" }
  else { next }
  for my $p (@{$cur{$b}}) {
    my $y = join(' | ', @{$y25{$new} // []}); my $c = join(' | ', @{$cur{$new} // []});
    print join("\t", $b, $p, $new, $y || '-', $c || '-'), "\n";
  }
}
