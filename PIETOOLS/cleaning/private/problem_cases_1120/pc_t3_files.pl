# pc_t3_files.pl PLAN VERIFY EXAMPLES_ROOT LIBCASES OUTCASES
# For the PDE example files among the TRADE3 cuts that every test verified (VERIFY: verify_EXVT3.csv), appends
# after the last st statement of the file (worktree copy only) the cut setting: st = lpisettings_sop(tier) and
# the statements, so libcheck runs the file with it. Writes the libcheck cases of these files to OUTCASES.
# Indexed files (Saba) are skipped: their st is set per case. Keeps the line ending of each file.
use strict; use warnings; use File::Find;
sub csv { my $l = shift; $l =~ s/\r?\n$//; return ($l =~ /("[^"]*"|[^,]*)(?:,|$)/g); }
my ($plan, $verify, $root, $libcases, $out) = @ARGV;
my %ok; open(my $v,'<',$verify) or die; <$v>; while (<$v>) { my @f = csv($_); $ok{$f[0]} = $f[5]; }
my %case; open(my $c,'<',$libcases) or die; <$c>; while (<$c>) { my @f = csv($_); $case{$f[1]} //= [@f]; }
my %path; find(sub { $path{$1} = $File::Find::name if /^(PIETOOLS_PDE_Ex_\w+)\.m$/ }, $root);
open(my $o,'>',$out) or die; print $o "id,file,index\n";
open(my $p,'<',$plan) or die; <$p>;
while (<$p>) {
    my @f = csv($_); my ($file, $rows, $tier, $stm) = @f[1,3,7,8]; $stm =~ s/"//g;
    next unless $file =~ /^PIETOOLS_PDE_Ex_/;
    my @ks = split / /, $rows; next if grep { ($ok{$_} // '') ne '1' } @ks;
    if ($file =~ /Saba/) { print STDERR "skip $file (per-case st)\n"; next; }
    my $fp = $path{$file} or die "no file $file\n";
    open(my $in,'<',$fp) or die; my @L = <$in>; close $in;
    my $last = -1; for my $i (0..$#L) { $last = $i if $L[$i] =~ /^\s*st(\s*=|\.)/; }
    die "no st in $file\n" if $last < 0;
    my $eol = $L[$last] =~ /\r\n$/ ? "\r\n" : "\n";
    splice(@L, $last+1, 0, "st = lpisettings_sop('$tier');  % worktree test 11/20/2026: the TRADE3 cut$eol", "$stm$eol");
    open(my $w,'>',$fp) or die; print $w @L; close $w;
    die "no library case for $file\n" unless $case{$file};
    print $o join(',', @{$case{$file}}[0..2]), "\n";
    print STDERR "set $file: $tier; $stm\n";
}
