# pc_lib39_edit.pl FILE STATEMENT...
# Worktree copy of the PDE 39 example (PIETOOLS_PDE_Ex_2D_NS_Antonelli.m) only: sets the stability form to 'P'
# (old line kept, commented, as (was)) and adds the STATEMENTs after the line st = lpisettings_sop('light');
# so libcheck runs the file with the setting certified in the harness. Keeps the CRLF endings of the file.
use strict; use warnings;
my ($f, @stm) = @ARGV;
open(my $in, '<', $f) or die "no $f\n"; my @L = <$in>; close $in;
my $eol = $L[0] =~ /\r\n$/ ? "\r\n" : "\n";
sub mark { my ($code, $sfx) = @_; my $m = "% MMP, 11/20/2026 (worktree test$sfx)";
           return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
my ($iform, $ist) = (-1, -1);
for my $i (0..$#L) {
    $iform = $i if $L[$i] =~ /^evalin\('base','stability_form = ''Q'';'\);/;
    $ist   = $i if $L[$i] =~ /^st = lpisettings_sop\('light'\);/;
}
die "form line not found\n" if $iform < 0;  die "st line not found\n" if $ist < 0;
splice(@L, $ist + 1, 0, map { mark($_, '') . $eol } @stm);
splice(@L, $iform, 1, mark("%evalin('base','stability_form = ''Q'';');", ', was') . $eol,
                      mark("evalin('base','stability_form = ''P'';');", '') . $eol);
open(my $w, '>', $f) or die; print $w @L; close $w;
print "edited $f: form line ", $iform + 1, ", st line ", $ist + 1, "\n";
