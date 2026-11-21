#!/usr/bin/perl
# restore_was.pl FILE TAG PAIRS - puts back (was) lines of 2025 commented code that the first
# version of cpass.pl apply dropped. PAIRS: "w1:r1,w2:r2,..." - HEAD line numbers of each
# (was) line and of the line it stands above. Each (was) line is reinserted, retagged
# "% TAG (was)" at column 77, right above the working-copy line whose text, marker removed,
# equals that of HEAD line r (unique match required). Keeps the line ending of that line.
use strict; use warnings;
my ($f, $tag, $pairs) = @ARGV;
my $REPO = '/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS';
my @head = map { s/\r?\n$//r } `git -C $REPO show "HEAD:PIETOOLS/$f"`;
open(my $h, '<:raw', "$REPO/PIETOOLS/$f") or die; my @work = <$h>; close $h;
my $MK = qr{\s*%\s*[A-Z]{2,3}(?:\s*,\s*[A-Z]{2,3})*\s*,?\s*\d{1,2}/\d{1,2}/\d{2,4}[^%]*$};
sub bare { my $t = shift; $t =~ s/\r?\n$//; $t =~ s/$MK//; $t =~ s/\s+$//; return $t }
my @ins;   # [working index, new line]
for my $p (split /,/, $pairs) {
    my ($w, $r) = split /:/, $p;
    my $cp = bare($head[$w - 1]);
    my $target = bare($head[$r - 1]);
    my @hit = grep { bare($work[$_]) eq $target } 0 .. $#work;
    die "no unique match for HEAD line $r ($target): @hit\n" unless @hit == 1;
    my ($eol) = $work[$hit[0]] =~ /(\r?\n)$/;
    my $line = length($cp) < 76 ? $cp . (' ' x (76 - length $cp)) . "% $tag (was)" : "$cp % $tag (was)";
    push @ins, [$hit[0], $line . $eol];
}
for my $i (sort { $b->[0] <=> $a->[0] } @ins) { splice @work, $i->[0], 0, $i->[1] }
open(my $o, '>:raw', "$REPO/PIETOOLS/$f") or die; print $o @work; close $o;
print "restore_was $f: ", scalar(@ins), " lines put back\n";
