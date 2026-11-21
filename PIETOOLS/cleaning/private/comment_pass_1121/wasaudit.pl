#!/usr/bin/perl
# wasaudit.pl FILE... - the (was) lines of HEAD that the pass must keep (their code is code or
# commented code of the 2025 version) and that the working copy no longer holds (compared by
# their code part). Prints "LOST <HEAD line no>| <line>" per lost line, or "OK  FILE".
use strict; use warnings;
my $REPO = '/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS';
my $S    = '/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/cpass1121';
my %p25;
open(my $c, '<', "$S/census.tsv") or die; while (<$c>) { chomp; my @a = split /\t/; $p25{$a[0]} = $a[1] } close $c;
sub code {
    my $s = shift; my ($out, $prev, $i, $n) = ('', '', 0, length $s);
    while ($i < $n) {
        my $ch = substr($s, $i, 1);
        last if $ch eq '%';
        if ($ch eq '.' && substr($s, $i, 3) eq '...') { $out .= '...'; last }
        if ($ch eq '"' || ($ch eq "'" && $prev !~ /[A-Za-z0-9_)\]}.']/)) {
            my $q = $ch; my $j = $i + 1;
            while ($j < $n) { my $d = substr($s, $j, 1);
                if ($d eq $q) { if (substr($s, $j + 1, 1) eq $q) { $j += 2; next } last } $j++ }
            $out .= substr($s, $i, $j - $i + 1); $prev = '"'; $i = $j + 1; next;
        }
        $out .= $ch; $prev = ($ch =~ /\s/) ? ' ' : $ch; $i++;
    }
    $out =~ s/\s+/ /g; $out =~ s/^ | $//g; return $out;
}
sub uncomment { my $y = shift; 1 while $y =~ s/^\s*%+\s?//; return code($y) }
my $DATE = qr{\d{1,2}/\d{1,2}/\d{2,4}};
my $INI  = qr{[A-Z]{2,3}(?:\s*\(CC\))?(?:\s*,\s*[A-Z]{2,3}(?:\s*\(CC\))?)*};
for my $f (@ARGV) {
    my @old  = map { s/\r?\n$//r } `git -C $REPO show "origin/PIETOOLS-2025:PIETOOLS/$p25{$f}"`;
    my @head = map { s/\r?\n$//r } `git -C $REPO show "HEAD:PIETOOLS/$f"`;
    my @work = do { open(my $h, '<:raw', "$REPO/PIETOOLS/$f") or die; map { s/\r?\n$//r } <$h> };
    my %R = map { $_ => 1 } @old;
    my %C; for (@old) { my $x = code($_); $C{$x} = 1 if $x ne ''; if (/^\s*%/) { my $y = uncomment($_); $C{$y} = 1 if $y ne '' } }
    my %W; for (@work) { if (/^(.*?\S)\s*%\s*$INI\s*,?\s*$DATE.*\(was\)/) { $W{uncomment($1)}++ } elsif (/^\s*%/) { $W{uncomment($_)}++ } }
    my @lost;
    for my $k (0 .. $#head) { my $t = $head[$k]; next if $R{$t};
        next unless $t =~ /^(.*?\S)\s*(%\s*$INI\s*,?\s*$DATE[^%]*\(was\).*)$/;
        my $cp = uncomment($1); next if $cp eq q{} || !$C{$cp};
        next unless $cp =~ /[=;()\[\]{}]/ || $cp =~ /^\w+$/;   # code, not prose of a revised help text
        if (($W{$cp} // 0) > 0) { $W{$cp}--; next }
        push @lost, sprintf("LOST %4d| %s", $k + 1, $t) }
    print @lost ? join("\n", "in $f:", @lost) . "\n" : "OK  $f\n";
}
