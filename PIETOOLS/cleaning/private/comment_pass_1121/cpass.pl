#!/usr/bin/perl
# cpass.pl MODE FILE [TAG] - helper of the comment pass of 11/21/2026 over files of the
# 2025 release. FILE is relative to the PIETOOLS folder; its 2025 version is read from
# origin/PIETOOLS-2025 (path from census.tsv). A line is "new" when no line of the 2025
# version equals it (CR stripped). TAG is the initials and date of the single note, e.g.
# "YP, DJ, MMP, 11/21/2026".
#   prep  FILE          report: authors, proposed TAG, new header entries, new markers with a
#                       verdict, new comment lines, the code diff against 2025
#   apply FILE TAG      drop new (was) lines whose code is not 2025 code; retag the other new
#                       markers to "% TAG" ("% TAG (was)"); markers at column 77
#   check FILE TAG      the rules of the pass; prints PASS or the failures
use strict; use warnings;
my ($mode, $f, $tag) = @ARGV;
my $REPO = '/c/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS';
my $S    = '/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/cpass1121';
my $wf   = $ENV{CPASS_WF} // "$REPO/PIETOOLS/$f";
my %ini  = ('MATTHEW M PEET' => 'MMP', 'mpeet' => 'MMP', 'Declan' => 'DJ', 'Yulia' => 'YP',
            'DanBraghini' => 'DB', 'Varshitha Janavi P' => 'VJ', 'Talitsky' => 'AT');

# 2025 path from the census
my $p25;
open(my $c, '<', "$S/census.tsv") or die "census: $!";
while (<$c>) { chomp; my @a = split /\t/; if ($a[0] eq $f) { $p25 = $a[1]; last } }
close $c;
die "not in census: $f\n" unless defined $p25;
my @old = map { s/\r?\n$//r } `git -C $REPO show "origin/PIETOOLS-2025:PIETOOLS/$p25"`;
my %R = map { $_ => 1 } @old;

# MATLAB code of a line with comments removed (port of codeident.awk), whitespace normalized
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
my %C = map { (code($_) ne '') ? (code($_) => 1) : () } @old;
# commented code of the 2025 version counts as 2025 code: it existed in the prior release
for (@old) { if (/^\s*%/) { my $y = $_; 1 while $y =~ s/^\s*%+\s?//; my $x = code($y); $C{$x} = 1 if $x ne '' } }

my $DATE = qr{\d{1,2}/\d{1,2}/\d{2,4}};
my $INI  = qr{[A-Z]{2,3}(?:\s*\(CC\))?(?:\s*,\s*[A-Z]{2,3}(?:\s*\(CC\))?)*};   # also "DJ (CC)"
my $ENTRY = qr{^%\s*(?:$INI\s*,?\s*$DATE\s*[:\-;,]|$DATE\s*-\s*$INI|$INI\s+$DATE\b)};
sub marker {   # trailing developer marker: (code part, marker) or ()
    my $l = shift;
    return ($1, $2) if $l =~ /^(.*?\S)\s*(%\s*$INI\s*,?\s*$DATE[^%]*(?:;\s*$INI\s*,?\s*$DATE[^%]*)*)$/;
    return ();
}
sub wascode { my $cp = shift; 1 while $cp =~ s/^\s*%+\s?//; return code($cp) }
# an entry dated before the 2025 release (12/29/2025) is a note of an earlier release, reworded
sub predates { my $t = shift; return 0 unless $t =~ m{(\d{1,2})/(\d{1,2})/(\d{2,4})};
               my ($m, $d, $y) = ($1, $2, $3); $y += 2000 if $y < 100;
               return sprintf('%04d%02d%02d', $y, $m, $d) lt '20251229' }
# initials|date of every header entry of the 2025 version: a new line with such a key is a 2025
# entry whose wording changed, not a new entry
sub ekey { my $t = shift; return ($t =~ /^%\s*($INI)\s*,?\s*($DATE)/) ? (($1 =~ s/\s//gr) . "|$2") :
                                 ($t =~ /^%\s*($DATE)\s*-\s*($INI)/) ? (($2 =~ s/\s//gr) . "|$1") : '' }
my %OLDKEY = map { my $e = ekey($_); ($e ne '') ? ($e => 1) : () } @old;

my @lines = do { open(my $h, '<:raw', $wf) or die "$wf: $!"; local $/; split /(?<=\n)/, <$h> };

if ($mode eq 'prep') {
    print "FILE $f\n2025 $p25\n";
    my @au = map { chomp; $_ } `git -C $REPO log --reverse --format=%an origin/PIETOOLS-2025..HEAD -- "PIETOOLS/$f"`;
    my (%seen, @ord); for (@au) { my $i = $ini{$_} // "?$_"; push @ord, $i unless $seen{$i}++ }
    my (%nseen, @nord);
    for my $l (@lines) { my $t = $l =~ s/\r?\n$//r; next if $R{$t};
        while ($t =~ /%\s*($INI)\s*,?\s*$DATE/g) { for (split /\s*,\s*/, $1 =~ s/\s*\(CC\)//gr) { push @nord, $_ unless $nseen{$_}++ } } }
    print "commit authors since 2025 (first change first): @ord\n";
    print "initials in new notes and markers: @nord\n";
    my %all = map { $_ => 1 } @ord; my @tag = @ord; for (@nord) { push @tag, $_ unless $all{$_}++ }
    print "proposed TAG: ", join(', ', grep { !/^\?/ } @tag), ", 11/21/2026\n\n";
    print "NEW HEADER ENTRIES (and following new comment lines):\n";
    for my $k (0 .. $#lines) { my $t = $lines[$k] =~ s/\r?\n$//r;
        next if $R{$t} || $t !~ $ENTRY;
        print sprintf("  %4d| %s\n", $k + 1, $t);
        for (my $j = $k + 1; $j <= $#lines; $j++) { my $u = $lines[$j] =~ s/\r?\n$//r;
            last if $R{$u} || $u !~ /^%/ || $u =~ $ENTRY; print sprintf("  %4d| %s\n", $j + 1, $u) } }
    print "\nNEW MARKED LINES (verdict):\n";
    for my $k (0 .. $#lines) { my $t = $lines[$k] =~ s/\r?\n$//r; next if $R{$t};
        my ($cp, $m) = marker($t); next unless defined $m;
        my $v = ($m =~ /\(was\)/) ? ($C{wascode($cp)} ? 'KEEP-WAS' : 'DROP-WAS') : 'RETAG';
        $v = 'REVIEW-BEGIN/END' if $t =~ /\b(BEGIN|END)\b/;
        print sprintf("  %-16s %4d| %s\n", $v, $k + 1, $t) }
    print "\nNEW COMMENT LINES WITHOUT A MARKER (internal notes?):\n";
    for my $k (0 .. $#lines) { my $t = $lines[$k] =~ s/\r?\n$//r; next if $R{$t} || $t =~ $ENTRY;
        next unless $t =~ /^\s*%/; my @m = marker($t); next if @m;
        print sprintf("  %4d| %s\n", $k + 1, $t) }
    print "\nCODE DIFF AGAINST 2025 (comments removed):\n";
    open(my $o1, '>', "$S/_old_$$.txt"); print $o1 map { (code($_) ne '') ? code($_) . "\n" : () } @old; close $o1;
    open(my $o2, '>', "$S/_new_$$.txt"); print $o2 map { my $t = s/\r?\n$//r; (code($t) ne '') ? code($t) . "\n" : () } @lines; close $o2;
    print `diff $S/_old_$$.txt $S/_new_$$.txt | head -120`;
}
elsif ($mode eq 'apply') {
    die "TAG needed\n" unless $tag;
    my ($drop, $re, $kw, @out) = (0, 0, 0);
    for my $l (@lines) {
        my ($t, $eol) = $l =~ /^(.*?)(\r?\n?)$/s;
        if (!$R{$t}) {
            my ($cp, $m) = marker($t);
            if (defined $m && $t !~ /\b(BEGIN|END)\b/) {
                if ($m =~ /\(was\)/) {
                    if (!$C{wascode($cp)}) { $drop++; next }
                    $t = "$cp % $tag (was)"; $kw++;
                } else { $t = "$cp % $tag"; $re++ }
                my ($a) = $t =~ /^(.*?\S)\s+% \Q$tag\E/;
                $t = (length($a) < 76 ? $a . (' ' x (76 - length $a)) : "$a ") . substr($t, index($t, "% $tag"));
            }
        }
        push @out, $t . $eol;
    }
    open(my $h, '>:raw', $wf) or die; print $h @out; close $h;
    print "apply $f: dropped $drop (was) lines, kept $kw (was) lines, retagged $re markers\n";
}
elsif ($mode eq 'check') {
    die "TAG needed\n" unless $tag;
    my @fail; my ($nent, $notelen, $inote) = (0, 0, 0);
    for my $k (0 .. $#lines) { my $t = $lines[$k] =~ s/\r?\n$//r;
        if ($inote) { if (!$R{$t} && $t =~ /^%\s{2,}\S/ && $t !~ $ENTRY) { $notelen++; next } $inote = 0 }
        next if $R{$t};
        if ($t =~ $ENTRY) {
            if ($t =~ /^% \Q$tag\E: Changes since PIETOOLS 2025: /) { $nent++; $notelen = 1; $inote = 1 }
            elsif ($t =~ /^% \Q$tag\E: (BEGIN|END)\b/) { }                 # block markers of a change
            elsif ($OLDKEY{ekey($t)} || predates($t)) { }                     # a 2025 entry, reworded
            else { push @fail, sprintf("line %d: new header entry other than the note: %s", $k + 1, $t) }
            next;
        }
        my ($cp, $m) = marker($t);
        if (defined $m) {
            push @fail, sprintf("line %d: marker is not the tag: %s", $k + 1, $m) unless $m =~ /^% \Q$tag\E( \(was\))?$/;
            push @fail, sprintf("line %d: (was) code not in 2025: %s", $k + 1, $cp) if $m =~ /\(was\)/ && !$C{wascode($cp)};
            my $col = index($t, $m) + 1; my $len = length($cp);
            push @fail, sprintf("line %d: marker at column %d", $k + 1, $col) unless $col == 77 || ($len >= 76 && $col == $len + 2);
        }
        push @fail, sprintf("line %d: internal wording? %s", $k + 1, $t)
            if $t =~ /%.*\b(maintainer|DECISIONS|D1\d\d\b|cleaning\/|ledger|session|Claude|agent|Sol job|comment pass)/i;
    }
    push @fail, "the single note occurs $nent times" unless $nent == 1;
    push @fail, "the note has $notelen lines (at most 9)" if $notelen > 9;
    print @fail ? join("\n", "CHECK $f:", @fail) . "\n" : "PASS  $f\n";
}
