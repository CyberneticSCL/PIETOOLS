#!/usr/bin/perl
# Applies h2_edits.txt, 11/26/2026. Each block "@@ path:N" names line N of a file; the "- " line
# must equal that line exactly (end of line excluded), and the "+ " lines replace it. All blocks
# of all files are checked before any file is written. New lines take the end of line of the
# line they replace, so a CRLF file stays CRLF and an LF file LF.
# Usage: perl h2_apply.pl <PIETOOLS folder> <h2_edits.txt> check|apply
use strict; use warnings;
my ($root,$tab,$mode) = @ARGV; die "mode check|apply\n" unless ($mode//'') =~ /^(check|apply)$/;
open(my $T,'<',$tab) or die "$tab: $!"; my (%blk,$cur);
while (my $l=<$T>) { $l =~ s/\r?\n$//;
  next if $l =~ /^#/ || $l eq '';
  if ($l =~ /^@@ (.+):(\d+)$/) { $cur = {f=>$1,n=>$2,old=>undef,new=>[]}; push @{$blk{$1}}, $cur; next }
  die "line outside a block: $l\n" unless $cur;
  if ($l =~ /^- (.*)$/s)  { die "two - lines at $cur->{f}:$cur->{n}\n" if defined $cur->{old}; $cur->{old}=$1; next }
  if ($l =~ /^\+ (.*)$/s) { push @{$cur->{new}}, $1; next }
  die "bad line: $l\n" }
close $T;
my ($nerr,$nblk,%out) = (0,0);
for my $f (sort keys %blk) {
  open(my $h,'<:raw',"$root/$f") or do { print "MISSING $f\n"; $nerr++; next };
  local $/; my $s = <$h>; close $h;
  my @L; while ($s =~ /\G([^\r\n]*)(\r\n|\n|\z)/gc) { push @L,[$1,$2]; last if $2 eq '' }
  pop @L if @L && $L[-1][0] eq '' && $L[-1][1] eq '';
  my %seen;
  for my $b (sort { $b->{n} <=> $a->{n} } @{$blk{$f}}) {
    $nblk++;
    die "two blocks on $f:$b->{n}\n" if $seen{$b->{n}}++;
    die "no + lines at $f:$b->{n}\n" unless @{$b->{new}};
    my $i = $b->{n}-1;
    if ($i > $#L || $L[$i][0] ne $b->{old}) {
      printf "MISMATCH %s:%d\n  want: %s\n  have: %s\n", $f, $b->{n}, $b->{old}, ($i<=$#L ? $L[$i][0] : '(past end)');
      $nerr++; next }
    my $eol = $L[$i][1] ne '' ? $L[$i][1] : ($s =~ /\r\n/ ? "\r\n" : "\n");
    my @rep = map { [$_,$eol] } @{$b->{new}};
    $rep[-1][1] = $L[$i][1];          # the last new line keeps the original end of line
    splice @L, $i, 1, @rep }
  $out{$f} = join('', map { $_->[0].$_->[1] } @L);
  printf "%-6s %s (%d blocks)\n", 'ok', $f, scalar @{$blk{$f}} }
print "blocks $nblk, errors $nerr\n";
exit 1 if $nerr;
if ($mode eq 'apply') { for my $f (sort keys %out) { open(my $o,'>:raw',"$root/$f") or die "$f: $!"; print $o $out{$f}; close $o } print "written ", scalar(keys %out), " files\n" }
