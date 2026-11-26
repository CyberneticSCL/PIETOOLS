#!/usr/bin/perl
# Classify each hit line: area (top folder group), file extension, code vs comment
# (token after the first % that is outside a quoted string = comment), and collect
# every identifier that contains copvar/cdopvar.
use strict; use warnings;
my (%area_files,%area_lines,%area_code,%area_cmt,%area_tok,%ident,%identfiles,%ext);
while (my $l = <STDIN>) {
  chomp $l; $l =~ s/\r$//;
  my ($path,$ln,$txt) = $l =~ /^([^:]+):(\d+):(.*)$/ or next;
  my $area;
  if    ($path =~ m{^PIETOOLS/cleaning/})        { $area = 'cleaning' }
  elsif ($path =~ m{^PIETOOLS/orphaned_file/})   { $area = 'orphaned_file' }
  elsif ($path =~ m{^PIETOOLS/testing_codes/})   { $area = 'testing_codes' }
  elsif ($path =~ m{^PIETOOLS/SOSTOOLS400/})     { $area = 'SOSTOOLS400' }
  elsif ($path =~ m{^PIETOOLS/([^/]+)/})         { $area = "toolbox:$1" }
  elsif ($path =~ m{^PIETOOLS/})                 { $area = 'toolbox:(root)' }
  else                                           { $area = "repo:$path" }
  my ($e) = $path =~ /\.([^.\/]+)$/; $e //= '(none)';
  $ext{"$area\t$e"}{$path}=1;
  $area_files{$area}{$path} = 1; $area_lines{$area}++;
  # find first comment % outside single-quoted strings (MATLAB); for non-.m files the whole line counts as text
  my $cpos = -1;
  if ($e eq 'm') {
    my $inq = 0; my $prev = '';
    for (my $i=0;$i<length $txt;$i++){
      my $c = substr($txt,$i,1);
      if ($c eq "'" ) {
        if ($inq) { $inq = 0 }
        elsif ($prev !~ /[\w\)\]\}\.'"]/ ) { $inq = 1 }   # quote opens a string unless it is a transpose
      } elsif ($c eq '"' && !$inq) { my $j = index($txt,'"',$i+1); $i = $j if $j>0; }
      elsif ($c eq '%' && !$inq) { $cpos = $i; last }
      $prev = $c if $c ne ' ';
    }
  }
  my $tpos = -1; if ($txt =~ /c(d)?opvar/ig) { $tpos = pos($txt) - length($&) }
  if ($e eq 'm' && ($cpos < 0 || $tpos < $cpos)) { $area_code{$area}++ } else { $area_cmt{$area}++ }
  pos($txt) = 0; while ($txt =~ /([A-Za-z0-9_]*c(?:d)?opvar[A-Za-z0-9_]*)/ig) { my $id=$1; $ident{$id}++; $identfiles{$id}{$path}=1; $area_tok{$area}++ }
}
open my $A, '>', $ARGV[0] or die; print $A "area\tfiles\tlines\tcode_lines(.m)\tcomment_or_text_lines\ttokens\n";
for my $a (sort keys %area_files) { printf $A "%s\t%d\t%d\t%d\t%d\t%d\n",$a,scalar keys %{$area_files{$a}},$area_lines{$a},$area_code{$a}//0,$area_cmt{$a}//0,$area_tok{$a}//0 }
close $A;
open my $I, '>', $ARGV[1] or die; print $I "identifier\toccurrences\tfiles\n";
for my $id (sort { $ident{$b} <=> $ident{$a} } keys %ident) { printf $I "%s\t%d\t%d\n",$id,$ident{$id},scalar keys %{$identfiles{$id}} }
close $I;
open my $E, '>', $ARGV[2] or die; for my $k (sort keys %ext) { printf $E "%s\t%d\n",$k,scalar keys %{$ext{$k}} } close $E;
