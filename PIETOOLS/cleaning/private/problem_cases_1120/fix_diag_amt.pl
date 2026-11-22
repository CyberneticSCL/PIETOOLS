# fix_diag_amt.pl IN > OUT: diag_amt also records exst and eta of lb_exec (header log, (was) lines,
# column-77 markers on the last line of each statement). Keeps the line ending of the file.
use strict; use warnings;
my ($h, $a, $b) = (0, 0, 0);
my $hdr_old = q{    fprintf(fid,'key,form,variant,tier,mods,status,fr,pinf,dinf,numerr,ndec,m,sec,note\n');  fclose(fid);};
my $hdr_new = q{    fprintf(fid,'key,form,variant,tier,mods,status,fr,pinf,dinf,numerr,ndec,m,sec,note,exst,eta\n');  fclose(fid);};
my $fmt_old = q{    fprintf(fid,'%s,%s,%s,%s,"%s",%s,%.4g,%g,%g,%g,%.0f,%.0f,%.1f,"%s"\n',key,form,[variant,sprintf(stt,v)],tier, ...};
my $fmt_new = q{    fprintf(fid,'%s,%s,%s,%s,"%s",%s,%.4g,%g,%g,%g,%.0f,%.0f,%.1f,"%s",%s,%.3g\n',key,form,[variant,sprintf(stt,v)],tier, ...};
my $row_old = q{        g(r,'sec'),note);};
my $row_new = q{        g(r,'sec'),note,f(r,'exst'),g(r,'eta'));};
sub mark { my ($code, $tag) = @_; my $m = "% MMP, 11/20/2026$tag"; return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
sub was  { (my $s = shift) =~ s/^(\s*)/$1% /; return $s; }
my @L = <>;
for (my $i = 0; $i < @L; $i++) {
    my ($body, $eol) = $L[$i] =~ /^(.*?)(\r?\n)$/ ? ($1, $2) : ($L[$i], '');
    if (!$h && $body eq '% MMP, 11/19/2026.') {
        print $body, $eol;
        print "% MMP, 11/20/2026: each row also records the executive status and the residual$eol";
        print "%                   eta of lb_exec (exst, eta): the verdict of the executive$eol";
        print "%                   flags beside the harness gate on the optimality measure.$eol";
        $h = 1; next;
    }
    if (!$a && $body eq $hdr_old) {
        print mark(was($hdr_old), ' (was)'), $eol, mark($hdr_new, ''), $eol; $a = 1; next;
    }
    if (!$b && $body eq $fmt_old) {         # a three-line statement: (was) block, then the new one
        my ($mid) = $L[$i+1] =~ /^(.*?)\r?\n$/;
        my ($end) = $L[$i+2] =~ /^(.*?)\r?\n$/;
        die "row end not found\n" unless $end eq $row_old;
        print was($fmt_old), $eol, was($mid), $eol, mark(was($row_old), ' (was)'), $eol;
        print $fmt_new, $eol, $mid, $eol, mark($row_new, ''), $eol;
        $i += 2; $b = 1; next;
    }
    print $L[$i];
}
die "anchors h=$h a=$a b=$b\n" unless $h && $a && $b;
