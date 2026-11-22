# fix_st_run.pl IN > OUT: st_run tag 'EXVT3' (the TRADE3 cuts, plan pc_plan_trade3.csv), with the header log,
# (was) lines and markers. Keeps the line ending of the file.
use strict; use warnings;
my ($h, $c, $p, $d) = (0, 0, 0, 0);
my $case_old = q{    case {'EXV1','EXV2','EXV1B','EXVX'}  % verification of the planned st and forms of each example};
my $case_new = q{    case {'EXV1','EXV2','EXV1B','EXVX','EXVT3'}  % verification of the planned st and forms of each example};
my $pf_end   = q{                      'example_settings_1119','phaseC_plan.csv');};
my $dim_old  = q{        if any(strcmp(tag,{'EXV1','EXV1B'})),  ids = ids(dims==1);  opts.tiers2d = {};  else,  ids = ids(dims==2);  end};
my $dim_new  = q{        if any(strcmp(tag,{'EXV1','EXV1B','EXVT3'})),  ids = ids(dims==1);  opts.tiers2d = {};  else,  ids = ids(dims==2);  end};
sub mark { my ($code, $tag) = @_; my $m = "% MMP, 11/20/2026$tag"; return length($code) < 76 ? sprintf("%-76s%s", $code, $m) : "$code  $m"; }
sub was  { (my $s = shift) =~ s/^(\s*)/$1% /; return $s; }
while (my $l = <>) {
    my ($body, $eol) = $l =~ /^(.*?)(\r?\n)$/ ? ($1, $2) : ($l, '');
    if (!$h && $body eq '% MMP, 11/18/2026.') {
        print $l;
        print "% MMP, 11/20/2026: tag 'EXVT3' verifies the cuts of the trade pass TRADE3$eol";
        print "%                   (plan pc_plan_trade3.csv of this folder; one-variable rows).$eol";
        $h = 1; next;
    }
    if (!$c && $body eq $case_old) { print mark(was($case_old),' (was)'), $eol, mark($case_new,''), $eol; $c = 1; next; }
    if ($c && !$p && $body eq $pf_end) {
        print $l;
        print mark("        % EXVT3: the cuts of the trade pass TRADE3 of 11/20/2026", ''), $eol;
        print mark("        if strcmp(tag,'EXVT3'),  pf = fullfile(here,'pc_plan_trade3.csv');  end", ''), $eol;
        $p = 1; next;
    }
    if ($p && !$d && $body eq $dim_old) { print mark(was($dim_old),' (was)'), $eol, mark($dim_new,''), $eol; $d = 1; next; }
    print $l;
}
die "anchors h=$h c=$c p=$p d=$d\n" unless $h && $c && $p && $d;
