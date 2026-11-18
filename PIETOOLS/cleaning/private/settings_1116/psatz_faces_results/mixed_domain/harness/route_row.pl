#!/usr/bin/perl
# route_row.pl FILE: the exec_tools_sop route shadow with the storage of
# every dimension from the proposal's row (PF_ROW_STORAGE, row in appdata
# PF_ROW). In place. SCRATCH, MMP, 11/16/2026.
use strict; use warnings;
my ($f) = @ARGV;
open(my $fh, '<', $f) or die $!;  my $s = do { local $/; <$fh> };  close $fh;
my @edits = (
 [q{    tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom);},
  q{%   tl.storage = @(prog,X,side) poslpivar_settings_sop(prog,X,st,'lf',side,PIE.dom); % MMP, 11/16/2026 (was; route_row)
    tl.storage = @(prog,X,side) pf_row_storage(prog,X,side,struct('vars',{reshape(PIE.vars(:,1).varname,1,[])},'dom',PIE.dom),getappdata(0,'PF_ROW')); % MMP, 11/16/2026 (route_row)}],
 [q{    tl.storage = @(prog,X,side) pos_def(prog,X,side,tl.dom);                % MMP, 11/13/2026},
  q{%   tl.storage = @(prog,X,side) pos_def(prog,X,side,tl.dom);                % MMP, 11/16/2026 (was; route_row); MMP, 11/13/2026
    tl.storage = @(prog,X,side) pf_row_storage(prog,X,side,tl.dom,getappdata(0,'PF_ROW')); % MMP, 11/16/2026 (route_row)}],
 [q{    tl.storage = @(prog,X,side) pos_nd(prog,X,side,S,tl.dom);               % MMP, 10/29/2026},
  q{%   tl.storage = @(prog,X,side) pos_nd(prog,X,side,S,tl.dom);               % MMP, 11/16/2026 (was; route_row); MMP, 10/29/2026
    tl.storage = @(prog,X,side) pf_row_storage(prog,X,side,tl.dom,getappdata(0,'PF_ROW')); % MMP, 11/16/2026 (route_row)}],
);
for my $e (@edits) {
    my ($o, $n) = @$e;
    my $i = index($s, $o);
    die "not found: $o" if $i < 0;
    substr($s, $i, length($o)) = $n;
}
open($fh, '>', $f) or die $!;  print $fh $s;  close $fh;
print "route_row: 3 storage lines replaced\n";
