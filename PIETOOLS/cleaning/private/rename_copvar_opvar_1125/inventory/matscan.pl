# Report which MAT-files (v5/v7) contain class-name byte strings in any data element,
# zlib-compressed elements (type 15) inflated first; MCOS class names live in subsystem data.
use strict; use warnings; use Compress::Zlib;
for my $path (@ARGV) {
  open my $fh, '<:raw', $path or die; local $/; my $d = <$fh>; close $fh;
  my %f; my $pos = 128;
  while ($pos + 8 <= length $d) {
    my ($typ,$nb) = unpack('VV', substr($d,$pos,8)); my $body;
    if ($typ >> 16) { $body = substr($d,$pos+4,4); $pos += 8; }
    else { $body = substr($d,$pos+8,$nb); $pos += 8 + $nb; $pos += (8 - $nb % 8) % 8 if $typ != 15; }
    if ($typ == 15) { $body = uncompress($body) // ''; }
    while ($body =~ /(?<![A-Za-z0-9_])(copvar|cdopvar|opvar|dopvar|sopvar|sdopvar|opvar2d|pde_struct|pie_struct|terms|sys|state)(?![A-Za-z0-9_])/g) { $f{$1}++ }
  }
  (my $s = $path) =~ s#.*?/PIETOOLS/PIETOOLS/##;
  print "$s\t", (%f ? join(',', map {"$_=$f{$_}"} sort keys %f) : '-'), "\n";
}
