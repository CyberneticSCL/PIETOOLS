Pre-change baseline of every file touched by the unpushed ndopvar work
======================================================================
Taken 2026-09-11.

WHAT THIS IS
  The files AS THEY STOOD BEFORE any of the unpushed work, so the original
  code can be recovered without needing git. The baseline is

      origin/ndopvar = e3a29ef0

  which is the last pushed commit. For the 9 files changed by the unpushed
  commits that is their state before those commits; for the 4 files modified
  only in the working tree, it is their state before those edits, since they
  were never committed.

  Current branch tip when this was taken: 420cd624 (9 commits ahead, plus
  those 4 uncommitted files).

  This does NOT store the modified versions. Those are the live working tree,
  and are in git history, so copying them here added nothing.

WHY IT LIVES HERE
  In the repository, so collaborators get it on a fetch -- but under a
  'private' folder, because 'pietools_path_update' is just
    addpath(genpath(pwd))
  and genpath SKIPS directories named 'private', '@...', '+...' and
  'resources'. It does NOT skip dot-directories, so a hidden folder would
  still have landed on the path; 'private' is the reliable choice. This was
  tested, not assumed.

  Three independent layers keep it off the path and inert:
    1. genpath never descends into 'private', so nothing here is on the path.
    2. Everything sits in a SUBfolder of 'private', so no file here becomes a
       private function of sopvar/ (only direct children of a private folder
       do that).
    3. Every source is stored as *.m.bak, so even if this folder were added
       to the path by hand, nothing is callable and nothing could shadow the
       live code. This also matters because the tree below preserves paths
       and so contains directories named @sopvar, @sdopvar and @mopvar; with
       no .m files in them they cannot register as duplicate class folders.

CONTENTS
  files_before_e3a29ef0/**.m.bak
        The 11 files that existed at e3a29ef0, repo-relative paths
        preserved. Each verified byte-for-byte against the blob at
        origin/ndopvar.
  files_before_e3a29ef0/NEW_FILES_NO_BASELINE.txt
        The 2 files created by the unpushed work, which therefore have no
        pre-change version:
            PIETOOLS/sopvar/misc/claude/canonical_var_order.m
            PIETOOLS/sopvar/misc/claude/monomial_gather.m
        To undo those, delete them.
  files_before_e3a29ef0/SHA256SUMS.txt     checksums of the 11 copies
  ndopvar-unpushed.bundle   git bundle of e3a29ef0..420cd624. Kept because
        the 9 commit messages carry the measurements behind each change
        (what came out bit-identical, what the end-to-end speedup was, at
        what variance) and plain file copies cannot preserve that.
  SHA256SUMS-bundle.txt     checksum of the bundle
  unpushed-commits.patch    combined before -> after diff of the 9 commits
  unpushed-commits.log      the 9 commit messages with per-file stats
  uncommitted-worktree.patch diff of the 4 uncommitted files only
  filelist.txt              all 13 repo-relative paths

THE 13 FILES
  with a pre-change version here (11):
    PIETOOLS/sopvar/@mopvar/mopvar.m                                    (uncommitted edit)
    PIETOOLS/sopvar/@sdopvar/mtimes.m
    PIETOOLS/sopvar/@sdopvar/mtimes_AT.m
    PIETOOLS/sopvar/@sdopvar/sdopvar.m
    PIETOOLS/sopvar/@sopvar/mtimes.m
    PIETOOLS/sopvar/@sopvar/private/leftShiftMonomials_SS.m
    PIETOOLS/sopvar/@sopvar/sopvar.m
    PIETOOLS/sopvar/@sopvar/sopvar2opvar2d.m                            (uncommitted edit)
    PIETOOLS/sopvar/Testfolder/sdopvar/claude_tests/test_canonical_multiplier.m (uncommitted edit)
    PIETOOLS/sopvar/Testfolder/sdopvar/rand_sopvar.m                    (uncommitted edit)
    PIETOOLS/sopvar/int_semisep.m
  created by the unpushed work, no baseline (2):
    PIETOOLS/sopvar/misc/claude/canonical_var_order.m
    PIETOOLS/sopvar/misc/claude/monomial_gather.m

THE COMMITS THIS IS THE BASELINE FOR
  The work was originally 11 commits and was then consolidated into 4
  before pushing, so the hashes below DO NOT resolve on the branch any
  more. They are listed because ndopvar-unpushed.bundle in this folder
  contains exactly them, with their original messages:

      420cd624  Vectorize int_2b's C3a assembly over the monomial-pair grid
      9d5fdba8  Memoise int_semisep's per-key results on the eight distinct factors
      22880a96  Build int_semisep's packed outputs from triplets, not by index assignment
      2f548d5b  Hoist the basis lift out of leftShiftMonomials_SS's coefficient loops
      28e0bd51  Sum the nMid diagonal blocks in mtimes_AT's PART 3
      516e502a  Sum the nMid diagonal blocks in mtimes instead of masking a full product
      bf4beffe  Keep the zero-block parameter shorthand working through the reorder
      f226ae30  Fix the live sdopvar composition path and harden canonical_var_order
      f40509e7  Align the right operand of an sopvar composition by variable name

  The 4 consolidated commits reproduce the identical final tree
  (55a8437cb32bb22f8db29e6af2840cdb1e2ec52b), so nothing in the content
  changed; only the grouping did. The baseline e3a29ef0 is unaffected,
  which is why every file copy here is still correct.

  The bundle plus the two backup commits (65c97468, 5a368623) are also
  reachable from the local tag 'claude-pre-squash' until it is deleted.

REVERTING TO THE BASELINE
  One file, from the repo root:
      cp PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0/<path>.m.bak <path>
  All 11, and remove the 2 new ones (bash, from the repo root). This
  overwrites live files -- be sure first:
      D=PIETOOLS/sopvar/private/backup_unpushed_2026-09-11/files_before_e3a29ef0
      ( cd "$D" && find . -name '*.m.bak' ) | sed 's|^\./||; s|\.bak$||' \
        | while read -r f; do cp "$D/$f.bak" "$f"; done
      rm -f PIETOOLS/sopvar/misc/claude/canonical_var_order.m \
            PIETOOLS/sopvar/misc/claude/monomial_gather.m
  Via git instead, which is cleaner if it is available:
      git checkout origin/ndopvar -- <path>

VERIFYING
      cd files_before_e3a29ef0 && sha256sum -c SHA256SUMS.txt   # 11x ": OK"
      sha256sum -c SHA256SUMS-bundle.txt
      git bundle verify ndopvar-unpushed.bundle
  All three passed at creation, and each .m.bak was additionally compared
  byte-for-byte against its origin/ndopvar blob.

NOT COVERED
  stash@{0} ("On master: !!GitHub_Desktop<master>") is excluded by request.
  It is untouched and still in the repo.
