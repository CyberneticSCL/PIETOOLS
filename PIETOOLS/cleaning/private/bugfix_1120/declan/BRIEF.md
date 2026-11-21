# Report for Declan, 11/21/2026 – brief for every agent

Maintainer, 11/21/2026: "draft the report for Declan". The report lists the defects that the bug-fix
pass of 11/20 held because the code is Declan Jagt's (DJ), or lies in PIESOS, which he maintains. The
items are in `cleaning/private/release_comment_pass/CANDIDATES.md`: read its header and its outcome
section. Its line numbers refer to 78641839; HEAD is now fa0deda4, and two commits of Declan came after
78641839: 0a1ed7f0 (algebraic constraints, 21 files) and 5b395a21 (pie_struct for non-square T, 6 files).
Paths are relative to the PIETOOLS folder `C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS`
(git paths carry a `PIETOOLS/` prefix). The classes opvar, opvar2d, dopvar and dopvar2d no longer exist
in the repository (retired 11/08/2026).

Records of the bug-fix pass: `cleaning/private/bugfix_1120/` – BRIEF.md (its section "MATLAB" holds the
lock runner and the run rules, which apply here unchanged) and the ledgers A_sopvar (C101), D_piesim
(C53, C58), E_examples (C84, C92) and the outcome section of CANDIDATES.md (C115).

## Rules

- Edit no repository file. Write only your ledger `cleaning/private/bugfix_1120/declan/<agent>_ledger.md`
  and files in your scratch folder
  `/c/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/<agent>/`.
  No commit, stash, reset, checkout or restore. Other sessions have uncommitted edits in this checkout;
  leave them alone.
- Never run code that changes state outside the MATLAB session: `pietools_update` (mkdir, startup.m,
  savepath), `savepath`, anything that writes outside your scratch folder. For such an item, run a copy
  in your scratch folder with the side-effecting calls replaced by stubs, or argue from the MATLAB rule.
- MATLAB only through the lock runner (`bash .../scratchpad/mlrun.sh "run('<script>')"`, Bash timeout
  600000; exit 75 = lock busy, retry). At most one batch MATLAB at a time across all agents; every run
  under 5 minutes; no SDP beyond a few seconds.

## Per item

1. Locate the cited lines at HEAD. If the defect is gone at HEAD, record ALREADY FIXED with the commit.
2. Author. Establish who wrote the cited lines: `git blame -w -C -C -M`, `git log -S'<text>'`,
   `git log -L`. Plain blame names Declan for lines he only imported or moved: the "PIETOOLS 2021b
   Release." commits of 12/31/2021 (e15f2543, 0af98826, ...), 4cefdf8f "Move old initializers/converters
   to legacy" (2022-08-18), "Add @pde_struct/convert" (2022-09-27). Record author, commit and date. If
   the lines are not by Declan, record NOT DJ with the evidence (except in PIESOS, see dS below).
3. Reproduce. For a reachable path, write a minimal script `<item>_repro.m` in your scratch folder. It
   starts with `cd('C:\Users\mpeet\OneDrive\Documents\GitHub\PIETOOLS\PIETOOLS'); pietools_path_update;`,
   asserts with `which` that each function under test resolves in the working tree, and writes what it
   shows (wrong value next to the right one, the error message, the missing output) to
   `<item>_repro.txt` beside it, so that the coordinator can rerun it unchanged. For a semantics-only
   item (dead code, unreachable branch, text), state the MATLAB rule or the text with file:line; no
   run. A hedged item ("may", "unless", "not determined by reading") is a hypothesis: if it does not
   reproduce, record NOT REPRODUCED. Keep measured and inferred apart.
4. Consequence class: silent wrong result / error on a reachable path / error on a rare or unreachable
   path / dead code or no effect / misleading text or message. Say how a user reaches the path.
5. Draft the report entry (format below). Suggest a fix only when it follows from the code, its help
   or the user manual (`PIETOOLS_documentation/`), in one line; otherwise none.

## Report entry format

```
**C26** `parser/@pde_struct/expand_delays.m:547`, local add_delay_state (DJ, 1a2b3c4d, 2024-05-01) –
<what happens, as measured>. <consequence and how a user reaches it>. Fix: <one line, optional>.
Reproduction: <at most 6 lines of MATLAB, or "by reading: <rule>">
```

Wording: terse and factual; no "I", no possessives ('s), no "we"; spaced en dashes ` – `, never em
dashes; name no retired class as if it existed. Address nobody: the maintainer forwards the report.

## Ledger

Per item: status (CONFIRMED (run) / CONFIRMED (reading) / NOT REPRODUCED / ALREADY FIXED / NOT DJ),
HEAD location, author evidence, reproduction (script, measured output), consequence, report entry.
Update it after each item, so that an interruption loses at most one item.

Final message: one line per item (status and class), then the list of scripts; at most 40 lines.

## Assignments

- **dP** (parser, GUI, root scripts): C13, C18, C19, C20, C24, C25, C26, C27, C28, C29, C30, C33, C34.
- **dC** (converters, PIESIM, sopvar messages): C37, C41, C42, C43, C45, C47; C53 part (a) and C58 (D_piesim
  ledger); C101 (A_sopvar ledger: authors and the proposed wording).
- **dS** (PIESOS): C62, C65, C72 (DJ); C63, C64, C66, C67, C68, C69, C70, C71, C73, C111 (SS, AT, CR).
  PIESOS is maintained by Declan: for the items of other authors record the author, and the status
  CONFIRMED / NOT REPRODUCED / ALREADY FIXED as usual (never NOT DJ).
- **dE** (examples, demos): C78, C80, C81, C82, C83, C88, C96, C97, C98, C99, C100; C84 (the Saba file),
  C92 (line 62 of CL/stability_8.m and the Das files), C115 (E_examples ledger, outcome section).
