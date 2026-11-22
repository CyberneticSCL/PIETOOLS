# orphaned_file/claude/

The markdown files of the toolbox that Claude sessions wrote, moved here on 11/22/2026 (maintainer:
"move all the .md files into orphaned_code/claude"). The path below this folder is the original path
below `PIETOOLS/`:

| File | Content |
|---|---|
| `PIETOOLS_demos/cuadmm/README.md` | the cuADMM bridge and its baseline suite |
| `PIETOOLS_demos/lowrank_2d_stability/README.md`, `GUIDE.md`, `ALGORITHM.md` | the low-rank 2-D stability package |
| `PIETOOLS_demos/sopvar_demos/README.md` | the heatNd benchmark |
| `executives_sopvar/README.md` | the container executives |
| `lpi_programming_sopvar/README.md` | the container LPI functions |
| `settings_sopvar/README.md` | the settings tiers |

Three markdown files stay in place:
- the repository `README.md`, the landing page of the repository;
- `SOSTOOLS400/README.md`, the README of the SOSTOOLS toolbox (2021);
- `.claude/agents/pietools-applications.md`, the definition of an agent type, read by Claude Code.

The move is committed as renames, so `git log --follow` traces each file to its original path.
`executives_sopvar/README.md` and `settings_sopvar/README.md` carry the edits of 11/22/2026 that removed
pointers into `cleaning/` and `testing_codes/`.
