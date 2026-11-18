# Fixed-point A/B runs, 11/17/2026 (DECISIONS D92-D94)

Scripts run from a harness folder (`cd <harness>; addpath(<this folder>); <script>('new'|'head')`):
- `ab_dde.m`: the light tier at four delays HEAD certifies and the working tree refuses, one
  tier knob varied per variant (`ab_dde_{new,head}.log`); V8 (lpi.var.deg{1} = [4 2 2]) is the
  variant that certifies. V1/V1b/V2/V3 on the new tree are invalid (nested cells from setfield).
- `ab_2d.m`: PDE 39 light, one solve per variant, no raise (`ab_2d_{new,head}.log`).
- `ab_dde_fix.m` with `ab_bisect.m` (lb_bisect of lb_run.m, extracted): the ten 1-D bisection
  rows under light with var.deg{1} = [4 2 2] (`ab_dde_fix_new.log`).
- `head_light_tier.m`: settings_sopvar_light.m of 33d30877, renamed (unused in the end).
- `run_ab.sh`: the slot-guarded runner (two MATLAB processes at most, watchdog).
Logs keep the `@@AB` result lines and the executive's `dw k:` classification lines.
MMP, 11/17/2026.
