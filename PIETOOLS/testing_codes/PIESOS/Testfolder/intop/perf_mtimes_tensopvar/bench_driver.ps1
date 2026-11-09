param(
    [Parameter(Mandatory=$true)][string]$Tag,          # prefix for output files
    [Parameter(Mandatory=$true)][string]$SpecsJson,    # JSON array of spec objects
    [int]$MaxMB = 5000,
    [int]$TimeoutSec = 600
)
# One fresh guarded MATLAB process per spec; collects result JSONs into $Tag.results.jsonl
$sp = Split-Path -Parent $MyInvocation.MyCommand.Path
$specs = $SpecsJson | ConvertFrom-Json
$out = Join-Path $sp "$Tag.results.jsonl"
if (Test-Path $out) { Remove-Item $out }
$i = 0
foreach ($s in $specs) {
    $i++
    $specFile = Join-Path $sp "$Tag.spec$i.json"
    $resFile  = Join-Path $sp "$Tag.res$i.json"
    $logFile  = Join-Path $sp "$Tag.log$i.txt"
    ($s | ConvertTo-Json -Depth 5) | Out-File -Encoding ascii $specFile
    if (Test-Path $resFile) { Remove-Item $resFile }
    $env:BENCH_SPEC = $specFile
    $env:BENCH_OUT  = $resFile
    $g = & (Join-Path $sp "run_guarded.ps1") -Launcher (Join-Path $sp "bench_run.m") -Log $logFile -MaxMB $MaxMB -TimeoutSec $TimeoutSec
    $line = [ordered]@{ idx = $i; spec = $s; guard = "$g" }
    if (Test-Path $resFile) {
        $line.result = (Get-Content $resFile -Raw | ConvertFrom-Json)
    } else {
        $line.result = $null
    }
    ($line | ConvertTo-Json -Depth 12 -Compress) | Add-Content -Encoding ascii $out
    $r = $line.result
    if ($r -and $r.ok) {
        "{0,3} N={1} d={2} q={3,-8} t_call={4,8:N3}s t_min={5,8:N3}s peak_before={6,6:N0} peak_after={7,6:N0} MB nnzA={8} nnzC={9} | {10}" -f $i,$s.N,$s.d,$s.q,$r.t_call,$r.t_min,$r.peak_before,$r.peak_after,$r.nnzA,$r.nnzC,$g
    } else {
        "{0,3} N={1} d={2} q={3,-8} FAILED/KILLED: {4} | {5}" -f $i,$s.N,$s.d,$s.q,$(if($r){$r.err}else{'(no result)'}),$g
    }
}
