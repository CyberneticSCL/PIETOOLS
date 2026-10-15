param(
    [Parameter(Mandatory=$true)][string]$Launcher,   # .m script to run via matlab -batch run('...')
    [Parameter(Mandatory=$true)][string]$Log,        # stdout log
    [int]$MaxMB = 5000,                              # kill MATLAB (this run only) above this working set
    [int]$TimeoutSec = 900                           # ... or after this many seconds
)
# Runs one MATLAB R2026b batch job, tracks its peak working set, and kills ONLY the processes
# it started if the cap is exceeded.
$p = Start-Process -FilePath "C:\Program Files\MATLAB\R2026b\bin\matlab.exe" `
    -ArgumentList @('-batch',"run('$Launcher')") `
    -RedirectStandardOutput $Log -RedirectStandardError "$Log.err" -PassThru -WindowStyle Hidden
$launcherPid = $p.Id
$sw = [Diagnostics.Stopwatch]::StartNew()
$childPid = $null
$peakMB = 0.0
$reason = 'finished'
while (-not $p.HasExited) {
    if (-not $childPid) {
        $c = Get-CimInstance Win32_Process -Filter "ParentProcessId=$launcherPid AND Name='MATLAB.exe'" -ErrorAction SilentlyContinue
        if ($c) { $childPid = ($c | Select-Object -First 1).ProcessId }
    }
    if ($childPid) {
        $pr = Get-Process -Id $childPid -ErrorAction SilentlyContinue
        if ($pr) {
            $pr.Refresh()
            $mb = $pr.PeakWorkingSet64 / 1MB
            if ($mb -gt $peakMB) { $peakMB = $mb }
            if (($pr.WorkingSet64 / 1MB) -gt $MaxMB) {
                $reason = "KILLED: working set above $MaxMB MB"
                Stop-Process -Id $childPid -Force -ErrorAction SilentlyContinue
                Stop-Process -Id $launcherPid -Force -ErrorAction SilentlyContinue
                break
            }
        }
    }
    if ($sw.Elapsed.TotalSeconds -gt $TimeoutSec) {
        $reason = "KILLED: timeout $TimeoutSec s"
        if ($childPid) { Stop-Process -Id $childPid -Force -ErrorAction SilentlyContinue }
        Stop-Process -Id $launcherPid -Force -ErrorAction SilentlyContinue
        break
    }
    Start-Sleep -Milliseconds 300
}
"GUARD_RESULT launcher=$launcherPid child=$childPid peak_MB=$([int]$peakMB) wall_s=$([int]$sw.Elapsed.TotalSeconds) status=$reason" | Out-File -Append -Encoding ascii $Log
Write-Output "GUARD_RESULT peak_MB=$([int]$peakMB) wall_s=$([int]$sw.Elapsed.TotalSeconds) status=$reason"
