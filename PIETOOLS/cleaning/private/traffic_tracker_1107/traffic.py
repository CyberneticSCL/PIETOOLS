#!/usr/bin/env python3
"""Append the GitHub traffic of a repository to CSV files that outlive the
14 days the GitHub traffic API retains.

Run once a day (by .github/workflows/traffic.yml). Each run fetches the last
14 days of views and clones, merges them by UTC day into the files below, and
records a dated snapshot of the quantities GitHub reports only as current
totals.

INPUT
  argv[1]            'str', directory of the record; created if absent.
  GITHUB_REPOSITORY  'str' environment variable, 'owner/name'.
  TRAFFIC_TOKEN      'str' environment variable, a fine-grained personal
                     access token with repository permission
                     'Administration: Read-only'. The workflow GITHUB_TOKEN
                     cannot read the traffic API.

OUTPUT, files in argv[1]
  views.csv      date,count,uniques: page views per UTC day.
  clones.csv     date,count,uniques: git clones per UTC day.
  referrers.csv  date,referrer,count,uniques: top 10 referring sites over
                 the 14 days ending on date.
  paths.csv      date,path,title,count,uniques: top 10 pages over the 14
                 days ending on date.
  downloads.csv  date,release,asset,downloads: cumulative downloads of each
                 release asset on date.
  repo.csv       date,stars,forks,watchers on date.
  README.md      monthly totals, latest release downloads, file index.
"""

import csv
import datetime
import json
import os
import sys
import urllib.request

API = "https://api.github.com"


def make_get(token):
    """Return get(path) -> decoded JSON of the GitHub REST API at path.
    INPUT  token 'str', bearer token.
    OUTPUT 'function'."""
    def get(path):
        req = urllib.request.Request(API + path, headers={
            "Authorization": "Bearer " + token,
            "Accept": "application/vnd.github+json",
            "X-GitHub-Api-Version": "2022-11-28",
            "User-Agent": "traffic-recorder"})
        with urllib.request.urlopen(req, timeout=60) as r:
            return json.load(r)
    return get


def read(path):
    """Rows of a CSV file as dicts; [] if the file does not exist.
    INPUT  path 'str'.  OUTPUT 'list' of 'dict'."""
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write(path, fields, rows):
    """Write rows (dicts keyed by fields) to a CSV file with '\\n' line ends.
    INPUT  path 'str'; fields 'list' of 'str'; rows 'list' of 'dict'."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def merge_daily(path, items):
    """Merge per-day counts into path. An API value replaces the stored row
    for its day: the newest day is partial until a later run.
    INPUT  path 'str'; items 'list' of 'dict' with fields timestamp
           (UTC midnight, ISO 8601), count, uniques, as the API returns."""
    rows = {r["date"]: r for r in read(path)}
    for it in items:
        day = it["timestamp"][:10]
        rows[day] = {"date": day, "count": it["count"],
                     "uniques": it["uniques"]}
    write(path, ["date", "count", "uniques"], [rows[d] for d in sorted(rows)])


def put_snapshot(path, fields, today, new):
    """Replace the rows dated today by new, so that a second run on one day
    leaves one snapshot.
    INPUT  path 'str'; fields 'list' of 'str' (first is 'date'); today 'str'
           YYYY-MM-DD; new 'list' of 'dict' keyed by fields."""
    rows = [r for r in read(path) if r["date"] != today] + new
    write(path, fields, rows)


def releases(get, repo):
    """All published releases of repo, following pagination.
    INPUT  get 'function' (make_get); repo 'str'.  OUTPUT 'list' of 'dict'."""
    out, page = [], 1
    while True:
        batch = get(f"/repos/{repo}/releases?per_page=100&page={page}")
        out += [r for r in batch if not r.get("draft")]
        if len(batch) < 100:
            return out
        page += 1


def summarize(out, repo, today):
    """Write README.md: totals per month and overall, and the latest
    snapshot of release downloads and repository counts.
    INPUT  out 'str', record directory; repo 'str'; today 'str' YYYY-MM-DD."""
    # Monthly sums of the daily series; uniques are summed per day.
    months = {}
    for name in ("views", "clones"):
        for r in read(os.path.join(out, name + ".csv")):
            m = months.setdefault(r["date"][:7], [0, 0, 0, 0])
            k = 0 if name == "views" else 2
            m[k] += int(r["count"])
            m[k + 1] += int(r["uniques"])
    days = [r["date"] for n in ("views", "clones")
            for r in read(os.path.join(out, n + ".csv"))]
    total = [sum(m[k] for m in months.values()) for k in range(4)]

    lines = [f"# Traffic record for {repo}", "",
             "Recorded daily by `.github/workflows/traffic.yml` on the default",
             "branch. GitHub retains traffic for 14 days; this branch keeps it.",
             "", f"Days recorded: {min(days) if days else '-'} to "
             f"{max(days) if days else '-'} (UTC). Updated {today}.", "",
             "| Month | Views | Unique views\\* | Clones | Unique clones\\* |",
             "|---|---:|---:|---:|---:|"]
    lines += [f"| {k} | {v[0]} | {v[1]} | {v[2]} | {v[3]} |"
              for k, v in sorted(months.items())]
    lines += [f"| **Total** | **{total[0]}** | **{total[1]}** | "
              f"**{total[2]}** | **{total[3]}** |", "",
              "\\* Sum of daily unique counts: a visitor active on k days "
              "counts k times.", ""]

    # Latest snapshot of release downloads (cumulative since upload).
    dl = read(os.path.join(out, "downloads.csv"))
    last = max((r["date"] for r in dl), default=None)
    lines += [f"## Release downloads, cumulative to {last or '-'}", "",
              "| Release | Asset | Downloads |", "|---|---|---:|"]
    lines += [f"| {r['release']} | {r['asset']} | {r['downloads']} |"
              for r in dl if r["date"] == last]
    lines += ["", "The GitHub API reports no count of source-archive "
              "(Download ZIP) downloads.", ""]

    # Latest repository counts.
    rp = read(os.path.join(out, "repo.csv"))
    if rp:
        r = rp[-1]
        lines += [f"Stars {r['stars']}, forks {r['forks']}, watchers "
                  f"{r['watchers']} on {r['date']}.", ""]

    lines += ["## Files", "",
              "- `views.csv`, `clones.csv`: count and uniques per UTC day.",
              "- `referrers.csv`, `paths.csv`: top 10 over the 14 days ending "
              "on the row date; overlapping windows, so do not sum rows.",
              "- `downloads.csv`: cumulative downloads per release asset.",
              "- `repo.csv`: stars, forks, watchers per day.", ""]
    with open(os.path.join(out, "README.md"), "w", encoding="utf-8",
              newline="\n") as f:
        f.write("\n".join(lines))


def main(out, repo, get, today):
    """Fetch the traffic of repo and merge it into the record in out.
    INPUT  out 'str'; repo 'str' owner/name; get 'function' (make_get);
           today 'str' YYYY-MM-DD, UTC date of the snapshot."""
    os.makedirs(out, exist_ok=True)
    t = f"/repos/{repo}/traffic"

    # Daily series: the API returns the last 14 UTC days.
    merge_daily(os.path.join(out, "views.csv"),
                get(t + "/views?per=day")["views"])
    merge_daily(os.path.join(out, "clones.csv"),
                get(t + "/clones?per=day")["clones"])

    # 14-day top-10 lists, kept as dated snapshots.
    put_snapshot(os.path.join(out, "referrers.csv"),
                 ["date", "referrer", "count", "uniques"], today,
                 [{"date": today, "referrer": r["referrer"],
                   "count": r["count"], "uniques": r["uniques"]}
                  for r in get(t + "/popular/referrers")])
    put_snapshot(os.path.join(out, "paths.csv"),
                 ["date", "path", "title", "count", "uniques"], today,
                 [{"date": today, "path": p["path"], "title": p["title"],
                   "count": p["count"], "uniques": p["uniques"]}
                  for p in get(t + "/popular/paths")])

    # Cumulative counts, kept as dated snapshots.
    put_snapshot(os.path.join(out, "downloads.csv"),
                 ["date", "release", "asset", "downloads"], today,
                 [{"date": today, "release": r["tag_name"],
                   "asset": a["name"], "downloads": a["download_count"]}
                  for r in releases(get, repo) for a in r["assets"]])
    info = get(f"/repos/{repo}")
    put_snapshot(os.path.join(out, "repo.csv"),
                 ["date", "stars", "forks", "watchers"], today,
                 [{"date": today, "stars": info["stargazers_count"],
                   "forks": info["forks_count"],
                   "watchers": info["subscribers_count"]}])

    summarize(out, repo, today)


if __name__ == "__main__":
    token = os.environ.get("TRAFFIC_TOKEN", "")
    if not token:
        sys.exit("TRAFFIC_TOKEN is empty: add a repository secret holding a "
                 "fine-grained token with 'Administration: Read-only'.")
    main(sys.argv[1], os.environ["GITHUB_REPOSITORY"], make_get(token),
         datetime.datetime.now(datetime.timezone.utc).date().isoformat())
