import os, requests
from collections import defaultdict

OWNER, REPO = "CyberneticSCL", "PIETOOLS"
token = os.getenv("GITHUB_TOKEN")  # optional

headers = {
    "Accept": "application/vnd.github+json",
    "X-GitHub-Api-Version": "2022-11-28",
    "User-Agent": "LANL-Sachin-PIETOOLS",  # REQUIRED
}
if token:
    headers["Authorization"] = f"Bearer {token}"

assets, page = [], 1
while True:
    url = f"https://api.github.com/repos/{OWNER}/{REPO}/releases?per_page=100&page={page}"
    r = requests.get(url, headers=headers, timeout=30)
    r.raise_for_status()
    batch = r.json()
    if not batch:
        break
    for rel in batch:
        for a in rel.get("assets") or []:
            assets.append({
                "tag": rel.get("tag_name"),
                "name": a.get("name"),
                "download_count": int(a.get("download_count") or 0),
            })
    page += 1

if not assets:
    print("No release assets found (check rate limit or headers).")
else:
    for row in sorted(assets, key=lambda x: (x["tag"] or "", x["name"] or "")):
        print(f"{row['tag']}: {row['name']}  -> {row['download_count']}")
    totals = defaultdict(int)
    for r in assets: totals[r["tag"]] += r["download_count"]
    print("\nPer-release totals:")
    for tag, total in sorted(totals.items(), key=lambda kv: kv[1], reverse=True):
        print(f"{tag}: {total}")
    print("\nGrand total:", sum(r["download_count"] for r in assets))
