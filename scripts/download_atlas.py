"""Explicit opt-in download; no large data are downloaded by the default workflow."""
import argparse
import csv
import os
import subprocess
from pathlib import Path
from urllib.parse import urlparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default="config/downloads.atlas.tsv")
    args = parser.parse_args()
    with open(args.manifest) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        url = urlparse(row["url"])
        target = Path(row["path"])
        if url.scheme != "https" or url.hostname != "cf.10xgenomics.com":
            raise ValueError("This downloader accepts only the documented 10x HTTPS host")
        if target.is_absolute() or ".." in target.parts or target.parts[0] != "data":
            raise ValueError("Download destinations must stay under data/")
        if target.exists():
            print(f"Exists, left unchanged: {target}")
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + ".partial")
        subprocess.run(["curl", "--fail", "--location", "--proto", "=https", "--proto-redir", "=https",
                        "--retry", "3", "--continue-at", "-", "--output", str(partial), row["url"]], check=True)
        if partial.stat().st_size == 0:
            raise ValueError(f"Empty download: {partial}")
        os.replace(partial, target)
        print(f"Downloaded {target}; inspect and pin actual checksums before analysis")


if __name__ == "__main__":
    main()
