#!/usr/bin/env python3
"""fetch_data.py -- download and unpack the Zenodo data archive.

Data are not stored in this git repository (see data/README.md for the expected
layout and full provenance notes); they are hosted on Zenodo:

    DOI: 10.5281/zenodo.21238983
    https://doi.org/10.5281/zenodo.21238983

(The 69 combinatorial ProteinGym assays bundled under data/proteingym/ are a
faithful subset of the official ProteinGym archive, 10.5281/zenodo.15293562.)

Usage:
    /usr/bin/python3 scripts/fetch_data.py
    make data

Safe to re-run: if data/cleaned/manifest.csv already exists, this is a no-op, so
it never clobbers a local development copy of the data. Uses only the standard
library (urllib, zipfile) -- no new runtime dependencies.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(ROOT, "data")
MANIFEST = os.path.join(DATA_DIR, "cleaned", "manifest.csv")

RECORD_ID = "21238983"
DOI = "10.5281/zenodo.21238983"
DOI_URL = f"https://doi.org/{DOI}"
API_URL = f"https://zenodo.org/api/records/{RECORD_ID}"

MANUAL_INSTRUCTIONS = f"""
Could not download the data automatically ({{reason}}).

Please fetch it manually:

  1. Open {DOI_URL} in a browser (Zenodo record {RECORD_ID}).
  2. Download the archive (e.g. compressible-data.zip) from the "Files" panel.
  3. Unzip it so that the result is:
       data/cleaned/manifest.csv
       data/cleaned/*.csv
       data/proteingym/proteingym_ref.csv
       data/proteingym/DMS_ProteinGym_substitutions/*.csv
  4. Re-run `make verify` / `make reproduce` once data/cleaned/manifest.csv exists.

See data/README.md for the full expected layout and provenance notes.
""".strip()


def already_present() -> bool:
    return os.path.isfile(MANIFEST)


def resolve_download_url():
    """Query the Zenodo REST API for the record's file download link.

    Returns (url, filename).
    """
    req = urllib.request.Request(API_URL, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        record = json.loads(resp.read().decode("utf-8"))
    files = record.get("files") or []
    if not files:
        raise RuntimeError("Zenodo record has no files listed")
    # Prefer a .zip archive; fall back to the first file.
    zip_files = [f for f in files if f.get("key", "").endswith(".zip")]
    chosen = zip_files[0] if zip_files else files[0]
    links = chosen.get("links", {}) or {}
    url = links.get("self") or links.get("download")
    if not url:
        raise RuntimeError("Could not find a download link in the Zenodo file record")
    return url, chosen.get("key", "data.zip")


def download(url: str, dest_path: str) -> None:
    req = urllib.request.Request(url, headers={"Accept": "application/octet-stream"})
    with urllib.request.urlopen(req, timeout=60) as resp, open(dest_path, "wb") as out:
        total = int(resp.headers.get("Content-Length", 0) or 0)
        chunk = 1 << 20
        downloaded = 0
        while True:
            block = resp.read(chunk)
            if not block:
                break
            out.write(block)
            downloaded += len(block)
            if total:
                pct = 100 * downloaded / total
                print(f"\r  downloading... {downloaded / 1e6:8.1f} / {total / 1e6:.1f} MB "
                      f"({pct:5.1f}%)", end="", flush=True)
            else:
                print(f"\r  downloading... {downloaded / 1e6:8.1f} MB", end="", flush=True)
    print()


def extract(zip_path: str) -> None:
    """Unzip into a scratch dir, locate cleaned/ and proteingym/ wherever they
    landed (the archive may or may not wrap them in a top-level directory), and
    copy them into data/.
    """
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp)

        found_cleaned = None
        found_proteingym = None
        for root, _dirs, files in os.walk(tmp):
            if os.path.basename(root) == "cleaned" and "manifest.csv" in files:
                found_cleaned = root
            if os.path.basename(root) == "proteingym":
                found_proteingym = root

        if found_cleaned is None:
            raise RuntimeError("Downloaded archive does not contain a cleaned/manifest.csv")

        os.makedirs(DATA_DIR, exist_ok=True)

        dest_cleaned = os.path.join(DATA_DIR, "cleaned")
        if os.path.exists(dest_cleaned):
            shutil.rmtree(dest_cleaned)
        shutil.copytree(found_cleaned, dest_cleaned)

        if found_proteingym is not None:
            dest_pg = os.path.join(DATA_DIR, "proteingym")
            if os.path.exists(dest_pg):
                shutil.rmtree(dest_pg)
            shutil.copytree(found_proteingym, dest_pg)
        else:
            print("  [warn] no 'proteingym' directory found in archive; ProteinGym analyses "
                  "(scripts/run_proteingym.py) will not be reproducible until it is added "
                  "manually -- see data/README.md.")


def main() -> int:
    if already_present():
        print("data already present (data/cleaned/manifest.csv exists) -- no-op.")
        return 0

    print(f"Fetching data archive from Zenodo record {RECORD_ID} ({DOI_URL}) ...")
    try:
        url, filename = resolve_download_url()
    except (urllib.error.URLError, urllib.error.HTTPError, OSError, RuntimeError, ValueError) as e:
        print(MANUAL_INSTRUCTIONS.format(reason=f"could not resolve download URL via the "
                                                  f"Zenodo API: {e}"))
        return 1

    with tempfile.TemporaryDirectory() as tmp:
        zip_path = os.path.join(tmp, filename)
        try:
            print(f"  file: {filename}")
            download(url, zip_path)
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            print(MANUAL_INSTRUCTIONS.format(reason=f"download failed: {e}"))
            return 1

        try:
            print("Extracting into data/ ...")
            extract(zip_path)
        except (zipfile.BadZipFile, RuntimeError, OSError) as e:
            print(MANUAL_INSTRUCTIONS.format(reason=f"extraction failed: {e}"))
            return 1

    if not already_present():
        print(MANUAL_INSTRUCTIONS.format(reason="data/cleaned/manifest.csv missing after extraction"))
        return 1

    print(f"OK: data/cleaned/manifest.csv present. Data fetched from {DOI_URL}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
