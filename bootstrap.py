"""Create the project folders and obtain available upstream inputs.

Python 3.10+, standard library only. Run from any working directory:
    python bootstrap.py --folders-only
    python bootstrap.py
    python bootstrap.py --manifest additional_downloads.json

Keep bootstrap_directories.json and bootstrap_manual_inputs.json beside this file.
An optional manifest is a JSON list of objects with url, path and optional sha256
and format fields. Supported formats: file, gzip, zip-gmt, tar-gzip.
Paths are relative to this script. Archive formats extract to the path's parent.

This prepares inputs; it does not run analyses or install Python/R packages.
Default reference URLs use current releases, not the publication's frozen versions.
Use archived URLs and SHA256 hashes in a manifest for exact reproducibility.
Existing nonempty files are preserved, and checksum mismatches are reported.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent
DOWNLOADS = [
    {"url": "https://purl.obolibrary.org/obo/go/go-basic.json", "path": "resources/go/go-basic.json"},
    {"url": "https://current.geneontology.org/annotations/rgd.gaf.gz", "path": "resources/go/rgd.gaf.gz"},
    {"url": "https://current.geneontology.org/annotations/mgi.gaf.gz", "path": "resources/go/mgi.gaf.gz"},
    {"url": "https://reactome.org/download/current/ReactomePathways.gmt.zip", "path": "resources/reactome/reactome_pathways.gmt", "format": "zip-gmt"},
    {"url": "https://reactome.org/download/current/ReactomePathways.txt", "path": "resources/reactome/reactome_pathways.txt"},
    {"url": "https://reactome.org/download/current/ReactomePathwaysRelation.txt", "path": "resources/reactome/reactome_pathways_relation.txt"},
    {"url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281833/suppl/GSE281833_count_matrix.tsv.gz", "path": "resources/interferon/gse28183_tec/source/GSE281833_count_matrix.tsv", "format": "gzip"},
    {"url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161664/suppl/GSE161664_RAW.tar", "path": "resources/interferon/gse161664_lung/source/GSE161664_RAW.tar", "format": "tar-gzip"},
]


def local_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    if not path.is_relative_to(ROOT) or path == ROOT:
        raise ValueError(f"Path must be inside the project: {value}")
    return path


def digest(path: Path) -> str:
    sha = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
    return sha.hexdigest()


def present(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def install(source, destination: Path) -> None:
    """Copy through a temporary sibling so interruptions leave no partial target."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=destination.parent, delete=False) as output:
        temporary = Path(output.name)
        try:
            shutil.copyfileobj(source, output)
        except BaseException:
            output.close()
            temporary.unlink(missing_ok=True)
            raise
    try:
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)


def unpack_tar(archive: Path, destination: Path) -> None:
    # Read only flat regular files; never extract links or trust archive paths.
    with tarfile.open(archive) as bundle:
        for member in bundle:
            if not member.isfile():
                continue
            name = Path(member.name)
            if name.name != member.name or "\\" in member.name:
                raise ValueError(f"Expected flat GEO archive; unexpected member: {member.name}")
            target = destination / name.name
            compressed = target.suffix == ".gz"
            if compressed:
                target = target.with_suffix("")
            if present(target):
                continue
            stream = bundle.extractfile(member)
            if stream is None:
                raise ValueError(f"Cannot read archive member: {member.name}")
            with stream:
                if compressed:
                    with gzip.GzipFile(fileobj=stream) as decoded:
                        install(decoded, target)
                else:
                    install(stream, target)


def obtain(item: dict) -> None:
    target = local_path(item["path"])
    kind = item.get("format", "file")
    if kind not in {"file", "gzip", "zip-gmt", "tar-gzip"}:
        raise ValueError(f"Unsupported format: {kind}")
    checksum = item.get("sha256")  # Hash of downloaded bytes, before extraction.
    if kind == "file" and present(target):
        if checksum and digest(target).lower() != checksum.lower():
            raise ValueError(f"Existing file has incorrect SHA256: {target}")
        print(f"Present: {item['path']}")
        return
    if kind in {"gzip", "zip-gmt"} and present(target) and not checksum:
        print(f"Present: {item['path']} (unverified existing extracted file)")
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    # The tar is retained so a second run can repair missing extracted files.
    with tempfile.TemporaryDirectory(dir=ROOT) as scratch:
        archive = Path(scratch) / "download"
        if kind == "tar-gzip" and present(target):
            archive = target
        else:
            print(f"Downloading: {item['url']}", flush=True)
            request = urllib.request.Request(item["url"], headers={"User-Agent": "research-bootstrap/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response:
                install(response, archive)
        sha = digest(archive)
        if checksum and sha.lower() != checksum.lower():
            raise ValueError(f"SHA256 mismatch for {item['path']}")
        print(f"  Download SHA256: {sha}")
        if kind == "file":
            with archive.open("rb") as stream:
                install(stream, target)
        elif kind == "gzip":
            if not present(target):
                with gzip.open(archive, "rb") as stream:
                    install(stream, target)
        elif kind == "zip-gmt":
            with zipfile.ZipFile(archive) as bundle:
                members = [name for name in bundle.namelist() if name.lower().endswith(".gmt")]
                if len(members) != 1:
                    raise ValueError(f"Expected one GMT file, found {members}")
                if not present(target):
                    with bundle.open(members[0]) as stream:
                        install(stream, target)
        else:
            unpack_tar(archive, target.parent)
            if archive != target:
                with archive.open("rb") as stream:
                    install(stream, target)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folders-only", action="store_true", help="Create directories without downloading")
    parser.add_argument("--manifest", type=Path, help="Additional downloads; matching paths override default sources")
    args = parser.parse_args()
    directories = json.loads((ROOT / "bootstrap_directories.json").read_text(encoding="utf-8"))
    for value in directories:
        local_path(value).mkdir(parents=True, exist_ok=True)
    print(f"Prepared {len(directories)} project directories under {ROOT}")
    failures = []
    if not args.folders_only:
        downloads = {item["path"]: item for item in DOWNLOADS}
        if args.manifest:
            additional = json.loads(args.manifest.read_text(encoding="utf-8"))
            for item in additional:
                downloads[item["path"]] = item
        print("Default reference URLs contain current releases; they may differ from the manuscript inputs.")
        for item in downloads.values():
            try:
                obtain(item)
            except Exception as error:
                failures.append(item["path"])
                print(f"FAILED: {item['path']}: {error}", file=sys.stderr)
    manual = json.loads((ROOT / "bootstrap_manual_inputs.json").read_text(encoding="utf-8"))
    missing = [value for value in manual if not present(local_path(value))]
    if missing:
        print("\nInputs without a recorded automatic source (supply these separately):")
        for value in missing:
            print(f"  {value}")
    print("\nBootstrap does not regenerate standardized data, enrichment results, or manuscript inputs.")
    print("Raw rodent CEL/RNA-seq inputs also require separate acquisition; they are not bundled here.")
    print(f"Summary: {len(failures)} failed downloads; {len(missing)} missing manual inputs.")
    return 1 if failures or (missing and not args.folders_only) else 0


if __name__ == "__main__":
    raise SystemExit(main())
