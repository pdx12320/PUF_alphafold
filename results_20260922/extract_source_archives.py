#!/usr/bin/env python3
"""Recover exact original analyses from locally retained uploads; no network or model execution."""
from pathlib import Path
import argparse, hashlib, json, zipfile

PREFIXES = {
    "PUF_classification_v4_dynamic_fullCP_20260921_complete.zip":
        "PUF_alphafold/retrain_classification_v4_dynamic_fullCP_20260921/",
    "designated_five_construct_held_out_20260922.zip":
        "designated_five_construct_held_out_20260922/",
    "retrospective_ten_construct_holdout_20260922_results.zip":
        "retrospective_ten_construct_holdout_20260922/",
    "PUF_Dry_Lab_Wiki_and_Retraining_20260922.zip": "PUF_Wiki_20260922/",
}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archives", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--all-members", action="store_true",
                        help="Also recover binaries and large arrays; never load pickle objects.")
    args = parser.parse_args()
    manifest = json.loads((Path(__file__).parent / "source_archives.json").read_text())
    reference = {r["file"]: r for r in manifest}
    args.destination.mkdir(parents=True, exist_ok=True)
    log = []
    for name, prefix in PREFIXES.items():
        source = args.archives / name
        if not source.is_file():
            log.append({"file": name, "status": "missing"})
            continue
        if digest(source) != reference[name]["sha256"]:
            raise ValueError(f"SHA256 mismatch: {source}")
        target = (args.destination / Path(name).stem).resolve()
        count = 0
        with zipfile.ZipFile(source) as z:
            for member in z.infolist():
                if member.is_dir() or not member.filename.startswith(prefix):
                    continue
                relative = Path(member.filename[len(prefix):])
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"Unsafe archive path: {member.filename}")
                if not args.all_members and relative.suffix.lower() not in {
                    ".csv", ".json", ".md", ".txt", ".py", ".yaml", ".yml"
                }:
                    continue
                destination = (target / relative).resolve()
                if not destination.is_relative_to(target):
                    raise ValueError("Path traversal detected")
                destination.parent.mkdir(parents=True, exist_ok=True)
                if destination.exists():
                    if destination.read_bytes() != z.read(member):
                        raise FileExistsError(destination)
                    continue
                destination.write_bytes(z.read(member))
                count += 1
        log.append({"file": name, "status": "verified", "extracted": count})
    print(json.dumps(log, indent=2))

if __name__ == "__main__":
    main()
