#!/usr/bin/env python3

"""Validate the managed DEG scientific export and its recorded provenance."""

import hashlib
import re
from pathlib import Path


source_record = Path("OMIX_MODULE_SOURCE.md").read_text()
expected_commit = "018ef0e10bf29c69166b38ee6d81d4d71d18aeca"
expected_hash = "e0e64e4de04b91419c77d1d44ac607f62fa172bc1b99fd0936ebb3e6ae8eb13e"
managed_path = Path("code/functions/OMIX_DEG_Analysis.R")

commit_match = re.search(
    r"Canonical source reference:\*\*\s*\[`([0-9a-f]{40})`\]", source_record
)
assert commit_match, "Canonical source reference is missing"
assert commit_match.group(1) == expected_commit

managed_files = sorted(
    path.relative_to("code/functions").as_posix()
    for path in Path("code/functions").rglob("*")
    if path.is_file()
)
assert managed_files == ["OMIX_DEG_Analysis.R"], managed_files

observed_hash = hashlib.sha256(managed_path.read_bytes()).hexdigest()
assert observed_hash == expected_hash, observed_hash
assert source_record.count(f"`{expected_hash}`") >= 2
assert "`R/OMIX_DEG_Analysis.R`" in source_record
assert "`code/functions/OMIX_DEG_Analysis.R`" in source_record

print("OMIX DEG managed-source provenance passed")
