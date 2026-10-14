#!/usr/bin/env python3

"""Validate the managed DEG scientific export and its recorded provenance."""

import hashlib
import re
from pathlib import Path


source_record = Path("OMIX_MODULE_SOURCE.md").read_text()
expected_commit = "11bda7fb2cbfd94c0641f864d1c3a49a678f5af1"
expected_hash = "44cfe51ea649ed72aa5c7c55ddfa6aed907517cad408c52a34a96ca2c692a7f2"
managed_path = Path("code/functions/OMIX_DEG_Analysis.R")

commit_match = re.search(
    r"Canonical source reference:\*\*\s*\[`([0-9a-f]{40})`\]", source_record
)
assert commit_match, "Canonical source reference is missing"
assert commit_match.group(1) == expected_commit
assert "**Canonical module version:** `0.4.2`" in source_record
assert "**Canonical interface version:** `2`" in source_record

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
