#!/usr/bin/env python3

import json
import re
from pathlib import Path

panel = json.loads(Path(".codeocean/app-panel.json").read_text())
assert panel["named_parameters"] is True
assert "help_text" not in Path(".codeocean/app-panel.json").read_text()

parameters = panel["parameters"]
names = [parameter.get("param_name") for parameter in parameters]
assert all(names), "Every App Panel parameter must define param_name"
assert len(names) == len(set(names)), "App Panel param_name values must be unique"

main_source = Path("code/main.R").read_text()
cli_names = set(re.findall(r'make_option\("--([a-z0-9_]+)"', main_source))
assert set(names).issubset(cli_names), sorted(set(names) - cli_names)

by_name = {parameter["param_name"]: parameter for parameter in parameters}
expected_defaults = {
    "input_type": "auto",
    "analysis_mode": "raw_counts",
    "batch_effect_columns": "auto",
    "filter_low_expression": "auto",
    "return_batch_corrected_values": "auto",
    "remove_donor_effect_for_downstream": "auto",
    "normalization_method": "auto",
    "write_normalization_diagnostics": "auto",
    "summarization_method": "auto",
}
for name, expected in expected_defaults.items():
    assert by_name[name].get("default_value") == expected, (name, by_name[name].get("default_value"))

assert by_name["analysis_mode"]["extra_data"] == [
    "raw_counts",
    "harmony_mean_expression",
    "sct_mean_expression",
]
assert by_name["normalization_method"]["extra_data"][:2] == ["auto", "None"]
assert by_name["pseudobulk_manifest"]["category"] == "advanced"

print("OMIX DEG App Panel bindings and mode-aware defaults passed")
