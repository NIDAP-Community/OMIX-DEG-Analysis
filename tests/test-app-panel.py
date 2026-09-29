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

expected_order = [
    "counts",
    "metadata",
    "pseudobulk_manifest",
    "moo",
    "input_type",
    "analysis_mode",
    "gene_names_column",
    "sample_names_column",
    "samples_to_include",
    "contrast_variable_columns",
    "contrasts",
    "covariate_columns",
    "batch_effect_columns",
    "donor_variable_column",
    "filter_low_expression",
    "normalization_method",
    "return_batch_corrected_values",
    "remove_donor_effect_for_downstream",
    "summarization_method",
    "write_normalization_diagnostics",
]
assert names == expected_order, (names, expected_order)

main_source = Path("code/main.R").read_text()
cli_names = set(re.findall(r'make_option\("--([a-z0-9_]+)"', main_source))
assert set(names).issubset(cli_names), sorted(set(names) - cli_names)
assert cli_names == set(expected_order) | {"output_dir"}, sorted(cli_names)

cli_defaults = dict(
    re.findall(
        r'make_option\("--([a-z0-9_]+)"[^\n]*?default\s*=\s*"([^"]*)"',
        main_source,
    )
)
assert cli_defaults["input_type"] == "table"
assert cli_defaults["analysis_mode"] == "raw_counts"
assert cli_defaults["counts"] == ""
assert cli_defaults["metadata"] == ""
assert cli_defaults["pseudobulk_manifest"] == ""
assert cli_defaults["moo"] == ""

by_name = {parameter["param_name"]: parameter for parameter in parameters}
expected_defaults = {
    "input_type": "table",
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
assert by_name["input_type"]["extra_data"] == ["table", "moo"]
assert by_name["normalization_method"]["extra_data"][:2] == ["auto", "None"]
assert by_name["pseudobulk_manifest"]["category"] == "advanced"

file_controls = {"counts", "metadata", "pseudobulk_manifest", "moo"}
assert all(by_name[name]["type"] == "file" for name in file_controls)
auto_boolean_controls = {
    "filter_low_expression",
    "return_batch_corrected_values",
    "remove_donor_effect_for_downstream",
    "write_normalization_diagnostics",
}
for name in auto_boolean_controls:
    assert by_name[name]["extra_data"] == ["auto", "true", "false"]
assert "counts_file" not in cli_names
assert "metadata_file" not in cli_names
assert "moo_file" not in cli_names

print("OMIX DEG App Panel bindings and mode-aware defaults passed")
