#!/usr/bin/env python3
"""Convert a results CSV (columns: name, time, y) to the column-oriented JSON
format used by result_chamber_sphere.json.

Usage:
    python csv_to_json.py input.csv [output.json]

If no output path is given, the input's extension is replaced with .json.
"""

import argparse
from pathlib import Path

import pandas as pd


def csv_to_json(input_path: Path, output_path: Path) -> None:
    df = pd.read_csv(input_path)
    # Column-oriented dict-of-dicts, matching pandas' default to_json layout
    # (e.g. {"name": {"0": ...}, "time": {"0": ...}, "y": {"0": ...}}).
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_json(output_path, orient="columns", indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Path to the input .csv file")
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        help="Path to the output .json file (defaults to input with .json suffix)",
    )
    args = parser.parse_args()

    output_path = args.output or args.input.with_suffix(".json")
    csv_to_json(args.input, output_path)
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
