"""Filter HMDB extract TSV by status values.

Default behavior keeps only rows with status in:
- quantified
- detected
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def filter_hmdb_extract_by_status(
    input_file: str,
    output_file: str,
    keep_statuses: tuple[str, ...] = ("quantified", "detected"),
) -> tuple[int, int]:
    """Filter HMDB TSV rows by status.

    Returns:
        (input_row_count, kept_row_count)
    """

    input_path = Path(input_file)
    output_path = Path(output_file)
    keep = {status.strip().lower() for status in keep_statuses if status.strip()}

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8", newline="") as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Input file has no header: {input_file}")

        status_column = "status"
        if status_column not in reader.fieldnames:
            raise ValueError(
                f"Input file must contain a '{status_column}' column. "
                f"Found columns: {reader.fieldnames}"
            )

        with output_path.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
            writer.writeheader()

            input_rows = 0
            kept_rows = 0
            for row in reader:
                input_rows += 1
                row_status = (row.get(status_column) or "").strip().lower()
                if row_status in keep:
                    writer.writerow(row)
                    kept_rows += 1

    return input_rows, kept_rows


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Filter HMDB extract TSV to selected status values.",
    )
    parser.add_argument(
        "input_file",
        nargs="?",
        default="data/hmdb_metabolites_extract.tsv",
        help="Input TSV file.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        default="data/hmdb_metabolites_extract_quantified_detected.tsv",
        help="Output TSV file.",
    )
    parser.add_argument(
        "--keep-status",
        nargs="+",
        default=["quantified", "detected"],
        help="Status values to keep (case-insensitive).",
    )
    return parser


def main() -> None:
    args = build_argument_parser().parse_args()
    input_rows, kept_rows = filter_hmdb_extract_by_status(
        input_file=args.input_file,
        output_file=args.output_file,
        keep_statuses=tuple(args.keep_status),
    )
    print(f"Read {input_rows} rows from {args.input_file}")
    print(f"Wrote {kept_rows} rows to {args.output_file}")


if __name__ == "__main__":
    main()
