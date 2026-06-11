#!/usr/bin/env python
"""
compare_token_dicts.py

Compare the token dictionaries of two Hugging Face datasets saved with
`dataset.save_to_disk`.  The JSON-Lines report contains:

    unique_count_dataset1  – number of unique tokens in ds1
    unique_count_dataset2  – number of unique tokens in ds2
    max_token_dataset1     – highest token-id seen in ds1
    max_token_dataset2     – highest token-id seen in ds2
    only_in_dataset1       – list of tokens only in ds1
    only_in_dataset2       – list of tokens only in ds2
    rogue_examples         – row locations of ds2-only tokens

Example
-------
python compare_token_dicts.py /path/ds1 /path/ds2 \
    --column input_ids --text-column text \
    --max-tokens 50 --max-examples 3 --output token_diff.jsonl
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Set, Iterable, Union

from datasets import Dataset, DatasetDict, load_from_disk


##############################################################################
# Helpers
##############################################################################

def iter_examples(ds: Union[Dataset, DatasetDict]) -> Iterable[tuple[str, int, dict]]:
    """Yield (split_name, row_idx, example) over all rows, flattening splits."""
    if isinstance(ds, DatasetDict):
        for split_name, split in ds.items():
            for i, ex in enumerate(split):
                yield split_name, i, ex
    else:
        for i, ex in enumerate(ds):
            yield "data", i, ex


def collect_tokens(ds: Union[Dataset, DatasetDict], column: str) -> Set[int]:
    """Return the union of token IDs from `column` across all splits."""
    if isinstance(ds, DatasetDict):
        return set().union(*(collect_tokens(s, column) for s in ds.values()))

    if column not in ds.column_names:
        raise ValueError(
            f"Column '{column}' not found. Available columns: {ds.column_names}"
        )

    vocab: Set[int] = set()
    for ex in ds:
        vocab.update(int(t) for t in ex[column])
    return vocab


def collect_rogue_examples(
    ds: Union[Dataset, DatasetDict],
    column: str,
    rogue_tokens: Set[int],
    max_tokens: int,
    max_examples: int,
    text_column: str | None,
) -> Dict[int, List[dict]]:
    """
    Map each rogue token → up to `max_examples` row references.
    Limits to `max_tokens` rogue tokens overall.
    """
    sample_map: Dict[int, List[dict]] = {
        tok: [] for tok in list(rogue_tokens)[:max_tokens]
    }
    active = set(sample_map)

    for split, idx, ex in iter_examples(ds):
        row_tokens = set(int(t) for t in ex[column]) & active
        if not row_tokens:
            continue

        for tok in row_tokens:
            lst = sample_map[tok]
            if len(lst) < max_examples:
                entry = {"split": split, "row_idx": idx}
                if text_column and text_column in ex:
                    entry["text"] = ex[text_column]
                lst.append(entry)

        if all(len(v) >= max_examples for v in sample_map.values()):
            break

    return sample_map


##############################################################################
# Main
##############################################################################

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare token dictionaries and locate rogue tokens."
    )
    parser.add_argument("ds1", type=Path, help="Path to first dataset directory")
    parser.add_argument("ds2", type=Path, help="Path to second dataset directory")
    parser.add_argument(
        "--column", default="input_ids",
        help="Column containing token sequences (default: 'input_ids')",
    )
    parser.add_argument(
        "--text-column", default=None,
        help="Optional raw-text column to include in rogue examples",
    )
    parser.add_argument(
        "--max-tokens", type=int, default=50,
        help="Store examples for at most this many rogue tokens (default 50)",
    )
    parser.add_argument(
        "--max-examples", type=int, default=3,
        help="Max example rows saved per rogue token (default 3)",
    )
    parser.add_argument(
        "--output", type=Path, default=Path("token_diff.jsonl"),
        help="Output JSON-Lines file (default token_diff.jsonl)",
    )
    args = parser.parse_args()

    print("🔄  Loading datasets …")
    ds1 = load_from_disk(args.ds1)
    ds2 = load_from_disk(args.ds2)

    print(f"🔍  Collecting tokens from column '{args.column}' …")
    vocab1 = collect_tokens(ds1, args.column)
    vocab2 = collect_tokens(ds2, args.column)

    n1, n2 = len(vocab1), len(vocab2)
    max1, max2 = (max(vocab1) if vocab1 else None,
                  max(vocab2) if vocab2 else None)
    only_in_1 = vocab1 - vocab2
    only_in_2 = vocab2 - vocab1

    # Console summary
    print(f"📊 Unique tokens — ds1: {n1:,} (max id {max1}) | "
          f"ds2: {n2:,} (max id {max2})")

    if not only_in_1 and not only_in_2:
        print("✅ Token dictionaries are identical. Nothing to write.")
        return

    print(f"❌ Dictionaries differ "
          f"({len(only_in_1):,} only in ds1, {len(only_in_2):,} only in ds2).")
    print("🔎  Locating rogue tokens in ds2 …")
    rogue_examples = collect_rogue_examples(
        ds2, column=args.column, rogue_tokens=only_in_2,
        max_tokens=args.max_tokens, max_examples=args.max_examples,
        text_column=args.text_column,
    )

    print(f"📝 Writing full report to {args.output} …")
    with args.output.open("w", encoding="utf-8") as f:
        json.dump(
            {
                "unique_count_dataset1": n1,
                "unique_count_dataset2": n2,
                "max_token_dataset1": max1,
                "max_token_dataset2": max2,
                "only_in_dataset1": sorted(only_in_1),
                "only_in_dataset2": sorted(only_in_2),
                "rogue_examples": rogue_examples,
            },
            f,
            ensure_ascii=False,
        )
        f.write("\n")

    print("✅ Done.")


if __name__ == "__main__":
    main()
