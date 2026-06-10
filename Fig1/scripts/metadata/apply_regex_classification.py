#!/usr/bin/env python3
"""Apply regex-based classification rules to library protocol and spot descriptor fields."""

import re
import pandas as pd
from typing import Dict, List, Optional, Tuple

class TechClassifier:
    """Regex-based classifier for sequencing technology and library layout."""

    def __init__(self):
        self.patterns = self._init_patterns()

    def _init_patterns(self) -> Dict[str, Dict[str, str]]:
        """Initialize regex patterns for classification."""
        return {
            # Single-cell technologies
            'sc_10x': {
                'pattern': r'(?i)(10x|10X|chromium).{0,80}(gem|droplet|v[23]|3[\'′].*gene|single\s*cell\s*3)',
                'confidence': 'high',
                'category': 'single_cell',
            },
            'sc_celseq': {
                'pattern': r'(?i)cel-?seq',
                'confidence': 'high',
                'category': 'single_cell',
            },
            'sc_mars': {
                'pattern': r'(?i)mars-?seq',
                'confidence': 'high',
                'category': 'single_cell',
            },
            'sc_smartseq': {
                'pattern': r'(?i)smart-?seq.{0,50}(ultra|v[234]|low\s*input)',
                'confidence': 'high',
                'category': 'single_cell',
            },
            'sc_generic': {
                'pattern': r'(?i)(drop-?seq|sci-?rna|geminim|split-?pool)',
                'confidence': 'high',
                'category': 'single_cell',
            },

            # Bulk/tissue technologies
            'bulk_trueseq': {
                'pattern': r'(?i)true-?seq.{0,50}(stranded|rna|total)',
                'confidence': 'high',
                'category': 'bulk',
            },
            'bulk_lexogen': {
                'pattern': r'(?i)lexogen|quant-?seq',
                'confidence': 'high',
                'category': 'bulk',
            },
            'bulk_nebnext': {
                'pattern': r'(?i)neb-?next.{0,50}(ultra|rna|directional)',
                'confidence': 'high',
                'category': 'bulk',
            },
            'bulk_bgi': {
                'pattern': r'(?i)(bgiseq|dnbseq)',
                'confidence': 'medium',
                'category': 'bulk',
            },

            # Bias indicators
            'bias_3prime': {
                'pattern': r'(?i)(3[\'′]|three\s*prime|oligo.?dt|quantseq.*3[\'′])',
                'confidence': 'high',
                'category': 'bias',
            },
            'bias_5prime': {
                'pattern': r'(?i)(5[\'′]|five\s*prime)',
                'confidence': 'high',
                'category': 'bias',
            },
            'bias_fullength': {
                'pattern': r'(?i)(full.?length|whole.?length|smart-?seq.*(?!3[\'′]))',
                'confidence': 'medium',
                'category': 'bias',
            },

            # Read configuration
            'layout_paired': {
                'pattern': r'(?i)(paired?.?end|pair.?end|2x\d+|pe\b)',
                'confidence': 'high',
                'category': 'layout',
            },
            'layout_single': {
                'pattern': r'(?i)(single.?end|se\b|1x\d+)',
                'confidence': 'high',
                'category': 'layout',
            },

            # Selection methods
            'selection_polya': {
                'pattern': r'(?i)(poly\s*a|polya|oligo.{0,5}dt|mrna.*select|select.*mRNA)',
                'confidence': 'high',
                'category': 'selection',
            },
            'selection_random': {
                'pattern': r'(?i)(random.*prim|random.*hexamer|whole.*transcr)',
                'confidence': 'high',
                'category': 'selection',
            },
            'selection_rrna_deplete': {
                'pattern': r'(?i)(ribo.?zero|rrna.*depl|deplete.*ribo)',
                'confidence': 'high',
                'category': 'selection',
            },
            'selection_small_rna': {
                'pattern': r'(?i)(small\s*rna|mirna|trna)',
                'confidence': 'high',
                'category': 'selection',
            },
        }

    def classify_protocol(self, text: str) -> Dict[str, List[Tuple[str, str, str]]]:
        """
        Classify library construction protocol text.

        Returns:
            Dict with categories as keys, list of (rule_name, match, confidence) as values
        """
        if not text or pd.isna(text):
            return {}

        text = str(text)
        results = {}

        for rule_name, rule_info in self.patterns.items():
            pattern = rule_info['pattern']
            category = rule_info['category']
            confidence = rule_info['confidence']

            if re.search(pattern, text):
                if category not in results:
                    results[category] = []
                results[category].append((rule_name, pattern[:50] + '...', confidence))

        return results

    def classify_spot_descriptor(self, xml_str: str) -> Dict[str, any]:
        """
        Parse spot descriptor XML for read structure.

        Returns:
            Dict with extracted features
        """
        if not xml_str or pd.isna(xml_str):
            return {}

        xml_str = str(xml_str)
        features = {}

        # Extract SPOT_LENGTH
        spot_len_match = re.search(r'<SPOT_LENGTH>(\d+)</SPOT_LENGTH>', xml_str)
        if spot_len_match:
            features['spot_length'] = int(spot_len_match.group(1))

        # Count READ_SPEC elements (indicates paired-end or multiplex)
        read_specs = re.findall(r'<READ_SPEC>', xml_str)
        features['num_reads'] = len(read_specs)

        # Check for technical reads (barcode/UMI indicators)
        has_technical = bool(re.search(r'<READ_CLASS>Technical\s*Read</READ_CLASS>', xml_str))
        features['has_technical_read'] = has_technical

        # Extract read types
        read_types = re.findall(r'<READ_TYPE>(\w+)</READ_TYPE>', xml_str)
        features['read_types'] = read_types

        # Infer layout from structure
        if len(read_specs) >= 2:
            features['inferred_layout'] = 'paired-end'
        else:
            features['inferred_layout'] = 'single-end'

        # Detect barcode signature: technical read + application read
        if has_technical and len(read_specs) >= 2:
            features['likely_barcode'] = True

        return features


def example_usage(protocol_file: str, spot_descriptor_file: str):
    """Example of loading and classifying real data."""
    classifier = TechClassifier()

    print("=" * 80)
    print("EXAMPLE: Classifying first 5 unique protocols")
    print("=" * 80)

    # Read first few unique protocols
    try:
        with open(protocol_file, 'r') as f:
            protocols = [line.strip() for line in f.readlines()[:6]]  # Skip header + 5 samples

        for i, protocol in enumerate(protocols[1:4], 1):  # Show 3 examples
            print(f"\n[Protocol {i}]")
            print(f"Text: {protocol[:150]}...")

            classification = classifier.classify_protocol(protocol)
            if classification:
                for category, matches in classification.items():
                    print(f"  {category.upper()}:")
                    for rule, pattern, confidence in matches:
                        print(f"    - {rule} ({confidence})")
            else:
                print("  [No matches]")

    except FileNotFoundError:
        print(f"Protocol file not found: {protocol_file}")

    print("\n" + "=" * 80)
    print("EXAMPLE: Parsing spot descriptors")
    print("=" * 80)

    try:
        with open(spot_descriptor_file, 'r') as f:
            descriptors = [line.strip() for line in f.readlines()[:10]]

        for i, desc in enumerate(descriptors[1:4], 1):
            print(f"\n[Descriptor {i}]")
            if desc.startswith('<'):
                features = classifier.classify_spot_descriptor(desc)
                print(f"  SPOT_LENGTH: {features.get('spot_length', 'N/A')}")
                print(f"  Layout: {features.get('inferred_layout', 'unknown')}")
                print(f"  Read types: {features.get('read_types', [])}")
                print(f"  Technical read: {features.get('has_technical_read', False)}")
                print(f"  Likely barcode: {features.get('likely_barcode', False)}")
            else:
                print(f"  Text: {desc}")
    except FileNotFoundError:
        print(f"Spot descriptor file not found: {spot_descriptor_file}")


def apply_to_dataframe(df: pd.DataFrame, protocol_col: str = 'experiment.library_construction_protocol') -> pd.DataFrame:
    """Apply classifier to a full dataframe."""
    classifier = TechClassifier()

    df['tech_classification'] = df[protocol_col].apply(
        lambda x: classifier.classify_protocol(x)
    )

    # Extract top category
    df['primary_category'] = df['tech_classification'].apply(
        lambda x: list(x.keys())[0] if x else 'unknown'
    )

    return df


if __name__ == '__main__':
    import sys

    protocol_file = 'unique_values/experiment.library_construction_protocol.txt'
    spot_descriptor_file = 'unique_values/experiment.spot_descriptor.txt'

    example_usage(protocol_file, spot_descriptor_file)
