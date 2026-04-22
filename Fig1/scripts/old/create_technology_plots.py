#!/usr/bin/env python3
"""
Script to create sequencing technology handling plots.
Joins accession list, processing outcomes, and manual annotations.
Creates multi-level pie charts showing processing outcomes by technology category.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import seaborn as sns
from collections import defaultdict
import argparse
from pathlib import Path

DRAW_OUTCOMES=True

def load_and_process_data():
    """Load and join the three input tables."""
    
    # Load accession list
    accessions_df = pd.read_csv('data/75k_unstable/ZF_SraEsearch-2025-06-22.csv')
    print(f"Loaded {len(accessions_df)} accessions")
    
    # Load seq-detective judgement data with correct header
    judgement_df = pd.read_csv('data/75k_unstable/seq-detective-judgement-summary-all.txt', 
                              sep='\t', header=None)
    # Set column names based on clarified structure: prefix, file1, file2, judgement1, judgement2, reason
    judgement_df.columns = ['prefix', 'file1', 'file2', 'judgement1', 'judgement2', 'reason']
    
    # Remove _subset suffix from prefix to get accession
    judgement_df['accession'] = judgement_df['prefix'].str.replace('_subset', '', regex=False)
    # Also try removing _subsample suffix in case that's the pattern
    judgement_df['accession'] = judgement_df['accession'].str.replace('_subsample', '', regex=False)
    
    # Create combined processing outcome (mate1+mate2 determination)
    judgement_df['outcome'] = judgement_df['judgement1'].astype(str) + judgement_df['judgement2'].astype(str)
    # Clean up NaN combinations
    judgement_df.loc[judgement_df['judgement2'].isna(), 'outcome'] = judgement_df['judgement1'].astype(str)
    
    print(f"Loaded {len(judgement_df)} judgement records")
    
    # Load manual annotations
    annotations_df = pd.read_csv('data/zf-core-v2-74K_problematic_with_bulk.kmers.allcols.csv')
    print(f"Loaded {len(annotations_df)} manual annotations")
    
    # Join accessions with judgement data
    merged_df = pd.merge(accessions_df, judgement_df[['accession', 'outcome', 'judgement1', 'judgement2']], 
                        left_on='Run', right_on='accession', how='left')
    
    # Join with manual annotations (include all hierarchy columns)
    annotation_cols = ['accession', 'resolution', 'group', 'technology', 'chemistry_ver', 'library_var']
    final_df = pd.merge(merged_df, annotations_df[annotation_cols], 
                       left_on='Run', right_on='accession', how='left', suffixes=('', '_manual'))
    
    print(f"Final merged dataset: {len(final_df)} records")
    
    return final_df

def categorize_technologies(df):
    """
    Categorize technologies using the hierarchy: resolution → group → technology → chemistry_ver → library_var
    """
    
    df = df.copy()
    
    # Flag entries without manual annotations at different levels
    df['has_resolution'] = ~df['resolution'].isna()
    df['has_group'] = ~df['group'].isna()
    df['has_technology'] = ~df['technology'].isna()
    df['has_chemistry_ver'] = ~df['chemistry_ver'].isna()
    df['has_library_var'] = ~df['library_var'].isna()
    
    # Create hierarchical categories with fallback propagation
    
    # Level 1: Resolution (highest level)
    df['tech_level1'] = df['resolution'].fillna('unknown')
    
    # Level 2: Group 
    df['tech_level2'] = df['group'].fillna(df['tech_level1'])
    
    # Level 3: Technology (medium resolution)
    df['tech_level3'] = df['technology'].fillna(df['tech_level2'])
    
    # Level 4: Chemistry version (high resolution) 
    df['tech_level4'] = df['chemistry_ver'].fillna(df['tech_level3'])
    
    # Level 5: Library variant (highest resolution)
    df['tech_level5'] = df['library_var'].fillna(df['tech_level4'])
    
    # Create simplified categories for visualization
    def simplify_resolution(res):
        if pd.isna(res) or res == 'unknown':
            return 'Unknown'
        return str(res).title()
    
    def simplify_group(group):
        if pd.isna(group) or group == 'unknown':
            return 'Other'
        # Group similar categories
        group_str = str(group)
        if 'FACS' in group_str:
            return 'FACS'
        elif 'Illumina' in group_str or 'Multiplexing' in group_str:
            return 'Illumina Multiplexing'
        elif 'Droplet' in group_str:
            return 'Droplet-based'
        elif 'Microwell' in group_str:
            return 'Microwell'
        elif 'RNA' in group_str or 'rRNA' in group_str:
            return 'RNA-focused'
        else:
            return 'Other'
    
    def simplify_technology(tech):
        if pd.isna(tech) or tech == 'unknown':
            return 'Other'
        tech_str = str(tech)
        if 'Smart-seq' in tech_str or 'SMART-seq' in tech_str:
            return 'Smart-seq'
        elif 'DeTCT' in tech_str:
            return 'DeTCT'
        elif '10x' in tech_str:
            return '10x Genomics'
        elif 'CEL-Seq' in tech_str:
            return 'CEL-Seq'
        elif 'CORALL' in tech_str:
            return 'CORALL'
        elif 'sci-RNA-seq' in tech_str:
            return 'sci-RNA-seq'
        elif 'inDrop' in tech_str:
            return 'inDrop'
        elif 'CLIP' in tech_str or 'RIP' in tech_str or 'ChIP' in tech_str:
            return 'Protein-RNA Interactions'
        else:
            return 'Other'
    
    # Apply simplifications
    df['resolution_simple'] = df['tech_level1'].apply(simplify_resolution)
    df['group_simple'] = df['tech_level2'].apply(simplify_group) 
    df['technology_simple'] = df['tech_level3'].apply(simplify_technology)
    
    # Flag entries that need propagation from lower resolution
    df['needs_technology_propagation'] = df['technology'].isna() & ~df['group'].isna()
    df['needs_group_propagation'] = df['group'].isna() & ~df['resolution'].isna()
    df['needs_any_annotation'] = df['technology'].isna() & df['group'].isna() & df['resolution'].isna()
    
    return df

def analyze_processing_outcomes(df):
    """Analyze processing outcomes by technology category."""
    
    # Filter out rows with missing outcomes
    df_clean = df[df['outcome'].notna()].copy()
    
    if len(df_clean) == 0:
        print("ERROR: No data left after filtering!")
        return pd.DataFrame(), df_clean
    
    # Count processing outcomes by technology at different resolution levels
    outcome_counts = df_clean.groupby(['resolution_simple', 'technology_simple', 'outcome']).size().reset_index(name='count')
    
    # Calculate totals by technology
    total_by_tech = df_clean.groupby(['resolution_simple', 'technology_simple']).size().reset_index(name='total')
    outcome_props = pd.merge(outcome_counts, total_by_tech, on=['resolution_simple', 'technology_simple'])
    outcome_props['proportion'] = outcome_props['count'] / outcome_props['total']
    
    return outcome_props, df_clean

def create_multilevel_pie_chart(outcome_props, df_clean, output_path='technology_outcomes_plot.svg'):
    """
    Create multi-level pie chart using polar coordinates showing processing outcomes by technology.
    """
    
    fig, ax = plt.subplots(figsize=(14, 12), subplot_kw=dict(projection='polar'))
    
    # Define colors for outcomes
    outcome_colors = {
        'BB': '#2E8B57',  # Both mates good - dark green
        'TB': '#FFD700',  # Mate1 bad, Mate2 good - gold
        'BT': '#FF8C00',  # Mate1 good, Mate2 bad - dark orange
        'TT': '#DC143C',  # Both mates bad - crimson
        'B': '#4169E1',   # Single end good - royal blue
        'T': '#8B0000',   # Single end bad - dark red
        'nan': '#808080'  # Missing - gray
    }
    
    # Get technology categories sorted by count
    tech_totals = df_clean.groupby('technology_simple').size().sort_values(ascending=False)
    tech_categories = tech_totals.index.tolist()
    
    # Define consistent color mapping for technology categories
    tech_color_map = {
        'Smart-seq': '#03185a',
        'DeTCT': '#fbcdfa', 
        '10x Genomics': '#818431',
        'CEL-Seq': '#226060',
        'CORALL': '#ee9f6f',
        'sci-RNA-seq': '#124364',
        'inDrop': '#fcb4b3',
        'Protein-RNA Interactions': '#4d734d',
        'Other': '#be9036'
    }
    
    # Import batlow colormap
    try:
        import matplotlib.colormaps as cm
        batlow = cm["batlowWS"]
        # Generate evenly spaced colors from batlow
        n_colors = len(tech_color_map)
        batlow_colors = [batlow(i / (n_colors - 1)) for i in range(n_colors)]
        
        # Update color map with batlow colors
        for i, tech_name in enumerate(tech_color_map.keys()):
            if tech_name not in tech_color_map:
                tech_color_map[tech_name] = batlow_colors[i]
    except:
        # Fallback to manual batlow-inspired colors if batlow not available
        batlow_colors = [
            '#03185a', '#0B2F8A', '#1E4B99', '#3568A6', '#4C84B1',
            '#63A0BC', '#7ABCC7', '#92D8D2', '#AAF4DD'
        ]
        for i, tech_name in enumerate(tech_color_map.keys()):
            if tech_name not in tech_color_map:
                tech_color_map[tech_name] = batlow_colors[i]
    
    # Calculate total samples
    total_all = len(df_clean)
    
    # Prepare data for polar bar chart
    gap_fraction = 0.02  # Fraction of total circle for gaps
    n_categories = len([cat for cat in tech_categories if tech_totals[cat]/total_all > 0.01])
    total_gap = gap_fraction * 2 * np.pi  # Total gap space
    gap_width = total_gap / max(n_categories, 1) if n_categories > 0 else 0
    available_space = 2 * np.pi - total_gap  # Space left for data
    
    theta_start = np.pi + gap_width  # Start with a gap before first category
    inner_bars = []
    outer_bars = []
    label_info = []  # Store info for drawing leader lines
    
    for i, tech_cat in enumerate(tech_categories):
        tech_data = outcome_props[outcome_props['technology_simple'] == tech_cat]
        tech_total = tech_totals[tech_cat]
        tech_proportion = tech_total / total_all
        tech_width = tech_proportion * available_space  # Use available space, not full circle
        
        if tech_width < 0.01:  # Skip very small categories
            continue
        
        # Inner ring - technology category
        tech_color = tech_color_map.get(tech_cat, '#CCCCCC')  # Use consistent color mapping
        theta_inner = theta_start + tech_width/2
        
        # Draw inner ring segment
        inner_bar = ax.bar(theta_inner, 1.45, width=tech_width, bottom=0.0,
                          color=tech_color, alpha=0.7, edgecolor='white', linewidth=2)
        inner_bars.append((inner_bar, tech_cat, tech_total, theta_inner))
        
        # Store label info for later drawing with leader lines
        if tech_width > 0.1 and np.isfinite(theta_inner):  # Only label reasonably sized segments
            label_info.append({
                'theta': theta_inner,
                'text': f"{tech_cat}\n({tech_total/total_all:,.2%})",
                'inner_radius': 1.0,
                'color': tech_color
            })
        
        # Outer ring - processing outcomes within this technology
        # Outcomes should span the same width as their parent category (with no internal gaps)
        outcome_start = theta_start
        tech_outcomes = tech_data.sort_values('count', ascending=False)

        if DRAW_OUTCOMES == True:
            for _, row in tech_outcomes.iterrows():
                # Scale outcome width to fit within the technology category width
                outcome_proportion_within_tech = row['count'] / tech_total
                outcome_width = tech_width * outcome_proportion_within_tech
                outcome_color = outcome_colors.get(row['outcome'], '#cccccc')
                outcome_theta = outcome_start + outcome_width/2
                
                # Draw outer ring segment (continuous within category, no gaps between outcomes)
                if outcome_width > 0.005:  # Only draw if visible
                    outer_bar = ax.bar(outcome_theta, 0.65, width=outcome_width, bottom=1.55,
                                       color=outcome_color, alpha=0.9, edgecolor='white', linewidth=0.5)
                
                    # Add outcome labels to large segments
                    if outcome_width > 0.15:  # Only label large outcome segments
                        ax.text(outcome_theta, 1.875, row['outcome'], 
                                ha='center', va='center', fontsize=12, color='white', 
                                fontweight='bold')
                    
                outcome_start += outcome_width
            
        theta_start += tech_width + gap_width  # Add gap after each category
    
    # Draw leader lines and labels outside the chart
    inner_r = 0.75  # Center of inner ring
    line_dist = 0.25 # Leader line extension from inner
    text_dist = 0.45 # Text distance from inner
    lasttheta = -1

    for label in label_info:
        theta = label['theta']
        
        # Check for valid theta
        if not np.isfinite(theta):
            continue

        
        if abs(theta - lasttheta) < 0.628:
            line_dist = -line_dist
            text_dist = -text_dist

        lasttheta = theta
        line_r = inner_r + line_dist  # Leader line extension
        text_r = inner_r + text_dist  # Text position
        
        # Use polar coordinate system properly
        ax.annotate(label['text'], 
                   xy=(theta, inner_r),  # Point at edge of inner ring (theta, radius)
                   xytext=(theta, text_r),  # Text position (theta, radius)
                   xycoords='data',
                   textcoords='data',
                   arrowprops=dict(arrowstyle='-', color='black', lw=2, alpha=0.8),
                   ha='center', va='center', fontweight='bold', fontsize=14,
                   bbox=dict(boxstyle="round,pad=0.3", facecolor='white', alpha=0.9, edgecolor='gray'))
    
    # Create legend for outcomes
    outcome_legend_elements = [plt.Rectangle((0,0),1,1, facecolor=color, label=f"{outcome} (n={df_clean[df_clean['outcome']==outcome].shape[0]:,})") 
                              for outcome, color in outcome_colors.items() 
                              if outcome in df_clean['outcome'].values]
    
    legend = ax.legend(handles=outcome_legend_elements, loc='center left', bbox_to_anchor=(1.2, 0.5),
                      title='Processing Outcomes\n(BB=both good, TB=mate1 bad/mate2 good,\nBT=mate1 good/mate2 bad, TT=both bad,\nB=single good, T=single bad)',
                      fontsize=12, title_fontsize=14)
    legend.get_title().set_fontweight('bold')
    
    # Customize polar plot
    ax.set_theta_zero_location('N')  # Start from top
    ax.set_theta_direction(-1)  # Clockwise
    ax.set_ylim(0, 3.0)  # Fill available space with larger rings
    ax.set_rticks([])  # Remove radial ticks
    ax.set_thetagrids([])  # Remove angular grid
    ax.grid(False)
    
    plt.title('Processing Outcomes by Sequencing Technology\n(Inner ring: Technology, Outer ring: Outcomes)', 
              fontsize=24, fontweight='bold', pad=40, y=1.05)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_path}")
    
    return fig, ax

def generate_summary_report(df, outcome_props):
    """Generate a summary report of the analysis."""
    
    total_accessions = len(df)
    missing_annotations = df['technology'].isna().sum()
    missing_outcomes = df['outcome'].isna().sum()
    
    print("\n" + "="*50)
    print("SEQUENCING TECHNOLOGY ANALYSIS SUMMARY")
    print("="*50)
    print(f"Total accessions processed: {total_accessions:,}")
    print(f"Missing manual annotations: {missing_annotations:,} ({missing_annotations/total_accessions*100:.1f}%)")
    print(f"Missing processing outcomes: {missing_outcomes:,} ({missing_outcomes/total_accessions*100:.1f}%)")
    
    print(f"\nTechnology categories found:")
    tech_counts = df['technology_simple'].value_counts()
    for tech, count in tech_counts.items():
        print(f"  {tech}: {count:,} ({count/total_accessions*100:.1f}%)")
    
    print(f"\nResolution level distribution:")
    res_counts = df['resolution_simple'].value_counts()
    for res, count in res_counts.items():
        print(f"  {res}: {count:,} ({count/total_accessions*100:.1f}%)")
    
    print(f"\nProcessing outcomes distribution:")
    outcome_counts = df['outcome'].value_counts()
    for outcome, count in outcome_counts.items():
        if pd.notna(outcome):
            print(f"  {outcome}: {count:,} ({count/total_accessions*100:.1f}%)")
    
    print(f"\nAccessions needing annotation propagation:")
    no_tech_count = df[df['needs_technology_propagation']].shape[0]
    no_group_count = df[df['needs_group_propagation']].shape[0]
    no_any_count = df[df['needs_any_annotation']].shape[0]
    print(f"  Need technology from group: {no_tech_count:,} ({no_tech_count/total_accessions*100:.1f}%)")
    print(f"  Need group from resolution: {no_group_count:,} ({no_group_count/total_accessions*100:.1f}%)")
    print(f"  Need any annotation: {no_any_count:,} ({no_any_count/total_accessions*100:.1f}%)")

def main():
    """Main execution function."""
    
    parser = argparse.ArgumentParser(description='Create sequencing technology handling plots')
    parser.add_argument('--output', '-o', default='technology_outcomes_plot.svg',
                       help='Output file path for the plot')
    parser.add_argument('--save-data', action='store_true',
                       help='Save processed data to CSV file')
    
    args = parser.parse_args()
    
    # Load and process data
    print("Loading and processing data...")
    df = load_and_process_data()
    
    # Categorize technologies
    print("Categorizing technologies...")
    df = categorize_technologies(df)
    
    # Analyze processing outcomes
    print("Analyzing processing outcomes...")
    outcome_props, df_clean = analyze_processing_outcomes(df)
    
    # Generate summary report
    generate_summary_report(df, outcome_props)
    
    # Create visualization
    print("\nCreating visualization...")
    fig, ax = create_multilevel_pie_chart(outcome_props, df_clean, args.output)
    
    # Create second plot without "Other" category
    print("\nCreating visualization without 'Other' category...")
    df_clean_no_other = df_clean[df_clean['technology_simple'] != 'Other'].copy()
    outcome_props_no_other = outcome_props[outcome_props['technology_simple'] != 'Other'].copy()
    
    if len(df_clean_no_other) > 0:
        output_no_other = args.output.replace('.svg', '_no_other.svg')
        fig2, ax2 = create_multilevel_pie_chart(outcome_props_no_other, df_clean_no_other, output_no_other)
        print(f"Plot without 'Other' category saved to {output_no_other}")
    else:
        print("No data available after removing 'Other' category")
    
    # Optionally save processed data
    if args.save_data:
        output_csv = args.output.replace('.svg', '_processed_data.csv')
        df.to_csv(output_csv, index=False)
        print(f"Processed data saved to {output_csv}")
        
        # Also save summary statistics
        summary_csv = args.output.replace('.svg', '_summary.csv')
        outcome_props.to_csv(summary_csv, index=False)
        print(f"Summary statistics saved to {summary_csv}")
    
    print(f"\nScript completed successfully!")

if __name__ == "__main__":
    main()