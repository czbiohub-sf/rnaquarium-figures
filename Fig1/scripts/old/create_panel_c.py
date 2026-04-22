#!/usr/bin/env python3
"""
Create Panel C: Pipeline filtering performance by seq-detective classification reason.

This script analyzes pipeline performance metrics across different seq-detective classification 
reasons to show how filtering effectiveness varies by the specific quality assessment criteria.
Uses horizontal violin plots to display the distribution of performance metrics.
"""

import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_and_merge_data():
    """Load stats and seq-detective judgement data, merge on sample ID."""
    
    # Load pipeline statistics
    stats_df = pl.read_csv('data/75k_unstable/stats-merged.csv')
    
    # Load seq-detective judgements
    judgement_df = pl.read_csv('data/75k_unstable/seq-detective-judgement-summary.txt', 
                              separator='\t', has_header=False,
                              new_columns=['sample_id', 'file1', 'file2', 'mate1_judge', 'mate2_judge', 'reason'])
    
    # Extract sample ID from judgement data (remove _subsample suffix)
    judgement_df = judgement_df.with_columns(
        pl.col('sample_id').str.replace('_subsample', '').alias('id')
    )
    
    # Use classification reason directly as category, with some cleaning
    judgement_df = judgement_df.with_columns(
        pl.col('reason').str.to_titlecase().alias('category')
    )
    
    # Merge datasets
    merged_df = stats_df.join(judgement_df.select(['id', 'category']), on='id', how='inner')
    
    return merged_df

def calculate_performance_metrics(df):
    """Calculate key pipeline performance metrics."""
    
    # Calculate all performance metrics using polars expressions
    metrics_df = df.with_columns([
        # FASTP filtering rate
        (pl.col('fastp_reads_after') / pl.col('fastp_reads_before')).alias('fastp_retention'),
        
        # Alignment rates for different aligners
        (pl.col('kallisto_aligned') / pl.col('kallisto_reads_before')).alias('kallisto_alignment_rate'),
        (pl.col('hisat2_aligned') / pl.col('hisat2_reads_before')).alias('hisat2_alignment_rate'),
        (pl.col('star_aligned_unique') / pl.col('star_reads_before')).alias('star_alignment_rate'),
        (pl.col('bowtie2_aligned') / pl.col('bowtie2_reads_before')).alias('bowtie2_alignment_rate'),
        
        # Deduplication rate
        (pl.col('dedup_reads_after') / pl.col('dedup_reads_before')).alias('dedup_retention'),
        
        # Overall pipeline efficiency (final reads / starting reads)
        (pl.col('final_reads') / pl.col('starting_reads')).alias('overall_efficiency')
    ])
    
    # Select only the metrics and category columns
    metrics_cols = [
        'fastp_retention', 'kallisto_alignment_rate', 'hisat2_alignment_rate',
        'star_alignment_rate', 'bowtie2_alignment_rate', 'dedup_retention', 
        'overall_efficiency', 'category'
    ]
    
    return metrics_df.select(metrics_cols)

def create_panel_c_plot(metrics_df):
    """Create Panel C visualization showing pipeline performance by seq-detective classification reason."""
    
    fig, axes = plt.subplots(6, 1, figsize=(15, 24))
    fig.suptitle('Panel C: Pipeline Filtering Performance by Classification Reason', fontsize=16, fontweight='bold')
    
    # Define metrics to plot
    plot_metrics = [
        ('fastp_retention', 'FASTP Read Retention'),
        ('kallisto_alignment_rate', 'Kallisto Alignment Rate'),
        ('hisat2_alignment_rate', 'HISAT2 Alignment Rate'),
        ('star_alignment_rate', 'STAR Alignment Rate'),
        ('dedup_retention', 'Deduplication Retention'),
        ('overall_efficiency', 'Overall Pipeline Efficiency')
    ]
    
    # Get all unique categories and sort by frequency for consistent ordering
    category_counts = metrics_df['category'].value_counts().sort('count', descending=True)
    all_categories = category_counts['category'].to_list()
    
    # Define color mapping based on specified groups
    def get_category_color(category):
        category_lower = category.lower()
        
        # Light blue: Usable Mapping Rate
        if 'usable mapping rate' in category_lower:
            return '#87CEEB'  # light blue
        
        # Pink: Nofeature rate, under 1.2% mapping rate, long read
        elif ('nofeature rate' in category_lower or 
              'under 1.2% mapping rate' in category_lower or 
              'long read' in category_lower):
            return '#FFB6C1'  # light pink
        
        # Dark blue: biological fallback assumption, mate1-mate2 similar, mate2-mate1 similar
        elif ('biological fallback assumption' in category_lower or 
              'mate1-mate2 similar' in category_lower or 
              'mate2-mate1 similar' in category_lower):
            return '#000080'  # dark blue
        
        # Red: mate1 technical, mate2 technical, mates < 9% mapping rate, sc-like readlen
        elif ('mate1 technical' in category_lower or 
              'mate2 technical' in category_lower or 
              'mates < 9% mapping rate' in category_lower or 
              'sc-like readlen' in category_lower):
            return '#DC143C'  # crimson red
        
        # Default gray for any unmatched categories
        else:
            return '#696969'
    
    category_colors = {cat: get_category_color(cat) for cat in all_categories}
    
    for i, (metric, title) in enumerate(plot_metrics):
        ax = axes[i]
        
        # Prepare data for horizontal violin plot
        data_for_plot = []
        categories = []
        
        for category in all_categories:
            category_data = metrics_df.filter(pl.col('category') == category)[metric].drop_nulls().to_numpy()
            if len(category_data) > 0:
                data_for_plot.append(category_data)
                categories.append(category)
        
        if data_for_plot:
            # Create horizontal violin plot
            parts = ax.violinplot(data_for_plot, positions=range(len(categories)), 
                                vert=False, showmeans=True, showmedians=True)
            
            # Color the violin plots
            for j, pc in enumerate(parts['bodies']):
                if j < len(categories):
                    color = category_colors.get(categories[j], '#696969')
                    pc.set_facecolor(color)
                    pc.set_alpha(0.7)
            
            # Set y-axis labels (categories) and rotate for readability
            ax.set_yticks(range(len(categories)))
            ax.set_yticklabels(categories, fontsize=10)
            ax.invert_yaxis()  # Put highest frequency category at top
            
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Rate' if 'rate' in metric.lower() else 'Retention/Efficiency')
        ax.grid(True, alpha=0.3, axis='x')
        ax.set_xlim(0, 1.1)
    
    plt.tight_layout()
    return fig

def generate_summary_statistics(metrics_df):
    """Generate summary statistics table."""
    
    summary_stats = metrics_df.group_by('category').agg([
        pl.col('fastp_retention').count().alias('count'),
        pl.col('fastp_retention').mean().round(3).alias('fastp_retention_mean'),
        pl.col('fastp_retention').std().round(3).alias('fastp_retention_std'),
        pl.col('kallisto_alignment_rate').mean().round(3).alias('kallisto_alignment_rate_mean'),
        pl.col('kallisto_alignment_rate').std().round(3).alias('kallisto_alignment_rate_std'),
        pl.col('overall_efficiency').mean().round(3).alias('overall_efficiency_mean'),
        pl.col('overall_efficiency').std().round(3).alias('overall_efficiency_std')
    ])
    
    return summary_stats

def main():
    """Main execution function."""
    
    print("Creating Panel C: Pipeline filtering performance by seq-detective classification reason")
    
    # Load and merge data
    print("Loading data...")
    merged_df = load_and_merge_data()
    print(f"Loaded {len(merged_df)} samples with classification reasons")
    
    # Calculate performance metrics
    print("Calculating performance metrics...")
    metrics_df = calculate_performance_metrics(merged_df)
    
    # Print classification reason distribution
    print("\nClassification reason distribution:")
    category_counts = metrics_df['category'].value_counts().sort('count', descending=True)
    print(category_counts)
    
    # Create visualization
    print("Creating visualization...")
    fig = create_panel_c_plot(metrics_df)
    
    # Save plot
    output_path = Path('figures/panel_c_pipeline_performance_by_classification_reason.png')
    output_path.parent.mkdir(exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved plot to {output_path}")
    
    # Generate and save summary statistics
    print("Generating summary statistics...")
    summary_stats = generate_summary_statistics(metrics_df)
    summary_path = Path('data/panel_c_summary_stats.csv')
    summary_stats.write_csv(summary_path)
    print(f"Saved summary statistics to {summary_path}")
    
    print("\nSummary Statistics by Classification Reason:")
    print(summary_stats)
    
    plt.show()

if __name__ == "__main__":
    main()