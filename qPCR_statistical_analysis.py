import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

def perform_statistical_analysis(input_file, output_prefix):
    # Read the CSV file produced by the main qPCR analysis script
    df = pd.read_csv(input_file)
    
    # Rename columns to remove special characters
    df = df.rename(columns={
        '∆∆Cq Expression': 'ddCq_Expression',
        '∆∆Cq Expression stdev': 'ddCq_Expression_stdev'
    })
    
    # Group the data by Target Gene
    grouped = df.groupby('Target Gene')
    
    for target_gene, group in grouped:
        # Perform one-way ANOVA
        model = ols('ddCq_Expression ~ C(Treatment)', data=group).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)
        anova_p_value = anova_table['PR(>F)'][0]
        
        # Perform t-tests (each treatment vs control)
        control_data = group[group['Treatment'] == 'Control']['ddCq_Expression']
        p_values = []
        for treatment in group['Treatment'].unique():
            if treatment != 'Control':
                treatment_data = group[group['Treatment'] == treatment]['ddCq_Expression']
                t_stat, p_value = stats.ttest_ind(control_data, treatment_data)
                p_values.append((treatment, p_value))
        
        # Create the plot
        plot_statistical_results(group, target_gene, p_values, anova_p_value, output_prefix)
        
        # Save statistical results
        save_statistical_results(target_gene, anova_p_value, p_values, output_prefix)

def plot_statistical_results(data, target_gene, p_values, anova_p_value, output_prefix):
    plt.figure(figsize=(10, 6))
    
    treatments = data['Treatment'].unique()
    expressions = data.groupby('Treatment')['ddCq_Expression'].mean()
    errors = data.groupby('Treatment')['ddCq_Expression_stdev'].mean()
    
    bars = plt.bar(treatments, expressions, yerr=errors, capsize=5)
    plt.title(f'Relative Expression of {target_gene}\nANOVA p-value: {anova_p_value:.4f}')
    plt.xlabel('Treatment Groups')
    plt.ylabel('∆∆Cq Expression')
    plt.xticks(rotation=45, ha='right')

    # Add significance stars
    for i, treatment in enumerate(treatments):
        if treatment != 'Control':
            p_value = dict(p_values).get(treatment)
            if p_value:
                height = expressions[treatment] + errors[treatment]
                plt.text(i, height, get_significance_symbol(p_value),
                         ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_{target_gene}_statistical_plot.png")
    plt.close()

def save_statistical_results(target_gene, anova_p_value, p_values, output_prefix):
    with open(f"{output_prefix}_{target_gene}_statistical_results.txt", 'w') as f:
        f.write(f"Statistical Results for {target_gene}\n")
        f.write(f"ANOVA p-value: {anova_p_value:.4f}\n\n")
        f.write("T-test results (vs Control):\n")
        for treatment, p_value in p_values:
            f.write(f"{treatment}: p-value = {p_value:.4f} {get_significance_symbol(p_value)}\n")

def get_significance_symbol(p_value):
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    else:
        return "ns"

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) != 3:
        print("Usage: python statistical_analysis.py <input_csv_file> <output_prefix>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    perform_statistical_analysis(input_file, output_prefix)