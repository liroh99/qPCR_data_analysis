import pandas as pd
import tkinter as tk
from tkinter import simpledialog, messagebox
import matplotlib.pyplot as plt
import sys
import numpy as np
import argparse
import yaml

class qPCRApp:

    def __init__(self, root, excel_path, output_path, skiplines, gui,show_plot):
        self.excel_path = excel_path
        self.output_path = output_path
        self.skiplines = skiplines
        self.gui = gui
        self.show_plot=show_plot

        if gui:
            self.root = root
            self.root.title("qPCR Data Analysis")
            self.create_widgets()

    def create_widgets(self):
        # Reference Gene
        self.ref_gene_label = tk.Label(self.root, text="Reference Gene:")
        self.ref_gene_label.grid(row=0, column=0, padx=10, pady=10)
        self.ref_gene_entry = tk.Entry(self.root)
        self.ref_gene_entry.grid(row=0, column=1, padx=10, pady=10)

        # Control Sample Names (NTC)
        self.ntc_label = tk.Label(self.root, text="Control Sample Names (NTC, comma-separated):")
        self.ntc_label.grid(row=1, column=0, padx=10, pady=10)
        self.ntc_entry = tk.Entry(self.root)
        self.ntc_entry.grid(row=1, column=1, padx=10, pady=10)

        # Water Template Sample Name
        self.water_label = tk.Label(self.root, text="Water Template Sample Name:")
        self.water_label.grid(row=2, column=0, padx=10, pady=10)
        self.water_entry = tk.Entry(self.root)
        self.water_entry.grid(row=2, column=1, padx=10, pady=10)

        # Number of Treated Groups
        self.num_groups_label = tk.Label(self.root, text="Number of Treated Groups:")
        self.num_groups_label.grid(row=3, column=0, padx=10, pady=10)
        self.num_groups_entry = tk.Entry(self.root)
        self.num_groups_entry.grid(row=3, column=1, padx=10, pady=10)

        # Replicates
        self.replicates_label = tk.Label(self.root, text="Replicates (2 for duplicates, 3 for triplicates):")
        self.replicates_label.grid(row=4, column=0, padx=10, pady=10)
        self.replicates_entry = tk.Entry(self.root)
        self.replicates_entry.grid(row=4, column=1, padx=10, pady=10)

        # Submit Button
        self.submit_button = tk.Button(self.root, text="Submit", command=self.submit)
        self.submit_button.grid(row=5, column=0, columnspan=2, padx=10, pady=10)

    def submit(self):
        ref_gene = self.ref_gene_entry.get()
        ntc = self.ntc_entry.get()
        water = self.water_entry.get()
        num_groups = self.num_groups_entry.get()
        replicates = self.replicates_entry.get()

        if not ref_gene or not ntc or not water or not num_groups.isdigit() or not replicates.isdigit():
            messagebox.showerror("Input Error", "Please fill in all fields correctly.")
            return

        control_samples = [sample.strip() for sample in ntc.split(',')]
        num_groups = int(num_groups)
        replicates = int(replicates)

        self.treated_groups = []

        for i in range(num_groups):
            group_name = simpledialog.askstring("Treated Group", f"Enter the name of treated group {i + 1}:")
            sample_names = simpledialog.askstring("Sample Names", f"Enter the sample names for treated group {i + 1} (comma-separated):")
            self.treated_groups.append((group_name, [sample.strip() for sample in sample_names.split(',')]))

        self.process_data(ref_gene, control_samples, water, self.treated_groups, replicates)


    def process_from_config(self, config_path):
        with open(config_path) as stream:
            config = yaml.safe_load(stream)

        ref_gene = config["ref_gene"]
        ntc = config["ntc"]
        water = config["water"]
        replicates = config["replicates"]
        self.treated_groups = []
        control_samples = [sample.strip() for sample in ntc.split(',')]

        for group in config["treated_groups"]:
            group_name = group["name"] 
            sample_names = group["samples"]
            self.treated_groups.append((group_name, [sample.strip() for sample in sample_names.split(',')]))

        self.process_data(ref_gene, control_samples, water, self.treated_groups, replicates)


    def process_data(self, ref_gene, control_samples, water, treated_groups, replicates):
        try:
            df = pd.read_excel(self.excel_path, sheet_name='Results', skiprows=self.skiplines)
        except Exception as e:
            print(f"Error reading the Excel file: {e}")
            sys.exit(1)

        relevant_columns = ['Sample Name', 'Target Name', 'CT']
        df_relevant = df[relevant_columns].copy()

        # Remove rows with "Undetermined" CT values
        df_relevant = df_relevant[df_relevant['CT'] != 'Undetermined']

        # Convert CT to numeric values
        df_relevant['CT'] = pd.to_numeric(df_relevant['CT'], errors='coerce')

        # Remove rows with NaN CT values (if any left after conversion)
        df_relevant = df_relevant.dropna(subset=['CT'])

        # Perform sanity checks
        self.check_water_template(df_relevant, water)

        # Ignore water template samples
        df_relevant = df_relevant[df_relevant['Sample Name'] != water]

        # Calculate mean Ct for each sample and target combination
        df_mean_ct = df_relevant.groupby(['Sample Name', 'Target Name'])['CT'].mean().reset_index()
        df_mean_ct.rename(columns={'CT': 'Ct Mean'}, inplace=True)

        # Prepare the data for table creation
        self.prepare_data(ref_gene, control_samples, df_mean_ct, treated_groups)

    def prepare_data(self, ref_gene, control_samples, df_mean_ct, treated_groups):
        all_tables = []
        plot_data = []

        for target_gene in df_mean_ct['Target Name'].unique():
            if target_gene == ref_gene:
                continue

            table = []
            for group_name, samples in treated_groups + [('Control', control_samples)]:
                group_data = []
                for sample in samples:
                    ref_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == ref_gene)]['Ct Mean']
                    target_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == target_gene)]['Ct Mean']
                    
                    if ref_ct_mean.empty or target_ct_mean.empty:
                        continue
                    
                    ref_ct_mean = ref_ct_mean.values[0]
                    target_ct_mean = target_ct_mean.values[0]
                    
                    deltaCq = target_ct_mean - ref_ct_mean
                    deltaCq_expression = 2 ** (-deltaCq)
                    group_data.append({
                        'Treatment': group_name,
                        'Sample': sample,
                        'Cq Reference Gene': ref_ct_mean,
                        'Cq Target': target_ct_mean,
                        'Target Gene': target_gene,
                        '∆Cq': deltaCq,
                        '∆Cq Expression': deltaCq_expression
                    })
                
                if group_data:
                    df_group = pd.DataFrame(group_data)
                    mean_deltaCq_expression = df_group['∆Cq Expression'].mean()
                    stdev_deltaCq_expression = df_group['∆Cq Expression'].std()
                    
                    for row in df_group.to_dict('records'):
                        row.update({
                            'Mean ∆Cq Expression': mean_deltaCq_expression,
                            '∆Cq Expression stdev': stdev_deltaCq_expression
                        })
                        table.append(row)

            df_table = pd.DataFrame(table)

            # Calculate ∆∆Cq expression and %KD
            control_mean_deltaCq_expression = df_table[df_table['Treatment'] == 'Control']['Mean ∆Cq Expression'].values[0]

            df_table['∆∆Cq Expression'] = df_table['Mean ∆Cq Expression'] / control_mean_deltaCq_expression
            df_table['∆∆Cq Expression stdev'] = df_table['∆Cq Expression stdev'] / control_mean_deltaCq_expression
            df_table['% KD'] = (1 - df_table['∆∆Cq Expression']) * 100

            all_tables.append(df_table)

            # Prepare data for plotting
            plot_data.append({
                'Target Gene': target_gene,
                'Treatment': df_table['Treatment'].tolist(),
                '∆∆Cq Expression': df_table['∆∆Cq Expression'].tolist(),
                '∆∆Cq Expression stdev': df_table['∆∆Cq Expression stdev'].tolist()
            })

        # Concatenate all tables
        final_table = pd.concat(all_tables, ignore_index=True)
        self.save_results(final_table)

        # Plot the data
        self.plot_data(plot_data)

    def save_results(self, final_table):
        table_path = self.output_path + "_table.csv"
        columns_order = ['Treatment', 'Sample', 'Target Gene', 'Cq Reference Gene', 'Cq Target', '∆Cq', 
                         '∆Cq Expression', 'Mean ∆Cq Expression', '∆Cq Expression stdev', 
                         '∆∆Cq Expression', '∆∆Cq Expression stdev', '% KD']
        final_table = final_table[columns_order]
        final_table.to_csv(table_path, index=False)
        print(f"Table saved to {table_path}")

    def plot_data(self, plot_data):
        treatment_groups = set()
        target_genes = set()
        
        for data in plot_data:
            treatment_groups.update(data['Treatment'])
            target_genes.add(data['Target Gene'])

        treatment_groups = sorted(treatment_groups)
        target_genes = sorted(target_genes)

        # Create DataFrames to hold the plotting data and error bars
        plot_df = pd.DataFrame(columns=['Treatment'] + target_genes)
        error_df = pd.DataFrame(columns=['Treatment'] + target_genes)

        for treatment in treatment_groups:
            row = {'Treatment': treatment}
            error_row = {'Treatment': treatment}
            for data in plot_data:
                if treatment in data['Treatment']:
                    index = data['Treatment'].index(treatment)
                    row[data['Target Gene']] = data['∆∆Cq Expression'][index]
                    error_row[data['Target Gene']] = data['∆∆Cq Expression stdev'][index]
            plot_df = pd.concat([plot_df, pd.DataFrame([row])], ignore_index=True)
            error_df = pd.concat([error_df, pd.DataFrame([error_row])], ignore_index=True)

        plot_df.set_index('Treatment', inplace=True)
        error_df.set_index('Treatment', inplace=True)
        plot_df = plot_df.astype(float)
        error_df = error_df.astype(float)

        # Plotting
        ax = plot_df.plot(kind='bar', yerr=error_df, figsize=(12, 6), capsize=5)
        ax.set_title('Relative Gene Expression Across Treatments')
        ax.set_xlabel('Treatment Groups')
        ax.set_ylabel('∆∆Cq Expression')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.legend(title='Genes')
        plot_path = self.output_path + ".jpg"
        plt.savefig(plot_path, format='jpeg')
        print(f"Plot saved to {plot_path}")

        if show_plot:
            plt.show()

    def check_water_template(self, df, water):
        water_samples = df[df['Sample Name'] == water]
        unclean_targets = water_samples[water_samples['CT'] <= 35]['Target Name'].tolist()
        if unclean_targets and self.gui:
            messagebox.showwarning("Water Template Check", f"Water template is not clean for targets: {', '.join(unclean_targets)}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Input your qPCR excel file for differential gene analysis")
    parser.add_argument('--input_path', '-i',  help="qPCR raw data", required=True)
    parser.add_argument('--output_path', '-o', help="Output file that includes all the calculations", required=True )
    parser.add_argument('--skiplines', type=int, default=47, help="How many lines to skip to reach the column names?")
    parser.add_argument('--config_path', '-config', help="path to yaml file to override gui")
    parser.add_argument('--show_plot','-show',action="store_true", help="Show plot for analysis")

    args = parser.parse_args()

    excel_path = args.input_path
    output_path = args.output_path
    skiplines = args.skiplines
    show_plot=args.show_plot

    if args.config_path:
        root = None
        gui = False
    else:
        root = tk.Tk()
        gui = True

    app = qPCRApp(root, excel_path, output_path, skiplines, gui,show_plot)
    if args.config_path:
        app.process_from_config(args.config_path)
    else: 
        root.mainloop()
