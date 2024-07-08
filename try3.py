import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import simpledialog, messagebox
import matplotlib.pyplot as plt

class qPCRApp:
    def __init__(self, root):
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

        # Submit Button
        self.submit_button = tk.Button(self.root, text="Submit", command=self.submit)
        self.submit_button.grid(row=4, column=0, columnspan=2, padx=10, pady=10)

    def submit(self):
        ref_gene = self.ref_gene_entry.get()
        ntc = self.ntc_entry.get()
        water = self.water_entry.get()
        num_groups = self.num_groups_entry.get()

        if not ref_gene or not ntc or not water or not num_groups.isdigit():
            messagebox.showerror("Input Error", "Please fill in all fields correctly.")
            return

        control_samples = [sample.strip() for sample in ntc.split(',')]
        num_groups = int(num_groups)

        self.treated_groups = []

        for i in range(num_groups):
            group_name = simpledialog.askstring("Treated Group", f"Enter the name of treated group {i + 1}:")
            sample_names = simpledialog.askstring("Sample Names", f"Enter the sample names for treated group {i + 1} (comma-separated):")
            self.treated_groups.append((group_name, [sample.strip() for sample in sample_names.split(',')]))

        self.process_data(ref_gene, control_samples, water, self.treated_groups)

    def process_data(self, ref_gene, control_samples, water, treated_groups):
        # Use the raw string method to specify the file path
        file_path = r"C:\Users\liron\Downloads\gzf1 ko imprinting influ new (2).xls"

        # Read the "Results" sheet using xlrd engine and skip the first 47 rows
        df = pd.read_excel(file_path, sheet_name='Results', engine='xlrd', skiprows=47)

        # Print the column names to diagnose the issue
        print("Column names in the DataFrame:")
        print(df.columns.tolist())  # Print column names for debugging

        relevant_columns = ['Sample Name', 'Target Name', 'CT']
        df_relevant = df[relevant_columns]

        # Replace "Undetermined" with 40
        df_relevant['CT'] = df_relevant['CT'].replace('Undetermined', 40)
        # Convert CT to numeric values
        df_relevant['CT'] = pd.to_numeric(df_relevant['CT'], errors='coerce')

        # Remove genes with Ct values over 32
        removed_genes = df_relevant[df_relevant['CT'] > 32]['Target Name'].unique()
        df_relevant = df_relevant[df_relevant['CT'] <= 32]

        if len(removed_genes) > 0:
            messagebox.showwarning("Removed Genes", f"The following genes were removed due to high Ct values: {', '.join(removed_genes)}")

        # Perform sanity checks
        self.check_water_template(df_relevant, water)

        # Calculate ∆Cq
        df_relevant['∆Cq'] = df_relevant.apply(
            lambda row: row['CT'] - df_relevant[(df_relevant['Sample Name'] == row['Sample Name']) & (df_relevant['Target Name'] == ref_gene)]['CT'].values[0],
            axis=1
        )

        # Calculate Expression (2^(-∆Cq))
        df_relevant['Expression'] = 2 ** (-df_relevant['∆Cq'])

        # Calculate Mean and Standard Deviation of Expression
        expression_stats = df_relevant.groupby(['Sample Name', 'Target Name'])['Expression'].agg(['mean', 'std']).reset_index()
        expression_stats.rename(columns={'mean': 'Mean Expression', 'std': 'Expression SD'}, inplace=True)

        # Merge with the main dataframe
        df_relevant = df_relevant.merge(expression_stats, on=['Sample Name', 'Target Name'])

        # Calculate ∆∆Cq for treated vs. control
        control_expression = df_relevant[df_relevant['Sample Name'].isin(control_samples)].groupby('Target Name')['Expression'].mean().reset_index()
        control_expression.rename(columns={'Expression': 'Control Expression'}, inplace=True)
        df_relevant = df_relevant.merge(control_expression, on='Target Name')
        df_relevant['∆∆Cq'] = df_relevant['Expression'] / df_relevant['Control Expression']

        # Calculate % KD (Knockdown percentage)
        df_relevant['% KD'] = (1 - df_relevant['∆∆Cq']) * 100

        # Generate Output Table
        output_columns = ['Sample Name', 'Target Name', 'CT', '∆Cq', 'Expression', 'Mean Expression', 'Expression SD', '∆∆Cq', '% KD']
        output_df = df_relevant[output_columns]
        output_file_path = r"C:\Users\liron\Downloads\qPCR_results_output.xlsx"
        output_df.to_excel(output_file_path, index=False)
        print(f"Output table saved to {output_file_path}")

        # Generate Graph
        self.generate_graph(df_relevant, control_samples)

    def check_water_template(self, df, water):
        water_samples = df[df['Sample Name'] == water]
        unclean_targets = water_samples[water_samples['CT'] <= 35]['Target Name'].tolist()
        if unclean_targets:
            messagebox.showwarning("Water Template Check", f"Water template is not clean for targets: {', '.join(unclean_targets)}")

    def generate_graph(self, df, control_samples):
        plt.figure(figsize=(10, 6))
        control_expression = df[df['Sample Name'].isin(control_samples)].groupby('Target Name')['Expression'].mean()
        treated_groups = df[~df['Sample Name'].isin(control_samples)].groupby(['Target Name', 'Sample Name'])['∆∆Cq'].mean().unstack()

        for target in treated_groups.index:
            plt.plot([target] * len(control_expression), control_expression[target], 'ko', label='Control' if target == treated_groups.index[0] else "")
            for group in treated_groups.columns:
                plt.plot(target, treated_groups.at[target, group], 'o', label=group)

        plt.xlabel('Target Name')
        plt.ylabel('∆∆Cq Expression')
        plt.title('qPCR ∆∆Cq Expression Analysis')
        plt.legend()
        plt.savefig(r"C:\Users\liron\Downloads\qPCR_results_graph.png")
        plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = qPCRApp(root)
    root.mainloop()
