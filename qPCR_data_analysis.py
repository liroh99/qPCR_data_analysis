import pandas as pd
import tkinter as tk
from tkinter import simpledialog, messagebox
import matplotlib.pyplot as plt
import sys

class qPCRApp:

    def __init__(self, root, excel_path, output_path):
        self.excel_path = excel_path
        self.output_path = output_path

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

    def process_data(self, ref_gene, control_samples, water, treated_groups, replicates):
        try:
            df = pd.read_excel(self.excel_path, sheet_name='Results', skiprows=47)
        except Exception as e:
            print(f"Error reading the Excel file: {e}")
            sys.exit(1)

        relevant_columns = ['Sample Name', 'Target Name', 'CT']
        df_relevant = df[relevant_columns].copy()

        # Replace "Undetermined" with 40
        df_relevant['CT'] = df_relevant['CT'].replace('Undetermined', 40)
        # Convert CT to numeric values
        df_relevant['CT'] = pd.to_numeric(df_relevant['CT'], errors='coerce')

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
        all_tables = []  # To store all tables
        plot_data = []  # To store data for plotting

        for target_gene in df_mean_ct['Target Name'].unique():
            if target_gene == ref_gene:
                continue

            table = []

            for group_name, samples in treated_groups:
                for sample in samples:
                    ref_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == ref_gene)]['Ct Mean']
                    target_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == target_gene)]['Ct Mean']
                    
                    if ref_ct_mean.empty or target_ct_mean.empty:
                        continue
                    
                    ref_ct_mean = ref_ct_mean.values[0]
                    target_ct_mean = target_ct_mean.values[0]
                    
                    deltaCq = target_ct_mean - ref_ct_mean
                    deltaCq_expression = 2 ** (-deltaCq)
                    table.append({
                        'Treatment': group_name,
                        'Sample': sample,
                        'Cq Reference Gene': ref_ct_mean,
                        'Cq Target Gene': target_ct_mean,
                        'Target Gene': target_gene,  # Add target gene name
                        'deltaCq': deltaCq,
                        'deltaCq Expression': deltaCq_expression
                    })

            for sample in control_samples:
                ref_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == ref_gene)]['Ct Mean']
                target_ct_mean = df_mean_ct[(df_mean_ct['Sample Name'] == sample) & (df_mean_ct['Target Name'] == target_gene)]['Ct Mean']
                
                if ref_ct_mean.empty or target_ct_mean.empty:
                    continue
                
                ref_ct_mean = ref_ct_mean.values[0]
                target_ct_mean = target_ct_mean.values[0]
                
                deltaCq = target_ct_mean - ref_ct_mean
                deltaCq_expression = 2 ** (-deltaCq)
                table.append({
                    'Treatment': 'Control',
                    'Sample': sample,
                    'Cq Reference Gene': ref_ct_mean,
                    'Cq Target Gene': target_ct_mean,
                    'Target Gene': target_gene,  # Add target gene name
                    'deltaCq': deltaCq,
                    'deltaCq Expression': deltaCq_expression
                })

            df_table = pd.DataFrame(table)

            # Calculate mean deltaCq expression for each treatment
            mean_deltaCq_expression = df_table.groupby('Treatment')['deltaCq Expression'].mean().reset_index()
            mean_deltaCq_expression.rename(columns={'deltaCq Expression': 'Mean deltaCq Expression'}, inplace=True)

            # Merge with the original table
            df_table = df_table.merge(mean_deltaCq_expression, on='Treatment')

            # Calculate delta deltaCq expression and %KD
            control_mean_deltaCq_expression = mean_deltaCq_expression[mean_deltaCq_expression['Treatment'] == 'Control']['Mean deltaCq Expression'].values[0]
            df_table['delta deltaCq Expression'] = df_table['Mean deltaCq Expression'] / control_mean_deltaCq_expression
            df_table['% KD'] = (1 - df_table['delta deltaCq Expression']) * 100

            all_tables.append(df_table)  # Append the table to all_tables

            # Prepare data for plotting
            plot_data.append({
                'Target Gene': target_gene,
                'Treatment': df_table['Treatment'].tolist(),
                'delta deltaCq Expression': df_table['delta deltaCq Expression'].tolist()
            })

        # Concatenate all tables
        final_table = pd.concat(all_tables, ignore_index=True)
        self.save_results(final_table)

        # Plot the data
        self.plot_data(plot_data)

    def save_results(self, final_table):
        table_path = self.output_path + "_table.csv"
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

        # Create a DataFrame to hold the plotting data
        plot_df = pd.DataFrame(columns=['Treatment'] + target_genes)

        for treatment in treatment_groups:
            row = {'Treatment': treatment}
            for data in plot_data:
                if treatment in data['Treatment']:
                    index = data['Treatment'].index(treatment)
                    row[data['Target Gene']] = data['delta deltaCq Expression'][index]
            plot_df = pd.concat([plot_df, pd.DataFrame([row])], ignore_index=True)

        plot_df.set_index('Treatment', inplace=True)
        plot_df = plot_df.astype(float)

        # Plotting
        ax = plot_df.plot(kind='bar', figsize=(12, 6))
        ax.set_title('Relative Gene Expression Across Treatments')
        ax.set_xlabel('Treatment Groups')
        ax.set_ylabel('Delta Delta Cq Expression')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.legend(title='Genes')
        plt.show()
        plot_path = self.output_path + ".png"
        plt.savefig(plot_path)
        print(f"Plot saved to {plot_path}")

    def check_water_template(self, df, water):
        water_samples = df[df['Sample Name'] == water]
        unclean_targets = water_samples[water_samples['CT'] <= 35]['Target Name'].tolist()
        if unclean_targets:
            messagebox.showwarning("Water Template Check", f"Water template is not clean for targets: {', '.join(unclean_targets)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <path_to_excel_file> <path_to_save_table_and_png>")
        sys.exit(1)
    excel_path = sys.argv[1]
    output_path = sys.argv[2]

    root = tk.Tk()
    app = qPCRApp(root, excel_path, output_path)
    root.mainloop()
