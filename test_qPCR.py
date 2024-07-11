import pandas as pd
import matplotlib.pyplot as plt
from qPCR_data_analysis_with_errorbars import qPCRApp  
import os

# Test function
def test_qpcr_analysis():
    # Setup
    excel_path = "qPCR_results.xls"  
    output_path = "test_output"
    skiplines = 47
    config_path = "config_file.yaml"
    show_plot = False

    # Create qPCRApp instance
    app = qPCRApp(None, excel_path, output_path, skiplines, gui=False, show_plot=show_plot)

    # Process data using config file
    app.process_from_config(config_path)

    # Check if output files are created
    table_path = output_path + "_table.csv"
    plot_path = output_path + ".jpg"
    
    assert pd.read_csv(table_path).shape[0] > 0, "Output table is empty"
    assert os.path.exists(plot_path), "Plot file was not created"

    # Check calculations
    df = pd.read_csv(table_path)
    
    # Check if all expected columns are present
    expected_columns = ['Treatment', 'Sample', 'Target Gene', 'Cq Reference Gene', 'Cq Target', 
                        '∆Cq', '∆Cq Expression', 'Mean ∆Cq Expression', '∆Cq Expression stdev', 
                        '∆∆Cq Expression', '∆∆Cq Expression stdev', '% KD']
    assert all(col in df.columns for col in expected_columns), "Some expected columns are missing"

    # Check if calculations are correct for a sample row
    sample_row = df.iloc[0]
    delta_cq = sample_row['Cq Target'] - sample_row['Cq Reference Gene']
    assert abs(sample_row['∆Cq'] - delta_cq) < 1e-6, "∆Cq calculation is incorrect"

    delta_cq_expression = 2 ** (-delta_cq)
    assert abs(sample_row['∆Cq Expression'] - delta_cq_expression) < 1e-6, "∆Cq Expression calculation is incorrect"

    control_expression = df[df['Treatment'] == 'Control']['Mean ∆Cq Expression'].values[0]
    delta_delta_cq_expression = sample_row['Mean ∆Cq Expression'] / control_expression
    assert abs(sample_row['∆∆Cq Expression'] - delta_delta_cq_expression) < 1e-6, "∆∆Cq Expression calculation is incorrect"

    kd_percentage = (1 - delta_delta_cq_expression) * 100
    assert abs(sample_row['% KD'] - kd_percentage) < 1e-6, "% KD calculation is incorrect"

