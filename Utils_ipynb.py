import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns

def append_csv_files_from_subfolders(root_folder):
    all_dfs = []  # List to hold all the dataframes

    for subdir, dirs, files in os.walk(root_folder):
        for file in files:
            if file == 'cell_stat.csv':
                file_path = os.path.join(subdir, file)
                # Read the file with header and store the header
                temp_df = pd.read_csv(file_path, header=0)
                all_dfs.append(temp_df)

    # Concatenate all dataframes in the list
    master_df = pd.concat(all_dfs, ignore_index=True)

    return master_df

def count_occurrences(dataframe, column_name, search_string):
    """
    Counts the occurrences of a search string in a specified column of a pandas DataFrame.

    Parameters:
    dataframe (pandas.DataFrame): The DataFrame to search in.
    column_name (str): The name of the column to search the string in.
    search_string (str): The string to search for.

    Returns:
    int: The count of occurrences of the search string in the specified column.
    """

    if column_name not in dataframe.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame")

    count = dataframe[column_name].str.contains(search_string, na=False).sum()
    print(f"Number of occurrences of '{search_string}' in the {column_name} column: {count}")
    return count


def normalize_image_name(name):
    # Define a normalization function if needed, e.g., strip file extensions
    return name.split('.czi')[0]


def remove_specific_rows(combined_df, excel_path):
    # Read the Excel file
    cells_to_remove_df = pd.read_excel(excel_path)

    # Normalize the Image_name columns in both dataframes
    combined_df['Image_name_normalized'] = combined_df['Image_name'].apply(normalize_image_name)
    cells_to_remove_df['Image_name_normalized'] = cells_to_remove_df['Image_name']

    # Perform a left merge with an indicator to identify rows to be removed
    merged_df = combined_df.merge(cells_to_remove_df[['Image_name_normalized', 'Cell_num']],
                                  on=['Image_name_normalized', 'Cell_num'],
                                  how='left',
                                  indicator=True)

    # Filter rows where '_merge' is 'left_only' (i.e., not to be removed)
    clean_df = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge', 'Image_name_normalized'])

    return clean_df

def extract_and_append_image_info(df, image_name_col='Image_name'):
    """
    Extracts specific information from the 'Image_name' column of a DataFrame and appends
    this information as new columns to the DataFrame, including 'Base_image_name'.

    The function first removes all white spaces from the 'Image_name' entries. Then, it parses
    each entry to extract the following details:
    - 'Base_image_name': Base name of the image
    - 'Processing': Image processing type (e.g., 'RAW', 'LSM')
    - 'Date': Date of the image capture
    - 'Time': Time point of the experiment (0, 24, 48 hours)
    - 'Cell Type': Type of the cell (e.g., 'KASH', 'KASH+doxy', 'MSC')
    - 'LIV': Indicates whether LIV is present (+LIV) or absent (-LIV)
    - 'Cisp': Cisplatin concentration or 'Control' if not present

    Parameters:
    df (pandas.DataFrame): The DataFrame containing the 'Image_name' column.
    image_name_col (str): The name of the column containing image names. Default is 'Image_name'.

    Returns:
    pandas.DataFrame: The original DataFrame with new columns appended immediately after the 'Image_name' column.
    """

    def extract_info(image_name):
        # Remove all white spaces
        image_name = image_name.replace(" ", "")

        # Extract 'Base_image_name'
        base_image_name = image_name.split("-LSM")[0]

        # Extract 'Processing'
        processing = 'LSM' if 'LSM' in image_name else 'RAW'

        # Extract 'Date'
        date_pattern = re.compile(r'\d+-\d+-\d+')
        date = date_pattern.search(image_name).group(0)

        # Extract 'Time'
        time_pattern = re.compile(r'0hr|24hr|48hr')
        time_match = time_pattern.search(image_name)
        time = time_match.group(0)[:-2] if time_match else 'Unknown'

        # Extract 'Cell Type'
        cell_type_pattern = re.compile(r'KASH(\+doxy)?|MSC')
        cell_type_match = cell_type_pattern.search(image_name)
        cell_type = cell_type_match.group(0) if cell_type_match else 'Unknown'

        # Extract 'LIV'
        liv = '+LIV' if '+LIV' in image_name else '-LIV'

        # Extract 'Cisp'
        cisp_pattern = re.compile(r'(\d+uM)|(-Cisplatin)')
        cisp_match = cisp_pattern.search(image_name)
        cisp = cisp_match.group(0) if cisp_match else 'Control'

        return pd.Series([base_image_name, processing, date, time, cell_type, liv, cisp],
                         index=['Base_image_name', 'Processing', 'Date', 'Time', 'Cell Type', 'LIV', 'Cisp'])

    # Apply the extraction function and create a temporary DataFrame
    temp_df = df[image_name_col].apply(extract_info)

    # Find the index where new columns should be inserted
    insert_loc = df.columns.get_loc(image_name_col) + 1

    # Insert the new columns into the DataFrame at the specified location
    for col in temp_df.columns:
        df.insert(loc=insert_loc, column=col, value=temp_df[col])
        insert_loc += 1

    return df



def filter_dataframe(df, processing=None, date=None, time=None, cell_type=None, liv=None, cisp=None):
    """
    Filters the DataFrame based on the specified criteria. Each parameter is optional.
    If a parameter is not specified, the filter for that criterion will not be applied.

    Parameters:
    df (pandas.DataFrame): The DataFrame to be filtered.
    processing, date, time, cell_type, liv, cisp (str): Filter criteria for the respective columns.

    Returns:
    pandas.DataFrame: The filtered DataFrame.
    """

    if processing is not None:
        df = df[df['Processing'] == processing]
    if date is not None:
        df = df[df['Date'] == date]
    if time is not None:
        df = df[df['Time'] == time]
    if cell_type is not None:
        df = df[df['Cell Type'] == cell_type]
    if liv is not None:
        df = df[df['LIV'] == liv]
    if cisp is not None:
        df = df[df['Cisp'] == cisp]

    return df


def get_correlation_matrices_plt(dataframes):
    """
    Plots correlation matrices for a list of DataFrames using specified columns.

    Parameters:
    dataframes (list of pandas.DataFrame): A list of DataFrames for which correlation matrices will be plotted.
    correlation_cols (list of str): Columns to be used for computing the correlation matrix.
    """
    # Define the columns for correlation matrix
    correlation_cols = [
        'Nucleus_volume',
        'Nucleus_length, micrometre', 'Nucleus_width, micrometre',
        'Nucleus_high, micrometre', 'Cy5-T1 av_signal_in_nuc_area_3D',
        'Cy5-T1 sum_pix_in_nuc_cylinder',
        'Cy5-T1 ring intensity coef', 'AF594-T2 av_signal_in_nuc_area_3D',
        'AF594-T2 sum_pix_in_nuc_cylinder',
        'AF594-T2 ring intensity coef', 'AF488-T3 av_signal_in_nuc_area_3D',
        'AF488-T3 sum_pix_in_nuc_cylinder',
        'AF488-T3 ring intensity coef'
    ]

    # Determine plot layout
    num_dfs = len(dataframes)
    nrows = int(num_dfs**0.5) + 1
    ncols = (num_dfs // nrows) + (num_dfs % nrows > 0)

    # Create subplots
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 8, nrows * 8))

    # Flatten axes array
    axes = axes.flatten()

    # Define text sizes
    title_size = 10  # Size of the title text
    label_size = 9  # Size of the label text
    annot_size = 6  # Size of the annotation text

    # Plot each DataFrame's correlation matrix
    for i, df in enumerate(dataframes):
        ax = sns.heatmap(df[correlation_cols].corr(), ax=axes[i], annot=True, fmt=".2f", cmap='coolwarm', annot_kws={'size': annot_size})
        axes[i].set_title(f'DF {i+1} Correlation Matrix', fontsize=title_size)
        ax.figure.axes[-1].yaxis.label.set_size(label_size)  # Colorbar label size

        # Adjusting the size of the tick labels
        for label in (axes[i].get_xticklabels() + axes[i].get_yticklabels()):
            label.set_fontsize(label_size)

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    return plt

def refine_and_filter_data(df, pixel_size):
    """
    Refines and filters the DataFrame for more efficient analysis. The function removes irrelevant columns,
    creates new columns for deeper insights, renames columns for clarity, performs specific calculations
    related to the nucleus volume and signal intensities, and adds a 'Group' column.

    Parameters:
    df (pandas.DataFrame): The DataFrame to be processed.
    pixel_size (float): The size of a pixel in cubic microns, used for volume calculation.

    Returns:
    pandas.DataFrame: The refined and filtered DataFrame with added calculations.
    """

    # Calculate Nucleus Volume in pixels and additional signal-related columns
    df['Nucleus_volume, pixels'] = (df['Nucleus_volume, cubic_micrometre'] / pixel_size).astype('float64')
    df['AF488-T3 total_signal_in_nuc_area_3D'] = (df['AF488-T3 av_signal_in_nuc_area_3D'] * df['Nucleus_volume, pixels']).astype('float64')
    df['AF488-T3 av_signal_in_nuc_cylinder'] = (df['AF488-T3 sum_pix_in_nuc_cylinder'] / df['Nucleus_cylinder, pixels_number']).astype('float64')

    # Remove unwanted columns
    columns_to_remove = ['Image_name', 'Nucleus_cylinder, pixels_number', 'Nucleus_length, micrometre',
                         'Nucleus_width, micrometre', 'Cy5-T1 av_signal_in_nuc_area_3D',
                         'Cy5-T1 sum_pix_in_nuc_cylinder', 'Cy5-T1 has ring', 'Cy5-T1 ring intensity coef',
                         'AF594-T2 sum_pix_in_nuc_cylinder',
                         'AF594-T2 has ring', 'AF488-T3 sum_pix_in_nuc_cylinder', 'AF488-T3 has ring',
                         'AF488-T3 ring intensity coef', 'Nucleus_volume, pixels', 'AF488-T3 av_signal_in_nuc_cylinder']
    df = df.drop(columns=columns_to_remove)

    # Rename columns
    columns_to_rename = {
        'Nucleus_volume, cubic_micrometre': 'Nucleus_volume',
        'Nucleus_high, micrometre': 'Nucleus_high',
        'AF594-T2 ring intensity coef': 'Ring_coefficient',
        'AF488-T3 av_signal_in_nuc_area_3D': 'Average_signal_488',
        'AF488-T3 total_signal_in_nuc_area_3D': 'Total_signal_488',
        'AF594-T2 av_signal_in_nuc_area_3D': 'Average_signal_594'
    }
    df = df.rename(columns=columns_to_rename)

    # Convert 'Time' from string to integer
    df['Time'] = df['Time'].str.replace('hr', '').astype(int)

    # Create 'Group' column
    df['Group'] = df['Cell Type'] + "_" + df['LIV'] + "_" + df['Cisp']+ "_" + df['Time'].astype(str) + "hr"

    return df



def plot_metrics_by_group_and_ring(type_to_analyse, lsm_df, RING_COEFF_CUT_OFF, color1="red", color2="blue"):
    # Filter DataFrame based on 'Cell Type'
    filtered_df = lsm_df[lsm_df['Cell Type'] == type_to_analyse].copy()

    # Ensure 'Time' is in the correct format for sorting
    filtered_df['Time'] = filtered_df['Time'].astype(int)

    # Sort 'Cisp' and 'LIV' as specified, then by 'Time'
    # Assuming 'Cisp' values correctly sort with "-Cisplatin" first, adjust if needed
    filtered_df.sort_values(by=['Cisp', 'LIV', 'Time'], inplace=True)

    # Create 'Group' column after sorting
    filtered_df['Group'] = filtered_df['Cisp'] + "_" + filtered_df['LIV'] + "_" + filtered_df['Time'].astype(str) + "hr"

    # Setup for plots, one per metric, all in a single column
    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12, 15), sharex=True)

    # Define the metrics to plot
    metrics = ['Nucleus_volume', 'Average_signal_488', 'Total_signal_488']
    titles = [
        f'Ring coeff {RING_COEFF_CUT_OFF} Cell Type: {type_to_analyse} \nNucleus Volume by Group and Ring',
        'Average Signal 488 by Group and Ring',
        'Total Signal 488 by Group and Ring'
    ]

    # Custom palette for hue based on 'Ring'
    palette = {False: color1, True: color2}

    # Create an explicit order for the x-axis based on the sorted 'Group'
    group_order = filtered_df['Group'].unique()

    for ax_idx, (ax, metric, title) in enumerate(zip(axs, metrics, titles)):
        sns.boxplot(x='Group', y=metric, hue='Ring', data=filtered_df, ax=ax, palette=palette, order=group_order)
        ax.set_title(title)
        ax.set_xlabel('Group' if metric == metrics[-1] else '')  # Only label the x-axis for the bottom plot
        ax.set_ylabel(metric)

        # For the last plot only, annotate each group with counts for False and True separately
        if metric == 'Total_signal_488':
            for i, group in enumerate(group_order):
                # Calculate counts for False and True within this group
                false_count = filtered_df[(filtered_df['Group'] == group) & (filtered_df['Ring'] == False)].shape[0]
                true_count = filtered_df[(filtered_df['Group'] == group) & (filtered_df['Ring'] == True)].shape[0]

                # Position the annotations below the group names on the x-axis
                # Adjust the positioning as necessary
                ax.text(i - 0.2, -0.1, f'no-ring:{false_count}', horizontalalignment='center', size='small', color='black', weight='semibold', transform=ax.get_xaxis_transform())
                ax.text(i + 0.2, -0.1, f'ring: {true_count}', horizontalalignment='center', size='small', color='black', weight='semibold', transform=ax.get_xaxis_transform())

    plt.tight_layout()
    return fig








