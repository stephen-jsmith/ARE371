import pandas as pd
import tkinter as tk
from tkinter import filedialog

def load_excel_file(file_path=None) -> pd.DataFrame:
    """
    Loads an Excel file into a pandas dataframe.

    Parameters:
        file_path (str): The path to the Excel file.

    Returns:
        pandas.DataFrame: The loaded dataframe.
    """
    root = tk.Tk()
    root.withdraw()
    if not file_path:
        file_path = filedialog.askopenfilename(
            title="Select an Excel file",
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")],
        )

        # Check if a file was selected
        if not file_path:
            print("No file selected.")
            exit()
        elif not file_path.endswith(".xlsx"):
            print("Invalid file type. Please select an Excel file.")
            exit()

    try:
        df = pd.read_excel(file_path)
        return df
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred while loading the Excel file: {e}")
        return None


# Example usage
df = load_excel_file()
if df is not None:
    print(df.head())
