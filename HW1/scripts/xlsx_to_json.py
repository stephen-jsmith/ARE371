import json
from load_xlsx import load_excel_file

def convert_xlsx_to_json(xlsx_file, json_file):
    data = load_excel_file(xlsx_file)
    with open(json_file, 'w') as f:
        json.dump(data, f)


if __name__ == "__main__":
    # Example usage
    xlsx_file = '/path/to/input.xlsx'
    json_file = '/path/to/output.json'
    convert_xlsx_to_json(xlsx_file, json_file)
