import json
from datetime import datetime

# Load data from JSON file
with open('data.json', 'r') as file:
    data = json.load(file)

# Get current date
current_date = datetime.now().strftime("%Y-%m-%d")

# LaTeX report template
latex_template = f"""
\\documentclass{{article}}
\\usepackage[utf8]{{inputenc}}

\\title{{{data['title']}}}
\\author{{{data['name']}}}
\\date{{{current_date}}}

\\begin{{document}}

\\maketitle

\\section{{Section 1}}
{data['section1']}

\\section{{Section 2}}
{data['section2']} 

\\section{{Section 3}}
{data['section3']} 

\\section{{Section 4}}
{data['section4']} 

\\section{{Section 5}}
{data['section5']} 

\\section{{Section 6}}
{data['section6']} 

\\section{{Section 7}}
{data['section7']} 

\\section{{Section 8}}
{data['section8']} 

\\end{{document}}
"""

# Write LaTeX report to file
with open('report.tex', 'w') as file:
    file.write(latex_template)

print("LaTeX report generated successfully.")