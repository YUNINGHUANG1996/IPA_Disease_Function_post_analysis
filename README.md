# IPA_Disease_Function_post_analysis
This is the python script to conduct basic and simple analysis over IPA output.
All the input files we need are:
1. IPA Output: Export the IPA Disease and Function results to Excel, remove the first row, and save it as a CSV file. Ensure the "B-H p-value" column (representing the adjusted p-value) is included in the output.
2. Reference Pathway File: This file links IPA’s main pathway categories (Molecular and Cellular Functions, Physiological System Development and Function, Diseases and Disorders) to their respective sub-pathways. A curated pathway reference file is also available for use.
3. Differential Expression Gene (DEG) Output: Include the output file from your differential expression analysis, detailing gene expression levels across conditions.

Running the Script
To use this script with your own data:

1. Update the file paths in ipa_analysis.py to point to your input files.
2. Ensure your Python environment includes the following libraries:
   import pandas as pd
   import matplotlib.pyplot as plt
   import numpy as np
   import seaborn as sns
   import os
   from openpyxl import Workbook
3. Execute the script by running:
   python ipa_analysis.py

