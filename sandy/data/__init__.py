import os
import pandas as pd

elements = pd.read_csv(os.path.join(__path__[0], "elements.csv")).set_index("Z")
metastates = pd.read_csv(os.path.join(__path__[0], "metastates.csv")).set_index("LISO")