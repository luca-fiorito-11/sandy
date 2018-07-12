import os
import pandas as pd

jaea_fns_175 = pd.read_csv(os.path.join(__path__[0], "JAEA_FNS_175.csv")).set_index("E")
