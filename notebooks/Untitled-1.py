import pandas as pd
diamonds = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/diamonds.csv')
print(diamonds.cut.isin(["Ideal","Good","Very good"]).head())
print(diamonds.head())
print(diamonds[diamonds.cut.isin(["Ideal","Good","Very good"])].head())


# 
