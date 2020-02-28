import numpy as np

def read_data(file):
  with open(file) as f:
    return np.array([row.replace("  ", " ").split(" ")[1:]
      for row in f.read().split("\n")[:-1]], dtype=np.float64)
