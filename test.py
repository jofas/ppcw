import os
import numpy as np

WARN_LEVEL = 0.001
ERROR_LEVEL = 0.05


def main():
  for file1 in os.listdir("test/"):
    file2 = "test/%s" % file1

    print("comparing %s with %s:" % (file1, file2))

    d1, d2 = read_data(file1), read_data(file2)

    assert not np.any(np.isnan(d1))
    assert not np.any(np.isnan(d2))

    diff = np.fabs(d1 - d2)
    sum = np.fabs(d1 + d2)

    np.where(sum != 0.0, sum, diff / sum)

    r = np.sum(diff, axis=1)

    max = 0.0
    for x in r:
      if x > WARN_LEVEL:
        print("warning!")
      if x > max:
        max = x

    if max > ERROR_LEVEL:
      print("FAILURE")
    else:
      print("SUCCESS")

    print("max = %s" % max)


def read_data(file):
  with open(file) as f:
    return np.array([to_float(row[2:].replace("  ", " ").split(" "))
      for row in f.read().split("\n")[:-1]])


def to_float(vec_of_str):
  return [float(x) for x in vec_of_str]


if __name__ == "__main__":
  main()
