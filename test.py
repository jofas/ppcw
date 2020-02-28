import os
import numpy as np


from util import read_data


WARN_LEVEL = 0.001
ERROR_LEVEL = 0.05


def main():
  for file in os.listdir("out/"):
    baseline = "data/%s" % file
    outfile = "out/%s" % file

    print("comparing %s:" % file)


    d1, d2 = read_data(outfile), read_data(baseline)

    diff = np.fabs(d1 - d2)
    sum = np.fabs(d1 + d2)

    np.where(sum != 0.0, sum, diff / sum)

    mx = max(np.sum(diff, axis=1))


    print("max = %s" % mx)

    assert not np.any(np.isnan(d1))
    assert not np.any(np.isnan(d2))
    assert mx <= ERROR_LEVEL

    if mx > WARN_LEVEL: print("WARNING!")


if __name__ == "__main__":
  main()
