import numpy as np
from statistics import mean

from fire import Fire

from util import read_data

def main(bench1, bench2):
  b1, b2 = read_data(bench1), read_data(bench2)

  print("stats for %s:" % bench1)
  print_stat(b1)

  print("\nstats for %s:" % bench2)
  print_stat(b2)


  print("\ncomparison:")

  perc = 1.0 - (mean(b2[:,-1]) / mean(b1[:,-1]))

  if perc < 0.0:
    print("REGRESSION: %s is %f%% slower than %s" % (bench2, perc * 100, bench1))
  else:
    print("IMPROVEMENT: %s is %f%% faster than %s" % (bench2, perc * 100, bench1))


def print_stat(b):
  avg_time_100_iter = np.sum(b[:,:-1]) / float(b[:,:-1].size)

  not_io = np.sum(b[:,:-1],axis=1)
  overall = b[:,-1]
  io = mean(overall - not_io)

  print("avg time per 100 iterations: %f" % avg_time_100_iter)
  print("avg time spent doing io:     %f" % io)
  print("avg percentage io:           %f" % mean((overall - not_io) / overall))
  print("avg overall time:            %f" % mean(overall))


if __name__ == "__main__":
  Fire(main)
