#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
 
if __name__ == '__main__':
    import sys
    sourceStr = sys.argv[1]

    inp = open(sourceStr, "r")
    for (index, line) in enumerate(inp):
        k = line.split(" ")
        args = [float(t) for t in k[:-1]]
        imp = int(k[-1])
        colors = ['b' for t in args]
        colors[imp] = 'r'
        labels = list(range(len(args)))

        fig = plt.figure()
        width = 0.5
        ind = np.arange(len(args))
        plt.bar(ind, args, color=colors)
        plt.xticks(ind + width / 2, labels)
        plt.savefig("%s.pdf" % str(index))
