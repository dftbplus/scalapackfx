import numpy as np
import sys

print(sys.argv[1])

fp = open(sys.argv[1], "r")
nn = int(fp.readline().strip())
matrix = np.fromfile(fp, dtype=float, count=nn*nn, sep=" ")
matrix.shape = ( nn, nn )
fp.close()

U, sigmas, Vt = np.linalg.svd(matrix)

print("SIGMA:")
for sigma in sigmas:
    print("{:23.15E}".format(sigma))
