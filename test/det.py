import numpy as np
import sys

fp = open(sys.argv[1], "r")
nn = int(fp.readline().strip())
ham = np.fromfile(fp, dtype=float, count=nn*nn, sep=" ")
ham.shape = ( nn, nn )
fp.close()

det = np.linalg.det(ham)

print("determinant: %e" % (det))
