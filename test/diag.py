import numpy as np
import sys

fp = open(sys.argv[1], "r")
nn = int(fp.readline().strip())
ham = np.fromfile(fp, dtype=float, count=nn*nn, sep=" ")
ham.shape = ( nn, nn )
fp.close()

fp = open(sys.argv[2], "r")
nn = int(fp.readline().strip())
over = np.fromfile(fp, dtype=float, count=nn*nn, sep=" ")
over.shape = ( nn, nn )
fp.close()

print("Cholesky")

print("TRANS:")
ll = np.linalg.cholesky(over)
htrans =  np.dot(np.dot(np.linalg.inv(ll), ham),
                 np.linalg.inv(np.transpose(ll)))
eigvals, eigvecs_trans = np.linalg.eig(htrans)
eigvecs = np.linalg.solve(np.transpose(ll), eigvecs_trans)

print("EIGENVALUES:")
for eigval in eigvals:
    print("{:23.15E}".format(eigval))

print("EIGVECS0:")
formstr = "{:14.6E}" * eigvecs.shape[0]
for vec in eigvecs.transpose():
    print(formstr.format(*vec))
#
print("DIFF:")
print(np.dot(ham, eigvecs[:,0]) - eigvals[0] * np.dot(over, eigvecs[:,0]))
