import sys
sys.path.append("../ckks")
from CKKSEncoder import CKKSEncoder
import numpy as np

M = 8
scale = 64
encoder = CKKSEncoder(M, scale)
z = np.array([0,1])
p = encoder.encode(z)
print(f"encoded\n{p}")
print(f"decoded\n{encoder.decode(p)}")