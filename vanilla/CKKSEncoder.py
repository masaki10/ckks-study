import numpy as np
from numpy.polynomial import Polynomial

class CKKSEncoder:
    """Basic CKKS encoder to encode complex vectors into polynomials."""
    
    def __init__(self, M: int):
        """Initialization of the encoder for M a power of 2. 
        
        xi, which is an M-th root of unity will, be used as a basis for our computations.
        """
        self.xi = np.exp(2 * np.pi * 1j / M)
        self.M = M
        
    @staticmethod
    def vandermonde(xi: np.complex128, M: int) -> np.array:
        """Computes the Vandermonde matrix from a m-th root of unity."""
        
        N = M //2
        matrix = []
        # We will generate each row of the matrix
        for i in range(N):
            # For each row we select a different root
            root = xi ** (2 * i + 1)
            row = []

            # Then we store its powers
            for j in range(N):
                row.append(root ** j)
            matrix.append(row)
        return matrix
    
    def sigma_inverse(self, b: np.array) -> Polynomial:
        """Encodes the vector b in a polynomial using an M-th root of unity."""

        # First we create the Vandermonde matrix
        A = CKKSEncoder.vandermonde(self.xi, self.M)

        # Then we solve the system (連立方程式)
        coeffs = np.linalg.solve(A, b)
        
        # Finally we output the polynomial
        p = Polynomial(coeffs)
        
        return p

    def sigma(self, p: Polynomial) -> np.array:
        """Decodes a polynomial by applying it to the M-th roots of unity."""

        outputs = []
        N = self.M //2

        # We simply apply the polynomial on the roots
        for i in range(N):
            root = self.xi ** (2 * i + 1)
            output = p(root)
            outputs.append(output)
        return np.array(outputs)

if __name__ == '__main__':
  # First we set the parameters
  M = 8
  N = M //2

  # We set xi, which will be used in our computations
  xi = np.exp(2 * np.pi * 1j / M)
#   print(xi)
  # (0.7071067811865476+0.7071067811865475j)

  # First we initialize our encoder
  encoder = CKKSEncoder(M)

  b = np.array([1, 2, 3, 4])
  
  # encode b
  p = encoder.sigma_inverse(b)
#   print(p)
  # (2.5+8.881784197001252e-16j) -
  # (2.3551386880256634e-16-0.7071067811865478j)·x¹ -
  # (5.13478148889135e-16-0.5000000000000002j)·x² -
  # (8.635508522760774e-16-0.7071067811865474j)·x³

  # decode b
  b_reconstructed = encoder.sigma(p)
#   print(b_reconstructed)
  # [1.+1.11022302e-16j 2.+5.82867088e-16j 3.+2.49800181e-16j
  # 4.+4.44089210e-16j]

  err = np.linalg.norm(b_reconstructed - b)
#   print(err)
  # 7.820967653285533e-16

  m1 = np.array([1, 2, 3, 4])
  m2 = np.array([1, -2, 3, -4])
  m_add = m1 + m2

  p1 = encoder.sigma_inverse(m1)
#   print(p1)
  # (2.5+8.881784197001252e-16j) -
  # (2.3551386880256634e-16-0.7071067811865478j)·x¹ -
  # (5.13478148889135e-16-0.5000000000000002j)·x² -
  # (8.635508522760774e-16-0.7071067811865474j)·x³

  p2 = encoder.sigma_inverse(m2)
#   print(p2)
  # (-0.49999999999999956-1.1102230246251565e-16j) -
  # (0.7071067811865471+1.0205600981444541e-15j)·x¹ +
  # (2.3453461395206428e-15-2.5j)·x² +
  # (0.7071067811865476+1.8841109504205327e-15j)·x³

  p_add = p1 + p2
#   print(p_add)
  # (2.0000000000000004+7.771561172376096e-16j) -
  # (0.7071067811865474-0.7071067811865468j)·x¹ +
  # (1.831867990631508e-15-1.9999999999999998j)·x² +
  # (0.7071067811865467+0.7071067811865492j)·x³

  decoded_p_add = encoder.sigma(p_add)
#   print(decoded_p_add)
  # [2.0000000e+00+3.06128380e-16j 0.0000000e+00+1.22124533e-15j
  # 6.0000000e+00+7.77156117e-16j 8.8817842e-16+4.44089210e-16j]

#   print(m_add - decoded_p_add)
  # [ 0.0000000e+00-3.06128380e-16j  0.0000000e+00-1.22124533e-15j
  # 0.0000000e+00-7.77156117e-16j -8.8817842e-16-4.44089210e-16j]

  # X^N + 1
  poly_modulo = Polynomial([1,0,0,0,1])
#   print(poly_modulo)
  # 1.0 + 0.0·x¹ + 0.0·x² + 0.0·x³ + 1.0·x⁴

  p_mult = p1 * p2 % poly_modulo
#   print(p_mult)
  # (-2.499999999999999-3.7816971776294395e-15j) -
  # (3.5355339059327346+0.7071067811865546j)·x¹ +
  # (1.1227130336521899e-14-7.5j)·x² +
  # (3.53553390593274-0.7071067811865394j)·x³

  decoded_p_mult = encoder.sigma(p_mult)
#   print(decoded_p_mult)
  # [  1.-4.51028104e-16j  -4.+1.04083409e-16j   9.+1.54737334e-15j
  # -16.-1.17753030e-14j]