import numpy as np
from numpy.polynomial import Polynomial

from utils import *
from Ciphertext import Ciphertext
from PublicKey import PublicKey

class CKKSEncrypter:
  def __init__(self, poly_degree, cipher_modulus, big_modulus, crt_context, public_key: PublicKey):
    self.poly_degree = poly_degree
    self.cipher_modulus = cipher_modulus
    self.big_modulus = big_modulus
    self.crt_context = crt_context
    self.public_key = public_key
    # self.secret_key = secret_key

  def encrypt(self, plaintext):
    p0 = self.public_key.p0
    p1 = self.public_key.p1

    random_poly = Polynomial(sample_from_triangle(self.poly_degree))
    err0 = Polynomial(sample_from_triangle(self.poly_degree))
    err1 = Polynomial(sample_from_triangle(self.poly_degree))

    poly_modulus = Polynomial([1] + [0] * (self.poly_degree-2) + [1])
    c0 = p0 * random_poly % poly_modulus
    c0 += self.cipher_modulus * err0
    c0 += plaintext
    # c0 %= self.big_modulus

    c1 = p1 * random_poly % poly_modulus
    c1 += self.cipher_modulus * err1
    # c1 %= self.big_modulus

    return Ciphertext(c0, c1)





