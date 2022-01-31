import random
import numpy as np
from numpy.polynomial import Polynomial
from utils import *
from PublicKey import PublicKey

class CKKSKeyGenerator:
  def __init__(self, poly_degree, big_modulus, cipher_modulus):
    self.poly_degree = poly_degree
    self.big_modulus = big_modulus
    self.cipher_modulus = cipher_modulus
    self.generate_secret_key()
    self.generate_public_key()

  def generate_secret_key(self):
    key = sample_hamming_weight_vector(self.poly_degree)
    self.secret_key = Polynomial(key)

  def generate_public_key(self):
    coeff = sample_from_uniform_distribution(0, self.big_modulus, self.poly_degree)
    pk_coeff = Polynomial(coeff)
    err = sample_from_triangle(self.poly_degree)
    pk_err = Polynomial(err)
    poly_modulus = Polynomial([1] + [0] * (self.poly_degree-2) + [1])
    p0 = pk_coeff * self.secret_key % poly_modulus
    p0 %= self.big_modulus
    p0 *= -1
    p0 += pk_err * self.cipher_modulus # ここ微妙
    p0 %= self.big_modulus
    p1 = pk_coeff
    self.public_key = PublicKey(p0, p1)