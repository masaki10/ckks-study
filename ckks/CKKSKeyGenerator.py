import random
import numpy as np
from numpy.polynomial import Polynomial

def sample_hamming_weight_vector(num_samples):
  sample = np.zeros(num_samples)
  total_weight = 0

  while total_weight < num_samples:
    index = random.randrange(0, num_samples)
    if sample[index] == 0:
      r = random.randrange(0, 2)
      if r == 0: sample[index] = -1
      else: sample[index] = 1
      total_weight += 1

  return sample

def sample_from_uniform_distribution(min_val, max_val, num_samples):
  assert(num_samples > 0 & isinstance(num_samples, int))

  return np.array([random.randrange(min_val, max_val) for _ in range(num_samples)])

def sample_from_triangle(num_samples):
  sample = np.zeros(num_samples)

  for i in range(num_samples):
    r = random.randrange(0, 4)
    if r == 0: sample[i] = -1
    elif r == 1: sample[i] = 1

  return sample


class CKKSKeyGenerator:
  def __init__(self, poly_degree, big_modulus):
    self.poly_degree = poly_degree
    self.big_modulus = big_modulus
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
    p0 += pk_err
    p0 %= self.big_modulus
    p1 = pk_coeff
    self.public_key = [p0, p1]