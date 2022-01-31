from email import message
import numpy as np
from numpy.polynomial import Polynomial
from Ciphertext import Ciphertext

class CKKSDecrypter:
  def __init__(self, poly_degree, secret_key, big_modulus, cipher_modulus):
    self.poly_degree = poly_degree
    self.big_modulus = big_modulus
    self.cipher_modulus = cipher_modulus
    self.secret_key = secret_key

  def decrypt(self, ciphertext):
    c0 = ciphertext.c0
    c1 = ciphertext.c1

    poly_modulus = Polynomial([1] + [0] * (self.poly_degree-2) + [1])
    sc1 = c1 * self.secret_key % poly_modulus
    message = c0 - sc1
    message %= self.big_modulus
    message %= self.cipher_modulus

    return message
