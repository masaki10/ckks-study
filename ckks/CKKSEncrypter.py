class CKKSEncrypter:
  def __init__(self, poly_degree, cipher_modulus, big_modulus, crt_context, public_key):
    self.poly_degree = poly_degree
    self.coeff_modulus = cipher_modulus
    self.big_modulus = big_modulus
    self.crt_context = crt_context
    self.public_key = public_key
    # self.secret_key = secret_key

  def encrypt(self, plaintext):
    pass

