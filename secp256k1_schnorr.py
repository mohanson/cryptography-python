import hashlib
import random
import secp256k1

# https://github.com/bitcoin/bips/blob/master/bip-0340.mediawiki
#
# Schnorr signature variant Elliptic Curve Schnorr signatures for message m and public key P generally involve a point
# R, integers e and s picked by the signer, and the base point G which satisfy e = hash(R || m) and s⋅G = R + e⋅P. Two
# formulations exist, depending on whether the signer reveals e or R:
#   1. ....
#   2. Signatures are pairs (R, s) that satisfy s⋅G = R + hash(R || m)⋅P. This supports batch verification, as there
#      are no elliptic curve operations inside the hashes. Batch verification enables significant speedups.

prikey = secp256k1.Fr(0x5f6717883bef25f45a129c11fcac1567d74bda5a9ad4cbffc8203c0da2a1473c)
pubkey = secp256k1.G * prikey

# Hash of messages.
with open('./secp256k1.py', 'rb') as f:
    m = int.from_bytes(hashlib.sha256(f.read()).digest(), 'little')
    m = secp256k1.Fr(m)
print(f'hash={m}')

# R = k ∗ G
# e = hash(R || m)
# s = k + e ∗ prikey
k = secp256k1.Fr(random.randint(0, secp256k1.N))
R = secp256k1.G * k
hasher = hashlib.sha256()
hasher.update(R.x.x.to_bytes(32, 'little'))
hasher.update(R.y.x.to_bytes(32, 'little'))
hasher.update(m.x.to_bytes(32, 'little'))
e = secp256k1.Fr(int.from_bytes(hasher.digest(), 'little'))
s = k + e * prikey
print(f'sign=(R={R}, s={s})')

# s ∗ G =? R + hash(R || m) ∗ P
verify = secp256k1.G * s == R + pubkey * e
print(f'verify={verify}')
