# Cryptography Python

This is a tutorial project about cryptography in python. Follow the steps and you will fully understand the magic of elliptic curves from scratch.

- Each file can be run independently!
- There are detailed paper references in the comments of the code

Based on this project, I implemented the SDKs of three famous blockchain projects: BTC, ETH and CKB.

- <https://github.com/mohanson/pybtc>: Almost fully supported.
- <https://github.com/mohanson/pyeth>: Almost fully supported.
- <https://github.com/mohanson/pyckb>: Almost fully supported.

## Secp256k1

Written according to this paper: <https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf>

- `secp256k1.py`: implement finite field and secp256k1 curve.
- `secp256k1_generate_public_key.py`: calculate bitcoin public key from private key.
- `secp256k1_sign.py`: signature messages and verify it.
- `secp256k1_extract_private_key.py`: extract the private key by two signatures that use the same k.
- `secp256k1_jacobian`: jacobian projective space.

## Polynomial arith

Polynomial manipulation in Python.

- `polynomial_numpy.py`: polynomial operation by numpy.
- `polynomial.py`: polynomial operation by hand writting python.

## Bn128

- `bn128.py`: bn128 curve, implements [eip-196](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-196.md), [eip-197](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-197.md).
- `bn128_ethereum_api.py`: bn128 pairing testsuite.

## Kzg

KZG polynomial commitment is a scheme for polynomial commitment proposed by Kate, Zaverucha, Goldberg and others in the paper Polynomial Commitments published in 2010.

- `kzg.py`: implement kzg proof based on bn128.
- `kzg_arith_relation.py`: a more complex proof.

## Secp256r1

Similar to secp256k1, but with different parameters.

- `secp256r1.py`: implement finite field and secp256r1 curve.
- `secp256r1_recovery_pubkey.py`: recover pubkey from signature.

## Ed25519

Implements the Ed25519 elliptic curve and eddsa signature algorithms.

- [ed25519.py](https://github.com/mohanson/pxsol/blob/master/pxsol/ed25519.py)
- [eddsa.py](https://github.com/mohanson/pxsol/blob/master/pxsol/eddsa.py)

# Licences

MIT.
