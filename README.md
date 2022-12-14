## secp256k1

paper: <https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf>

- `secp256k1.py`: implement finite field and secp256k1 curve.
- `secp256k1_generate_public_key.py`: calculate bitcoin public key from private key.
- `secp256k1_sign.py`: signature messages and verify it.
- `secp256k1_extract_private_key.py`: extract the private key by two signatures that use the same k.
- `secp256k1_jacobian`: jacobian projective space.

## polynomial arith

- `polynomial_numpy.py`: polynomial operation by numpy.
- `polynomial.py`: polynomial operation by hand writting python.

## bn128

- `bn128.py`: bn128 curve, implements [eip-196](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-196.md), [eip-197](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-197.md).
- `bn128_ethereum_api.py`: bn128 pairing testsuite.
