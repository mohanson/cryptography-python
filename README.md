# Cryptography Python

This is a tutorial project about cryptography in python. Follow the steps and you will fully understand the magic of elliptic curves from scratch.

- Each file can be run independently!
- There are detailed paper references in the comments of the code.

Based on this project, I implemented the SDKs of four famous blockchain projects: BTC, ETH, SOL and CKB.

- <https://github.com/mohanson/pabtc>: Almost fully supported.
- <https://github.com/mohanson/pleth>: Almost fully supported.
- <https://github.com/mohanson/pxsol>: Almost fully supported.
- <https://github.com/mohanson/pyckb>: Almost fully supported.

## Chapter 1: SECP256k1

SECP256k1 is the elliptic curve used by Bitcoin, Ethereum and many other cryptocurrencies, It's the most popular elliptic curve in the world. You can learn this algorithm through the original paper, but I recommend the following PDF tutorial instead.

<https://www.cs.miami.edu/home/burt/learning/Csc609.142/ecdsa-cert.pdf>

Try to run and read the following files one by one.

- `secp256k1.py`: implement finite field and SECP256k1 curve.
- `secp256k1_generate_public_key.py`: calculate bitcoin public key from private key.
- `secp256k1_sign.py`: sign messages and verify it.
- `secp256k1_extract_private_key.py`: extract the private key by two signatures that use the same k. This is the most serious and well-known attack method used by ECDSA.
- `secp256k1_jacobian`: jacobian projective space. This is a method used to speed up SECP256k1 and ECDSA.
- `secp256k1_schnorr.py`: bitcoin taproot upgrade (2021) uses the schnorr signature algorithm, which is also based on SECP256k1.

## Chapter 2: SECP256r1

Similar to SECP256k1, but with different parameters. In short, some scholars believe that SECP256k1 may contain a backdoor, so they modified some parameters of SECP256k1 and renamed it as SECP256r1. However, it is worth noting that:

- Currently, we have no evidence to prove that secp256k1 contains a backdoor.
- Currently, we have no evidence to prove that secp256r1 does not contain a backdoor.

Try to run and read the following files one by one.

- `secp256r1.py`: implement finite field and SECP256r1 curve.
- `secp256r1_recovery_pubkey.py`: recover pubkey from signature.

## Chapter 3: Ed25519

Implements the Ed25519 elliptic curve and eddsa signature algorithms. Because SECP256k1 has some security issues: it is vulnerable to attacks, so Edwards invented the ed25519 elliptic curve algorithm and the eddsa signature algorithm. This algorithm is currently used by Solana.

Read this article to understand why the SECP256k1 is not "perfect": <https://github.com/mohanson/pxsol/blob/master/doc/en/markdown/content/prikey_crypto_issue.md>

Try to run and read the following files one by one.

- [ed25519.py](https://github.com/mohanson/pxsol/blob/master/pxsol/ed25519.py)
- [eddsa.py](https://github.com/mohanson/pxsol/blob/master/pxsol/eddsa.py)

## Chapter 4: Polynomial arith

Polynomial manipulation in Python. In order to perform more complex calculations on elliptic curves, such as zero-knowledge proofs, we need to learn some basic knowledge first.

Try to run and read the following files one by one.

- `polynomial_numpy.py`: polynomial operation by numpy.
- `polynomial.py`: polynomial operation by hand writting python.

## Chapter 5: BN128

The zk algorithm used on Ethereum.

Try to run and read the following files one by one.

- `bn128.py`: bn128 curve, implements [eip-196](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-196.md), [eip-197](https://github.com/ethereum/EIPs/blob/master/EIPS/eip-197.md).
- `bn128_ethereum_api.py`: bn128 pairing testsuite.

## Chapter 6: KZG

KZG polynomial commitment is a scheme for polynomial commitment proposed by **K**ate, **Z**averucha, **G**oldberg and others in the paper Polynomial Commitments published in 2010.

Try to run and read the following files one by one.

- `kzg.py`: implement kzg proof based on bn128.
- `kzg_arith_relation.py`: a more complex proof.

It's difficult to explain the principles behind this algorithm; you'll need to read and study extensively on your own. This is the limit of my understanding of cryptography, therefore this tutorial ends here.

# Licences

MIT.
