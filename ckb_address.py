import hashlib

import segwit_addr as sa
import secp256k1

# Double checked by https://ckb.tools/generator

prikey = secp256k1.Fr(0xd5d8fe30c6ab6bfd2c6e0a940299a1e01a9ab6b8a8ed407a00b130e6a51435fc)
pubkey = secp256k1.G * prikey
assert pubkey.x.x == 0x97202631ccab00b8669e0b1fcc376f082513f22593c5e99fbf76ab02e8911d2e
assert pubkey.y.x == 0xeae37bf649d45e0cf83c5c057de60d685ece29e9b7e58959a638845d3d0659c6

if pubkey.y.x & 1 == 0:
    pubkey_byte = bytes([0x02]) + pubkey.x.x.to_bytes(32, byteorder='big')
else:
    pubkey_byte = bytes([0x03]) + pubkey.x.x.to_bytes(32, byteorder='big')
assert pubkey_byte.hex() == '0297202631ccab00b8669e0b1fcc376f082513f22593c5e99fbf76ab02e8911d2e'

args = hashlib.blake2b(pubkey_byte, digest_size=32, person=b'ckb-default-hash').digest()[:20]
assert args.hex() == 'e5126d9d897e5d5249607760f9da024119f9e296'

# https://github.com/nervosnetwork/rfcs/blob/master/rfcs/0021-ckb-address-format/0021-ckb-address-format.md
# https://github.com/rev-chaos/ckb-address-demo/blob/master/ckb_addr_test.py
payload = b''
payload += b'\x00'
# Append secp256k1 code hash
payload += bytes.fromhex('9bd7e06f3ecf4be0f2fcd2188b23f1b9fcc88e5d4b65a8637b17723bbda3cce8')
payload += b'\x01'
payload += args

cktaddr = sa.bech32_encode('ckt', sa.convertbits(payload, 8, 5), sa.Encoding.BECH32M)
assert cktaddr == 'ckt1qzda0cr08m85hc8jlnfp3zer7xulejywt49kt2rr0vthywaa50xwsq09zfkemzt7t4fyjcrhvrua5qjpr8u799s6se0vv'
