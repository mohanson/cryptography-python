# 多项式承诺

承诺(commitment)是许多密码协议的基本组成部分. 一个安全的承诺⽅案允许提交者发布⼀个称为承诺的值, 该值绑定到⼀条消
息上⽽不泄露该消息. 稍后, 它可能会打开承诺并将承诺的消息透露给验证者, 验证者可以检查消息是否与承诺⼀致.

多项式承诺(polynomial commitment)利用一个嵌入度为 d 的多项式隐藏所谓的秘密, 这个秘密可能是多项式的系数, 也可能就是 d + 1 个点. 当多项式的嵌入度远远小于验证域的阶时候, 密码学多年的研究已经证明, 可以通过只验证点值的方式就能实现计算意义上的完备性(computational soundness).

显而易见, 如果只验证少数几个点值对, 计算工作量远远小于验证 d + 1 个点值对插值求多项式, 这也提供了良好的简洁性, 即 succinct. 但是简单的利用多项式来实现"通用", "零知识"等重要的密码学特征, 还需要做很多工作.

## Motivation

0. Given a list of data yi, i=0, ...n-1, commitment is a "compressed" digest of the data
0. Merkle tree commitment = root hash(32 bytes)
0. Prover: given i, generate proof pi, yi
0. Verifier: given pi, yi, return true of false
0. Cons: log(n) proof size

KZG

0. commitment: 48 bytes
0. Single proof size: 48 bytes(no matter n is)
0. multiple proof size: 48 bytes

Applications:

0. Stateless(Verkle tree) validator 不需要保存全部账本对区块数据进行验证. 把所有交易触碰到账本的数据证明这些数据是真的
0. DAS

## Polynomial

f(x) = a₀ + a₁x + ... + aₙxⁿ

- Degree(嵌入度) deg(f(x)) = n
- aₙ != 0

## Encoding data into polynomial using Lagrange Interpolation

Given (xi, yi), xi != xj, 构造一个多项式 f(xi) = yi degree is n - 1

<https://zh.wikipedia.org/wiki/%E6%8B%89%E6%A0%BC%E6%9C%97%E6%97%A5%E6%8F%92%E5%80%BC%E6%B3%95>

```py
import scipy

x = [1,  2,  3,  4]
y = [4, 15, 40, 85]
ret = scipy.interpolate.lagrange(x, y)

print(ret) # x³ + x² + x + 1
```

```py
def h(x,y,a):
    ans=0.0
    for i in range(len(y)):
        t=y[i]
        for j in range(len(y)):
            if i !=j:
                t*=(a-x[j])/(x[i]-x[j])
        ans +=t
    return ans
x=[1,0]
y=[0,2]
print(h(x,y,2))
```

任意 n 个 point 都能唯一确定 degree is n - 1 的多项式. 例如: 2 点确定一条直线.

## Polynomial on an elliptic curve G1

G: generate point

f(x) * G = ai(xi * G) 先算 G1 和 xi 次方的结果

## Polynomial commitment with trusted setup

Now we have secret s belong to Fq, such that

- Nobody knows s (private key of the "god")
- [si]G1, i=1, ... is known to every body("god" public key)

Then we have the commitment as

C = [f(s)]1 = sum(ai * si * G) 不知道 s 和 f(s), 但知道在椭圆曲线上的值

finding another g(x) such that g(s) = f(s) is almost impossible.

ai 是原信息.

## Single proof

Given xi, yi, want to prove f(xi) = yi.

f(x) - yi = q(x)(x-xi)

-> [f(s) - yi] * G1 = [(s-xi)q(s)]G1
-> C - [yi]*G1 =

e: G1 x G2 -> GT

e(C - [yi]*G1, G2) = e([q(s)]xG1, [(s - xi)] x G2)

where [q(s)]*G1 is the proof (48 bytes as a point on an elliptic curve)

## 参考

https://www.youtube.com/watch?v=n4eiiCDhTes
https://www.youtube.com/watch?v=NVvNHe_RGZ8
https://copper-witch-857.notion.site/Polynomial-KZG-or-Kate-Commitment-DappLearning-Notes-fc426c8cb9a14878840852506865f13b
