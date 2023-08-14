# poly_arithmetic
Polynomial arithmetic (Verilog)

## Polynomial operations

* Multiplication in GF(2)
* Division in GF(2)
* Modular multiplication in GF(2)
* Multiplication in GF(2^m)
* Division in GF(2^m)
* Greatest Common Divisor in GF(2)
* Least Common Multiple in GF(2)

## Polynomial representation

### in GF(2)
`if` p(x) = a[0] + a[1]*x + a[2]*x^2 + ... + a[k]*x^k (a[i] is in GF(2), k <= POLY_P_WIDTH-1)
`then` px[POLY_X_WIDTH-1:0] = {a[k], ..., a[2], a[1], a[0]}

### in GF(2^m)
`if` p(x) = b[0] + b[1]*x + b[2]*x^2 + ... + b[k]*x^k (b[i] is in GF(2^m), k <= POLY_P_WIDTH-1) and width of b[i] is M_VALUE, bit-order of b[i] is little-endian
`then` px[POLY_X_WIDTH-1:0] = {b[k], ..., b[2], b[1], b[0]}

### Note: deg(p(x)) = (POLY_P_WIDTH - 1)

## Multiplication in GF(2)
### Operation: c(x) = a(x) * b(x)
### Parameters
| parameter | |
|-|-|
|POLY_A_WIDTH|Width of polynomial A|
|POLY_B_WIDTH|Width of polynomial B|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pc|Polynomial C (product)|

## Division in GF(2)
### Operation: a(x) / b(x) -> a(x) = b(x) * q(x) + r(x)
### Parameters
| parameter | |
|-|-|
|POLY_A_WIDTH|Width of polynomial A|
|POLY_B_WIDTH|Width of polynomial B|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pq|Polynomial Q (quotient)|
|pr|Polynomial R (remainder)|

## Modular multiplication in GF(2)
### Operation: c(x) = a(x)*b(x) (mod p(x))
### Parameters
| parameter | |
|-|-|
|POLY_A_WIDTH|Width of polynomial A|
|POLY_B_WIDTH|Width of polynomial B|
|POLY_P_WIDTH|Width of Polynomial P|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pp|Polynomial P (must be irreducible over GF(2)|
|pc|Polynomial C (product)|

## Multiplication in GF(2^m)
### Operation: c(x) = a(x) * b(x)
### Parameters
| parameter | |
|-|-|
|POLY_A_WIDTH|Width of polynomial A|
|POLY_B_WIDTH|Width of polynomial B|
|M_VALUE|Degree of 2 (m)|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pp|Polynomial P (must be irreducible over GF(2)|
|pc|Polynomial C (product)|

## Division in GF(2^m)
### Operation: a(x) / b(x) -> a(x) = b(x) * q(x) + r(x)
### Requirements: b(x) must be type of b[0]+b[1]*x+b[2]*x^2+...+b[k-1]*x^(k-1)+1*x^k, where k <= (POLY_B_WIDTH-1) and b[i] is element of GF(2^m)
### Parameters
| parameter | |
|-|-|
|POLY_A_WIDTH|Width of polynomial A|
|POLY_B_WIDTH|Width of polynomial B|
|M_VALUE|Degree of 2 (m)|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pp|Polynomial P (must be irreducible over GF(2)|
|pq|Polynomial Q (quotient)|
|pr|Polynomial R (remainder)|

## Greatest Common Divisor in GF(2)
### Operation: c(x) = GCD(a(x), b(x))
### Parameters
| parameter | |
|-|-|
|POLY_WIDTH|Width of Polynomials|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pc|Polynomial C (product)|

## Least Common Multiple in GF(2)
### Operation: c(x) = LCM(a(x), b(x))
### Parameters
| parameter | |
|-|-|
|POLY_WIDTH|Width of Polynomials|
### Ports
| port | |
|-|-|
|pa|Polynomial A|
|pb|Polynomial B|
|pc|Polynomial C (product)|
