# NumJS
A JavaScript library for numerical and scientific calculations.

## Number Theory
###### mod(n,m)
Modulus operator. There is a bug in the built-in modulus operator “%” for javascript where it doesn’t   return the correct value for negative numbers.  This function will return the correct value for both positive and negative numbers.  

Usage  
numjs.mod(-74, 31) //will result in 19  
numjs.mod(101, 31) //will result in 8  

###### isPrime(n)
Returns true if the number supplied is prime and false if it is not.

Usage  
numjs.isPrime(13) //returns true  
numjs.isPrime(6857) //returns true  
numjs.isPrime(15) //returns false  

###### gcd(a,b)
Returns the largest divisor common to a and b.

Usage  
numjs.gcd(20, 30)  //returns 10  
numjs.gcd(6857, 7919) //return 1  

###### extendedEuclidean(a,b)
For two integers, a and b, returns the greatest common divisor as well as the integers x and y such that  
ax + by = gcd(a,b).

Usage  
numjs.extendedEuclidean(89, 43)  //returns {x:-14,y:29,gcd:1}   
numjs.extendedEuclidean(99, 78)  //returns {x:-11,y:14,gcd:3}  
numjs.extendedEuclidean(240, 46) //returns {x:-9,y:47,gcd:2}  

###### nthPrime(n)
Return the nth prime.  

Usage  
numjs.nthPrime(1) //returns 2  
numjs.nthPrime(168) //returns 997  

###### primesToN(n)
Returns a list of primes less than or equal to n.

Usage  
numjs.primesToN(100) //returns [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97];  

###### primeFactors(n)
Returns the factorization of a number into its constituent primes.

Usage  
numjs.primeFactors(984) //returns [{ "power": 3, "prime": 2 }, { "power": 1, "prime": 3 }, { "power": 1, "prime": 41 }]  
numjs.primeFactors(499) //returns [{ "power": 1, "prime": 499 }]  

###### factorCount(n)
Returns the number of factors(divisors) of n.  For example, 28 has six divisors: 1,28,2,14,4,7.

Usage  
numjs.factorCount(28)  //returns 6  

## Combinatorics

## binomialCoefficient(n,k)
Returns the The binomial coefficient of n,k.  The binomial coefficient represents the number of ways of picking k unordered outcomes from n possibilities, also known as a combination or combinatorial number.

Usage  
numjs.binomialCoefficient(4, 2) //returns 6  
numjs.binomialCoefficient(40, 20) //returns 137846528820  

## Linear Algebra
#### RVector(n)  
Initialize a vector of real numbers of size n.  

##### Properties
dimension:  Returns the dimension of the vector.  
DELTA:  Returns the delta or error value for comparing the vector to another vector using isDeltaEqual  

##### Methods
add(B): Adds vector B (of size n or less) to the current vector.  
dotProduct(B):  Finds the dot product of the current vector and another vector B.  
get(i):  Returns the element and index i.  
set(i,value):  Sets the value of the element at index i to value.  
innerProductNorm(W): Returns the inner product norm given a matrix W  
matrixMultiply(W):  Multiplies the current vector by a matrix of real numbers (W) and return the resulting vector  
subtract(B):  Subtracts a vector B from the current vector.  
normalize():  Normalizes the vector (with respect to the Euclidean norm).  
norm():  Returns the Euclidean (2-norm) of the vector.  
norm1():  Returns the 1-norm of the vector.  
normP():  Returns the p-norm of the vector.  
normInf():  Returns the infinity-norm.  
scalarMultiply(x):  Multiply the vector by a scalar value x.  
isStrictEqual(v2):  Compare the values of the vector to another vector (v2).  If any of the values are not equal("==") return false.  
isDeltaEqual(v2):  Compare the values of the vector to another vector (v2).  If the absolute value of the difference is greater than DELTA return false.  
toString(): A string representation of the vector.  
copy(): Returns a copy of the current vector.  

#### RMatrix(n,m)
Initialize a matrix of real numbers of n rows by m columns

##### Methods
add(B):  Adds matrix B to the current matrix.  
copy():  Returns a copy of the matrix.  
getColumnVector(col):  Return a RVector object of the specified column.  
getRowVector(row):  Return a RVector object of the specified row.  
multiply(B):  Return a RMatrix object of the AxB product.  
subtract(B):  Subtract the RMatrix B from the current matrix, A-B.  
transpose():  Return the transpose (A<sup>T</sup>) of the current matrix.
trace():  Return the trace of the current matrix.  

