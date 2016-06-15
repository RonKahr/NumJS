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
