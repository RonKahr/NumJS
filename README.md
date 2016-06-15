# NumJS
A JavaScript library for numerical and scientific calculations.

## Number Theory
######  mod(n,m)
	Modulus operator. There is a bug in the built-in modulus operator “%” for javascript where it doesn’t return the correct value for negative numbers.  This function will return the correct value for both positive and negative numbers.

  	Usage
  		numjs.mod(-74, 31) //will result in 19
  		numjs.mod(101, 31) //will result in 8


