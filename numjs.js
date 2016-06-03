(function(window){
    'use strict';
    function define_numjs() {
        var numjs = {};

        //NumberTheory Functions 
        numjs.mod = function (n, m) {
            return ((n % m) + m) % m;
        }
        numjs.isPrime = function(n){
			if (n == 1) { return false; }
			if (n % 2 == 0 || n % 3 == 0) { return false; }
			for (var  i = 5; i * i <= n; i += 6) {
				if (n % i == 0 || n % (i + 2) == 0) { return false; }
			}
			return true;
        };
		numjs.gcd = function(a,b)
		{
		    var tmp = 0;
			if (b==0) {
				return a;
			}
			return numjs.gcd(b, a % b);
		};
		numjs.extendedEuclidean = function (a, b) {
		    var factors = {};
		    
		    if (b == 0) {
		        factors = { gcd: a, x: 1, y: 0 };
		        return factors;
		    }
		    else
		    {
		        var factors1 = numjs.extendedEuclidean(b, a % b);
		        factors = {gcd:factors1.gcd,x:factors1.y,y:factors1.x-Math.floor(a/b)*factors1.y};
		        return factors;
		    }
		};
		numjs.modularLES = function (a, b, n) {
		    var results = [];
		    var calc = 0;
		    var x0 = 0;
		    var factors = numjs.extendedEuclidean(a, n);
		    if (b % factors.gcd == 0)
		    {
		        x0 = numjs.mod((factors.x * (b / factors.gcd)), n);
		        for(var i=0;i<factors.gcd;i++)
		        {
		            calc = (x0 + i * (n / factors.gcd)) % n;
		            results.push(calc);
		        }
		    }
		    return results;

		};
		numjs.nthPrime = function (n) {
		    if (n == 1)
		        return 2;
		    var count = 2;
		    var candidate = 3;
		    while (count<n)
		    {
		        candidate += 2;
		        if (numjs.isPrime(candidate))
		            count++
		    }
		    return candidate;
		};
		numjs.primesToN = function(n)
        {
            var primes=[];
            var sievebound = Math.ceil(n/2);
            var sieve = new Array(sievebound), i = 0;
            while(i < sievebound)
            {		
                sieve[i++] = false;
            }
            var crosslimit=(Math.floor(Math.sqrt(n))-1)/2;
            for(var i=1;i<=crosslimit;i++)
            {
                if(!sieve[i])
                {
                    for(var j=2*i*(i+1);j<=sievebound;j+=2*i+1)
                    {
                        sieve[j] = true;
                    }
                }
            }
            primes.push(2);
            for(var k=1;k<sievebound;k++)
            {
                if(!sieve[k])
                {
                    primes.push(2*k+1);
                }
            }
            return primes;
		};
		numjs.primeFactors = function (n) {
		    var primeFactors = [];
		    var result = n;
		    var limit = Math.floor(Math.sqrt(n));
		    var primesToLimit = numjs.primesToN(limit);
		    var factorPower = 1;
		    for (var i = 0; i < primesToLimit.length; i++) {
		        if (result % primesToLimit[i] == 0)
		        {
		            factorPower = 0;
		            primeFactors.push({ prime: primesToLimit[i], power: factorPower })
		            while (result % primesToLimit[i] == 0)
		            {
		                result = result / primesToLimit[i];
		                factorPower++;

		            }
		            primeFactors[primeFactors.length - 1].power = factorPower;
		        }
		        if (result == 1)
		            break;
		    };
		    if (result > 1)
		        primeFactors.push({ prime: result, power: 1 });
		    return primeFactors;
		};
		numjs.eulersTotient = function (n) {
            //Implemented to avoid rounding errrors
		    var primeFactors = numjs.primeFactors(n);
		    var totient = n;
		    var numerator = 1;
		    var denomenator = 1;
		    for (var i = 0; i < primeFactors.length; i++) {
		        //This calculation uses the definintion of the totient function
		        //but may generate rounding errors.
                //Instead, multiply all the numerators and denomenators separately.
		        //totient *= (1-(1/primeFactors[i].prime))
		        numerator *= (primeFactors[i].prime - 1);
		        denomenator *= primeFactors[i].prime;
		    }
		    totient = (totient * numerator) / denomenator;
		    return totient;
		};
		numjs.factorCount = function(n) {
			var primeFactors = [];
			var numfactors = 1;
			primeFactors = numjs.primeFactors(n);
			for (var i = 0; i < primeFactors.length; i++)
			{
			    numfactors *= (primeFactors[i].power + 1);
			}
			return numfactors;
		};

        //Combinatorics
		numjs.binomialCoefficient = function(n,k)
		{	
			var coef = 1;
			
			for(var i=1;i<=k;i++)
			{
				coef *= (n+1-i)/i;
			}
			return coef;
		}
        
        //Linear Algebra
        //RVector
		numjs.RVector = function (n) {
			if (isNaN(n))
				throw new Error("Vector dimension must be a number and greater than 0.");
			var N = n;
			var vector = [];
			for(var i=0;i<n;i++)
			{
				vector[i]=0;
			}

			this.DELTA = 1e-15;
			this.dimension = function () {
			    return vector.length;
			};

			this.get = function(i) {
			    return vector[i];
			};

			this.set = function (i, value) {
			    vector[i] = value;
			};

		};
		numjs.RVector.prototype.add = function(B) {
			if (B.dimension()>this.dimension())
				throw new Error("Vector to be added cannot have a dimension greater than the current vector.");
			var value=0;
			for (var i=0; i< this.dimension(); i++)
			{
				value = this.get(i) + B.get(i);
				this.set(i,value);
			}
		};
		numjs.RVector.prototype.dotProduct = function(B) {
			if (B.dimension()!=this.dimension())
				throw new Error("Vector dimensions do not match.");
	
			var C = 0;
		
			for (var i=0; i<this.dimension();i++){
				C+=this.get(i)*B.get(i);
			}
		
			return C;
		};
		numjs.RVector.prototype.subtract = function(B) {
			if (B.dimension()!=this.dimension())
				throw new Error("Vector to be subtract must be of equal dimension to the current vector.");
		
			var value=0;
		
			for (var i=0; i< this.dimension(); i++)
			{
				value = this.get(i) - B.get(i);
				this.set(i,value);
			}
		};
		numjs.RVector.prototype.scalarMultiply = function(scalar) {
			var value=0;
		
			for (var i=0; i< this.dimension(); i++)
			{
				value = scalar * this.get(i);
				this.set(i,value);
			}
		};
		//Return the Euclidean (2-norm) of the vector
		numjs.RVector.prototype.norm = function()
		{
			var norm = 0;
			
			for (var i=0; i<this.dimension();i++)
			{
				norm += this.get(i)*this.get(i);
			}
			
			norm = Math.sqrt(norm);
			return norm;
		};
		numjs.RVector.prototype.norm1 = function()
		{
			var sum = 0;
			
			for(var i=0;i<this.dimension();i++){
				sum += Math.abs(this.get(i));
			}
			
			return sum;
		}
		numjs.RVector.prototype.normP = function(p){
			if (p<1)
			{
				throw new Error("p for p-norm must be greater than zero");
			}
			if (p==1){
				return this.norm1();
			}
			if (p==2){
				return this.norm();
			}
			
			var sum =0;
			for (var i=0; i<this.dimension();i++){
				sum += Math.pow(Math.abs(this.get(i)),p);
			}
			
			return Math.pow(sum, 1.0/p);
		}
		numjs.RVector.prototype.normInf = function()
		{
			var max=0;
			for(var i=0; i<this.dimension();i++){
				max = Math.max(max, Math.abs(this.get(i)));
			}
			
			return max;
		}
		//Need to implement RMatrix class first
		numjs.RVector.prototype.matrixMultiply = function (W) {
		    if (W.columnCount() > this.dimension())
		        throw new Error("The number of columns in the matrix must match the number of rows in the vector.");

		    var result = new numjs.RVector(W.rowCount());
		    var value = 0;

		    for (var i = 0; i < W.rowCount() ; i++) {
		        value = 0;
		        for (var j = 0; j < W.columnCount() ; j++) {
		            value += W.get(i, j) * this.get(j);
		        }
		        result.set(i, value);
		    }

		    return result;
		};
		//Need to implement RMatrix class first
		numjs.RVector.prototype.innerProductNorm = function (W) {
		    var result = 0;

		    var wx = this.matrixMultiply(W);

		    result = wx.dotProduct(wx);

		    result = Math.sqrt(result);

		    return result;
		};
		numjs.RVector.prototype.normalize = function () {
		    var norm = this.norm();

		    if (norm > 1E-13) {
		        this.scalarMultiply(1 / norm);
		    }
		};
        //Return a copy of the vector
		numjs.RVector.prototype.copy = function()
	    {
	        var vectorCopy = new numjs.RVector(this.dimension());
		
	        for(var i=0;i<this.dimension();i++)
	        {
	            vectorCopy.set(i, this.get(i));
	        }
	        return vectorCopy;
	    };
		numjs.RVector.prototype.toString = function () {
		    var vstring = "[";

		    for (var i = 0; i < this.dimension() ; i++) {
		        vstring += this.get(i);
		        if (i<this.dimension()-1)
                    vstring += ","
		    }

            vstring += "]"
		    return vstring;
		};
		numjs.RVector.prototype.isStrictEqual = function (v2) {
		    if (this.dimension() != v2.dimension())
		        return false;
		    for(var i=0;i<v2.dimension();i++)
		    {
		        if (this.get(i)!=v2.get(i))
		        {
		            return false;
		        }
		    }
		    return true;
                
		};
		numjs.RVector.prototype.isDeltaEqual = function (v2) {
		    if (this.dimension() != v2.dimension())
		        return false;
		    for (var i = 0; i < v2.dimension() ; i++) {
		        if (Math.abs(this.get(i) - v2.get(i)>this.DELTA)) {
		            return false;
		        }
		    }
		    return true;
		}

        //RMatrix
		numjs.RMatrix = function (n, m) {
			if (n<1)
				throw new Error("Matrix row length must be greater than 0.");
			if (m<1)
				throw new Error("Matrix column length must be greater than 0.");
			
			var rows=n;
			var cols=m;
			var matrix = [];

			this.DELTA = 1e-15;

			for (var i = 0; i < n; i++)
			{
			    matrix[i] = [];
				for(var j=0;j<m;j++)
				{
					matrix[i][j]=0;
				}
			}

			this.rowCount = function () {
			    return matrix.length;
			};

			this.columnCount = function () {
			    //even though this is a jagged array the columns of each row should match so
                //just get the number of columns in the first row
			    return matrix[0].length;
			};

		    //add these two methods to the constructor to make them public but don't use prototype
            //otherwise I would have to make the underlying matrix array public
			this.set = function (i, j, value) {
			    if (i < 0 || i > this.rowCount() - 1)
			        throw new Error("The row index must be non-negative and less than the row dimension of the matrix.");
			    if (j < 0 || j > this.columnCount() - 1)
			        throw new Error("The column index must be non-negative and less than the column dimension of the matrix. ");

			    matrix[i][j] = value;
			};

			this.get = function (i, j) {
			    if (i < 0 || i > this.rowCount() - 1)
			        throw new IndexOutOfBoundsException("The row index must be non-negative and less than the row dimension of the matrix.");
			    if (j < 0 || j > this.columnCount() - 1)
			        throw new IndexOutOfBoundsException("The column index must be non-negative and less than the column dimension of the matrix. ");

			    return matrix[i][j];
			};
		
		};
        //Return a copy of the vector
		numjs.RMatrix.prototype.copy = function () {
		    var matrixCopy = new numjs.RMatrix(this.rowCount(),this.columnCount());

		    for (var i = 0; i < this.rowCount() ; i++) {
		        for (var j = 0; j < this.columnCount() ; j++) {
		            matrixCopy.set(i, j, this.get(i,j));
		        }
		    }
		    return matrixCopy;
		};
		//Add matrix B to the current matrix
		numjs.RMatrix.prototype.add = function (B) {
		    if (B.rowCount() != this.rowCount())
				throw new Error("Matrix rows do not match.");
		    if (B.columnCount() != this.columnCount())
				throw new Error("Matrix columns do not match.");
			
			var value = 0;
			
			for (var i = 0; i < this.rowCount() ; i++)
			{
			    for(var j=0;j<this.columnCount();j++)
				{
					value = this.get(i,j) + B.get(i, j);
					this.set(i, j, value);
				}
			}
		};
		numjs.RMatrix.prototype.getColumnVector = function (col) {
		    if (col < 0 || col > this.rowCount() - 1)
		        throw new Error("The column index must be non-negative and less than the column dimension of the matrix.");

		    var vector = new numjs.RVector(cols);

		    for (var i = 0; i < this.rowCount() ; i++) {
		        vector.set(i, this.get(i,col));
		    }

		    return vector;
		};
        //Return a matrix row as a vector
		numjs.RMatrix.prototype.getRowVector = function (row) {
		    if (row < 0 || row > this.columnCount() - 1)
		        throw new Error("The row index must be non-negative and less than the row dimension of the matrix.");

		    var vector = new numjs.RVector(cols);

		    for (var j = 0; j < this.columnCount() ; j++) {
		        vector.set(j, this.get(row,j));
		    }

		    return vector;
		};
        //Return the resulting matrix of A*B
		numjs.RMatrix.prototype.multiply = function (B) {
		    if (this.columnCount() != B.rowCount())
		        throw new Error("The number of columns of matrix A must match rows of matrix B.");

		    var result = new numjs.RMatrix(this.rowCount(), B.columnCount());
		    var value = 0;

		    for (var i = 0; i < this.rowCount(); i++) {
		        for (var j = 0; j < B.columnCount() ; j++) {
		            value = 0.0;
		            for (var k = 0; k < this.columnCount(); k++) {
		                value += this.get(i,k) * B.get(k, j);
		            }
		            result.set(i, j, value);
		        }
		    }

		    return result;
		};
        //Subtract matrix B from the current matrix: A-B
		numjs.RMatrix.prototype.subtract = function(B)
		{
		    if (B.rowCount() != this.rowCount())
		        throw new Error("Matrix rows do not match.");
		    if (B.columnCount() != this.columnCount())
		        throw new Error("Matrix columns do not match.");

		    var value = 0;

		    for (var i = 0; i < this.rowCount() ; i++) {
		        for (var j = 0; j < this.columnCount(); j++) {
		            value = this.get(i,j) - B.get(i, j);
		            this.set(i, j, value);
		        }
		    }
		};
		numjs.RMatrix.prototype.transpose = function () {
		    var columnCount = this.columnCount();
		    var rowCount = this.rowCount();

		    var transpose = new numjs.RMatrix(columnCount, rowCount);
		    for (var i = 0; i < columnCount; i++) {
		        for (var j = 0; j < rowCount; j++) {
		            transpose.set(i, j, this.get(j,i));
		        }
		    }
		    return transpose;
		};
        //Returns the trace of the matrix
		numjs.RMatrix.prototype.trace = function () {
		    var trace = 0;
		    var minRowsCols = Math.min(this.rowCount(), this.columnCount());
		    for (var k = 0; k < minRowsCols; k++) {
		        trace += this.get(k,k);
		    }
		    return trace;
		};
        //Return the submatrix of A by eliminating the indicated row and column
		numjs.RMatrix.prototype.submatrix = function (row, col) {
		    if (row < 0 || row > this.rowCount() - 1)
		        throw new Error("The row index must be non-negative and less than the row dimension of the matrix.");
		    if (col < 0 || col > this.columnCount() - 1)
		        throw new Error("The column index must be non-negative and less than the column dimension of the matrix. ");

		    var rowCount = 0;
		    var colCount = 0;
		    var subm = new numjs.RMatrix(this.rowCount() - 1, this.columnCount() - 1);
		    for (var i = 0; i < this.rowCount() ; i++) {
		        if (i != row) {
		            colCount = 0;
		            for (var j = 0; j < this.columnCount(); j++) {
		                if (j != col) {
		                    subm.set(rowCount,colCount,this.get(i,j));
		                    colCount++;
		                }
		            }
		            rowCount++;
		        }
		    }
		    return subm;
		};
		numjs.RMatrix.prototype.toString = function () {
		    var mstring = "[";
		    for (var i = 0; i < this.rowCount() ; i++) {
		        mstring += "[";
		        for (var j = 0; j < this.columnCount() ; j++) {
		            mstring += this.get(i, j);
		            if (j < this.columnCount() - 1) {
		                mstring += ",";
		            }
		        }
		        mstring += "]";
		    }
		    mstring += "]"

		    return mstring;
		};

		numjs.RMatrix.prototype.isStrictEqual = function (B) {
		    if (this.rowCount()!=B.rowCount() || this.columnCount() != B.columnCount())
		    {
		        return false;
		    }
		    for(var i=0;i<this.rowCount();i++)
		    {
		        for(var j=0;j<this.columnCount();j++)
		        {
		            if (this.get(i, j) != B.get(i, j))
		                return false;
		        }
		    }
		    return true;
		};
		numjs.RMatrix.prototype.isDeltaEqual = function (B) {
		    if (this.rowCount() != B.rowCount() || this.columnCount() != B.columnCount()) {
		        return false;
		    }
		    for (var i = 0; i < this.rowCount() ; i++) {
		        for (var j = 0; j < this.columnCount() ; j++) {
		            if (Math.abs(this.get(i, j) - B.get(i, j))>this.DELTA)
		                return false;
		        }
		    }
		    return true;
		};

		numjs.RMatrix.prototype.LUPDecomposition = function () {
		    var pivotCount = 0;
		    var Acopy = this.copy();

		    var rowLength = Acopy.rowCount() - 1;
		    /*
             * pivot represents the permutation matrix.  It is implemented
             * as an array whose value indicates which column the 1 would appear.
             */
		    var pivot = [];
		    var p = 0;
		    var kp = 0;
		    var pik = 0;
		    var pikp = 0;
		    var aki = 0;
		    var akpi = 0;
		    var lup = [];

		    //Initialize the permutation matrix as the identity matrix
		    for (var j = 0; j <= rowLength; j++) {
		        pivot[j] = j;
		    }

		    for (var k = 0; k <= rowLength; k++) {
		        /*
                 * We want to find the permutation matrix that allows us 
                 * to avoid dividing by zero as well as avoid numerical instability
                 * from dividing by a relatively small number.
                 * We find the element with the largest absolute value of those in the
                 * current column (column k).  If all elements in the current column are zero
                 * then the matrix is singular and we throw an error
                 */
		        p = 0;
		        for (var i = k; i <= rowLength; i++) {
		            if (Math.abs(Acopy.get(i, k)) > p) {
		                p = Math.abs(Acopy.get(i, k));
		                kp = i;
		            }
		        }
		        if (p == 0) {
		            throw new Error("Matrix is singular.");
		        }
		        /*
                 * These lines update the pivot array
                 * by exchanging pivot[k] and pivot[kp].
                 */
		        pik = pivot[k];
		        pikp = pivot[kp];
		        pivot[k] = pikp;
		        pivot[kp] = pik;

		        /*
                 * Exchange rows k and kpi as determined by the pivot
                 */
		        for (var i = 0; i <= rowLength; i++) {
		            aki = Acopy.get(k, i);
		            akpi = Acopy.get(kp, i);
		            Acopy.set(k, i, akpi);
		            Acopy.set(kp, i, aki);
		            pivotCount++;
		        }

		        /*
                 * Compute the Schur complement
                 */
		        for (var i = k + 1; i <= rowLength; i++) {
		            Acopy.set(i, k, Acopy.get(i, k) / Acopy.get(k, k));
		            for (var j = k + 1; j <= rowLength; j++) {
		                Acopy.set(i, j, Acopy.get(i, j) - (Acopy.get(i, k) * Acopy.get(k, j)));
		            }
		        }
		    }
		    lup[0] = Acopy;
		    lup[1] = pivot;
		    lup[2] = pivotCount;
		    return lup;

		};
		numjs.RMatrix.prototype.CholeskyDecompositionGetL = function () {
		    var rowCount = this.rowCount();
		    var columnCount = this.columnCount();

		    if (rowCount != columnCount)
		        throw new Error("Matrix is not a square matrix.")

		    //L is the lower triangular matrix of the Cholesky decomposition
		    var L = new numjs.RMatrix(rowCount, columnCount);
		    var value = 0.0;
		    var ldiagonal = 0.0;
		    var lvalue = 0.0;

		    //lprimes are the vectors used in our computations in step 3 and 4. We just overwrite
		    //these values as needed as we iterate through the algorithm
		    var lprimes = new numjs.RVector(rowCount);

		    //Step 1
		    /*
             * All elements of L above the main diagonal will be equal to 0.
             * Set l00 equal to the square root of a00.  The remainder of the first column of L is 
             * the first column of A divided by l00.
             */
		    value = Math.sqrt(this.get(0, 0));
		    L.set(0, 0, value);
		    for (var i = 1; i < rowCount; i++) {
		        value = this.get(i, 0) / L.get(0, 0);
		        L.set(i, 0, value);
		    }

		    //Steps 2 through 5
		    for (var k = 1; k < rowCount; k++) {
		        //Step 2
		        //define all L primes
		        //L'm is a column vector of dimension k whose components are the first k-1
		        //elements of the mth row of L.
		        for (var m = k; m < rowCount; m++) {
		            var lm = new numjs.RVector(k);
		            for (var o = 0; o < k; o++) {
		                lm.set(o, L.get(m, o));
		            }
		            lprimes[m] = lm;
		        }
		        //Step 3
		        //compute the diagonal component
		        ldiagonal = Math.sqrt(this.get(k, k) - lprimes[k].dotProduct(lprimes[k]));
		        L.set(k, k, ldiagonal);
		        //Step 4
		        //compute all values below the diagonal
		        for (var p = k + 1; p < rowCount; p++) {
		            lvalue = (this.get(p, k) - lprimes[p].dotProduct(lprimes[k])) / ldiagonal;
		            L.set(p, k, lvalue);
		        }
		        ldiagonal = 0;

		    }

		    return L;

		};
		numjs.RMatrix.prototype.QRDecomposition = function () {
		    var n = this.rowCount();
		    var m = this.columnCount();
		    if (this.rowCount() < this.columnCount()) {
		        throw new Error("The number of rows must be greater than or equal to the number of columns.")
		    }

		    var Q = new numjs.RMatrix(n, m);
		    var R = new numjs.RMatrix(n, m);

		    //xvectors are the set of column vectors in A
		    //we perform the modified Gram-Schmidt process to these columns
		    var xvectors = [];
		    //qvectors are the orthonormal vectors resulting from the modified Gram-Schmidt process
		    //on the vectors of A
		    var qvectors = [];

		    //get the array of column vectors from A
		    for (var c = 0; c < m; c++) {
		        xvectors[c] = this.getColumnVector(c);
		    }


		};
		numjs.RMatrix.prototype.determinant = function () {
		    var result = 1.0;
		    var lup = this.LUPDecomposition();
		    var lu = lup[0];
		    var pivotCount = lup[2];

		    for (var i = 0; i < this.rowCount() ; i++) {
		        result *= lu.get(i, i);
		    }

		    result *= Math.pow(-1.0, pivotCount);

		    return result;

		};

        //Initialize an instance of the RMatrix class from a multi-dimensional array 
		numjs.convertArrayToRMatrix = function (data) {
		    if (data.length < 1)
		        throw new Error("The first dimension of the multidimensional array must be greater than 0.");
		    if (data[0].length < 1)
		        throw new Error("The second dimension of the multidimensional array must be greater than 0.");
		    for (var i = 0; i < data.length; i++) {
		        if (data[i].length != data[0].length)
		            throw new Error("The dimension of the each array in the multidimensional array must match.");
		    }

		    var rows = data.length;
		    var cols = data[0].length;
		    var matrix = new numjs.RMatrix(rows, cols);

		    for (var i = 0; i < data.length; i++) {
		        for (var j = 0; j < data[0].length; j++) {
		            matrix.set(i, j, data[i][j]);
		        }
		    }

		    return matrix;
		};
        //Return an identity matrix of size n by m
		numjs.identityMatrix = function (n, m) {
		    if (n < 1)
		        throw new Error("Matrix row length must be greater than 0.");
		    if (m < 1)
		        throw new Error("Matrix column length must be greater than 0.");

		    var A = new numjs.RMatrix(n, m);
		    for (var i = 0; i < n; i++) {
		        for (var j = 0; j < m; j++) {
		            if (i == j)
		                A.set(i, j, 1.0);
		            else
		                A.set(i, j, 0);
		        }
		    }

		    return A;
		};

		return numjs;
    }
    //define globally if it doesn't already exist
    if(typeof(numjs) === 'undefined'){
        window.numjs = define_numjs();
    }
    else{
        console.log("numjs already defined.");
    }
})(window);