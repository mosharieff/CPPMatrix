# C++ Matrix Algebra Library :microscope:
I wrote a simple matrix algebra library in C++ in this repository. Some of its features include matrix multiplication, inverse matrix, transpose, diagonalizing, exponent matrix, and matrix exponent. Additionally there is another class called 'statpack' in this repository which generates an ANOVA table for linear regression, calculates a min-variance/target-rate portfolio weights, along with calculating mean, variance, and standard deviation matrix formulas.


## Compiling Program :zap:
```sh
g++ test.cpp -std=c++11 or ./o.sh
```

<br/>

## StatPack Functions :computer:
Imports a .csv file and stores it in a string double vector
```C
vector<vector<string>> DataFrame(string filename);
```
Extracts a double vector from the dataframe given a column name
```C
vector<double> PullData(vector<vector<string>> df, string column);
```
Using the extracted double vector, a new dataframe could be generated using this
```C
vector<vector<double>> BuildDF(vector<vector<double>> x, vector<double> b);
```
Calculates percent change in a dataframe
```C
vector<vector<double>> RateOfReturn(vector<vector<double>> x, string side=["new","old"]);
```
Calculates percent change in a vector
```C
vector<double> PRateOfReturn(vector<double> x, string side=["new","old"]);
```

Adds a 1.0 to the beginning of each row of a data frame vector to be used to compute the intercept in a linear regerssion
```C
vector<vector<double>> RegIntercept(vector<vector<double>> x);
```
Computes the beta values for your linear/non-linear regression
```C
vector<vector<double>> Regression(vector<vector<double>> x, vector<double> y);      
```
Calculates the mean of each column in a matrix
```C
vector<vector<double>> Mean(vector<vector<double>> x);       
```
Calculates the covariance or correlation matrix
```C
vector<vector<double>> Variance(vector<vector<double>> x, string vtype="Covariance or Correlation");     
```
Generates an ANOVA table for a regression with R2, adjR2, FStat, and test statistic significant variables
```C
void ANOVA(vector<vector<double>> x, vector<double> y);
```
#### Portfolio Optimization :closed_lock_with_key:
Generates weights for a target rate portfolio and a min-variance portfolio
```C
vector<vector<double>> TargetRatePortfolio(vector<vector<double>> x, double r);
vector<vector<double>> MinVarPortfolio(vector<vector<double>> x);
```



<br/>

## Matrix Algebra Functions :books:
Print out a matrix or column vector
```C
void PRINTM(vector<vector<double>> x);
```
Generate the Identity matrix
```C
vector<vector<double>> IDENTITY(int n);
```
Matrix Multiplication Function
```C
vector<vector<double>> MMULT(vector<vector<double>> x, vector<vector<double>> y);
```
Calculate the Inverse of a matrix
```C
vector<vector<double>> INVERSE(vector<vector<double>> x);
```
Generate the Transpose of a matrix
```C
vector<vector<double>> TRANSPOSE(vector<vector<double>> x);
```
Diagonalize a matrix
```C
vector<vector<double>> DIAGONALIZE(vector<vector<double>> x, direction="Left or Right");
```
Calculate the determenant of a matrix
```C
double DETERMENANT(vector<vector<double>> x);
```
Calculate the trace of a matrix
```C
double TRACE(vector<vector<double>> x);
```
Raise a matrix to a power
```C
vector<vector<double>> POWER(vector<vector<double>> x, double power);
```
Raise number to matrix power
```C
vector<vector<double>> MPOWER(vector<vector<double>> x, double value);
```
Multiplies a coefficent to a matrix
```C
vector<vector<double>> COEF(double a, vector<vector<double>> x);
```
Adds or Subtracts a matrix by a vector
```C
vector<vector<double>> ADDMATRIX/SUBMATRIX(vector<vector<double>> x, vector<double> y);
```
Divides a matrix by a matrix
```C
vector<vector<double>> DIVISOR(vector<vector<double>> x, vector<vector<double>> y);
```
Converts a vector of length n into a matrix with dimensions (n x 1)
```C
vector<vector<double>> LINEAR(vector<double> x);
```
Converts a (n x 1) matrix to a single vector of length n
```C
vector<double> SINGLE(vector<vector<double>> x);
```
Generates the diagonal of a matrix
```C
vector<double> DIAGONAL(vector<vector<double>> x);
```
Generates a single vector with a series of set numbers
```C
vector<double> nums(double number, int n);
```
Calculates the mean of a single vector
```C
double mean(vector<double> y);
```

