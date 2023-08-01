# C++ Matrix Algebra Library
I wrote a simple matrix algebra library in C++ in this repository. Some of its features include matrix multiplication, inverse matrix, transpose, diagonalizing, exponent matrix, and matrix exponent. Additionally there is another class called 'statpack' in this repository which generates an ANOVA table for linear regression, calculates a min-variance/target-rate portfolio weights, along with calculating mean, variance, and standard deviation matrix formulas.
<br/>
## Matrix Algebra Functions
Print out a matrix or column vector
```C
PRINTM(vector<vector<double>> x);
```
Generate the Identity matrix
```C
IDENTITY(int n);
```
Matrix Multiplication Function
```C
MMULT(vector<vector<double>> x, vector<vector<double>> y);
```
Calculate the Inverse of a matrix
```C
INVERSE(vector<vector<double>> x);
```
Generate the Transpose of a matrix
```C
TRANSPOSE(vector<vector<double>> x);
```
Diagonalize a matrix
```C
DIAGONALIZE(vector<vector<double>> x, direction="Left or Right");
```
Calculate the determenant of a matrix
```C
DETERMENANT(vector<vector<double>> x);
```
Calculate the trace of a matrix
```C
TRACE(vector<vector<double>> x);
```
Raise a matrix to a power
```C
POWER(vector<vector<double>> x, double power);
```
Raise number to matrix power
```C
MPOWER(vector<vector<double>> x, double value);
```
