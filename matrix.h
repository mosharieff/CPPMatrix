#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>
#include <string>
#include <math.h>

class matrix {

    public:

        // Prints out a matrix
        void PRINTM(std::vector<std::vector<double>> x)
        {
            int n = x.size(), m = x[0].size();
            std::cout << "\nMatrix Dimensions: " << n << " x " << m << std::endl;

            // Loops through the double vector in order to print out output
            for(auto & i : x){
                for(auto & j : i){
                    std::cout << j << "\t";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        // Generates Identity Matrix based on inputed length
        std::vector<std::vector<double>> IDENTITY(int n)
        {
            std::vector<std::vector<double>> res;
            std::vector<double> temp;

            // Makes sure all diagonal values = 1 and other = 0
            for(int i = 0; i < n; ++i){
                temp.clear();
                for(int j = 0; j < n; ++j){
                    if(i == j){
                        temp.push_back(1.0);
                    } else {
                        temp.push_back(0.0);
                    }
                }
                res.push_back(temp);
            }

            return res;
        }

        // Multiplies Two Matrices
        std::vector<std::vector<double>> MMULT(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y)
        {
            std::vector<std::vector<double>> res;
            std::vector<double> temp;
            double total = 0;
            int m = x.size();
            int n = x[0].size();
            int o = y[0].size();

            // Multiplies matrix row i by second matrix column j
            for(int i = 0; i < m; ++i){
                temp.clear();
                for(int j = 0; j < o; ++j){
                    total = 0;
                    for(int k = 0; k < n; ++k){
                        total += x[i][k]*y[k][j];
                    }
                    temp.push_back(total);
                }
                res.push_back(temp);
            }

            return res;
        }

        // Inverse Matrix Function
        std::vector<std::vector<double>> INVERSE(std::vector<std::vector<double>> x)
        {
            std::vector<std::vector<double>> I = IDENTITY(x.size());

            double A, B;

            int n = x.size();

            // Gaussian Elimination on bottom half
            for(int i = 1; i < n; ++i){
                for(int j = 0; j < i; ++j){
                    A = x[i][j];
                    B = x[j][j];
                    for(int k = 0; k < n; ++k){
                        x[i][k] -= (A/B)*x[j][k];
                        I[i][k] -= (A/B)*I[j][k];
                    }
                }
            }

            // Gaussian Elimination on top half
            for(int i = 1; i < n; ++i){
                for(int j = 0; j < i; ++j){
                    A = x[j][i];
                    B = x[i][i];
                    for(int k = 0; k < n; ++k){
                        x[j][k] -= (A/B)*x[i][k];
                        I[j][k] -= (A/B)*I[i][k];
                    }
                }
            }

            // Divides identity matrix values from x diagonal values
            for(int i = 0; i < n; ++i){
                for(int j = 0; j < n; ++j){
                    I[i][j] /= x[i][i];
                }
            }

            return I;
        }

        // Transpose Matrix
        std::vector<std::vector<double>> TRANSPOSE(std::vector<std::vector<double>> x)
        {
            std::vector<std::vector<double>> res;
            std::vector<double> temp;

            int m = x.size(), n = x[0].size();

            // Shifts the rows and columns to be row=column, column=row
            for(int i = 0; i < n; ++i){
                temp.clear();
                for(int j = 0; j < m; ++j){
                    temp.push_back(x[j][i]);
                }
                res.push_back(temp);
            }

            return res;
        }

        // Diagonalize Matrix [Left or Right]
        std::vector<std::vector<double>> DIAGONALIZE(std::vector<std::vector<double>> x, std::string direction)
        {
            int n = x.size();

            double A, B;

            if(direction == "Right"){
                x = TRANSPOSE(x);
            }
            
            // Gaussian Elimination to create a diagonal matrix
            for(int i = 1; i < n; ++i){
                for(int j = 0; j < i; ++j){
                    A = x[i][j];
                    B = x[j][j];
                    for(int k = 0; k < n; ++k){
                        x[i][k] -= (A/B)*x[j][k];
                    }
                }
            }

            if(direction == "Right"){
                x = TRANSPOSE(x);
            }

            return x;
        }

        // Determenant of Matrix
        double DETERMENANT(std::vector<std::vector<double>> x)
        {
            x = DIAGONALIZE(x, "Left");
            x = DIAGONALIZE(x, "Right");
            double total = 1.0;
            // Takes the product of diagonalized matrix to compute determenant
            for(int i = 0; i < x.size(); ++i){
                total *= x[i][i];
            }
            return total;
        }

        // Calculates the trace of a matrix
        double TRACE(std::vector<std::vector<double>> x)
        {
            double total = 0;
            // Takes the sum of the diagonals
            for(int i = 0; i < x.size(); ++i){
                total += x[i][i];
            }
            return total;
        }

        // Raises matrix to a power
        std::vector<std::vector<double>> POWER(std::vector<std::vector<double>> x, double power)
        {
            // Raises a matrix to an inputted power
            for(int i = 0; i < x.size(); ++i){
                for(int j = 0; j < x[0].size(); ++j){
                    x[i][j] = pow(x[i][j], power);
                }
            }
            return x;
        }

        // Raises number to the power of a matrix
        std::vector<std::vector<double>> MPOWER(std::vector<std::vector<double>> x, double value)
        {
            // Creates a matrix exponent for the inputted value to be given the power of
            for(int i = 0; i < x.size(); ++i){
                for(int j = 0; j < x[0].size(); ++j){
                    x[i][j] = pow(value, x[i][j]);
                }
            }
            return x;
        }

        // Multiplies matrix by coefficent
        std::vector<std::vector<double>> COEF(double a, std::vector<std::vector<double>> x)
        {
            int n = x.size();
            int m = x[0].size();

            // Makes sure each element of the matrix is multiplied by the factor
            for(int i = 0; i < n; ++i){
                for(int j = 0; j < m; ++j){
                    x[i][j] *= a;
                }
            }

            return x;
        }

        // Adds a matrix and a vector
        std::vector<std::vector<double>> ADDMATRIX(std::vector<std::vector<double>> x, std::vector<double> y)
        {
            // Adds a single vector to a matrix
            for(int i = 0; i < x.size(); ++i){
                for(int j = 0; j < y.size(); ++j){
                    x[i][j] += y[j];
                }
            }
            return x;
        }

        // Subtracts a matrix and a vector
        std::vector<std::vector<double>> SUBMATRIX(std::vector<std::vector<double>> x, std::vector<double> y)
        {
            // Subtracts a single vector from a matrix
            for(int i = 0; i < x.size(); ++i){
                for(int j = 0; j < y.size(); ++j){
                    x[i][j] -= y[j];
                }
            }
            return x;
        }

        // Divides a matrix by a matrix
        std::vector<std::vector<double>> DIVISOR(std::vector<std::vector<double>> x, std::vector<std::vector<double>> y)
        {
            std::vector<std::vector<double>> res;
            std::vector<double> temp;

            // Divides matrix x by matrix y (used in correlation matrix computation)
            for(int i = 0; i < x.size(); ++i){
                temp.clear();
                for(int j = 0; j < x[0].size(); ++j){
                    temp.push_back(x[i][j]/y[i][j]);
                }
                res.push_back(temp);
            }
            return res;
        }

        // Converts single vector to linear vector
        std::vector<std::vector<double>> LINEAR(std::vector<double> x)
        {
            std::vector<std::vector<double>> res;
            for(auto & t : x){
                res.push_back({t});
            }
            return res;
        }

        // Converts a linear vector to a single vector
        std::vector<double> SINGLE(std::vector<std::vector<double>> x)
        {
            std::vector<double> res;
            for(auto & t : x){
                res.push_back(t[0]);
            }
            return res;
        }

        // Takes the diagonal of a matrix
        std::vector<double> DIAGONAL(std::vector<std::vector<double>> x)
        {
            std::vector<double> res;
            for(int i = 0; i < x.size(); ++i){
                res.push_back(x[i][i]);
            }
            return res;

        }

        // Generates a set of numbers
        std::vector<double> nums(double number, int n)
        {
            std::vector<double> z;
            for(int i = 0; i < n; ++i){
                z.push_back(number);
            }
            return z;
        }

        // Calculates the mean of a single vector (non linearized)
        double mean(std::vector<double> y)
        {
            double total = 0;
            for(auto & r : y){
                total += r;
            }
            return total / y.size();
        }
};



#endif