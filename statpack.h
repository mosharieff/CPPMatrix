#ifndef STATPACK_H
#define STATPACK_H
#define USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include "matrix.h"
#include <sstream>
#include <fstream>


class statpack {

    private:

        // Cumulative Normal Distribution function
        double N(double x)
        {
            auto fx = [](double x){
                return exp(-pow(x, 2)/2.0)/(sqrt(2*M_PI));
            };
            double p = 0, x0 = 0, C;
            int n = 101;
            x = abs(x);
            
            double dX = (x - x0)/((double) n - 1);
            
            for(int i = 0; i < n; ++i){
                if(i == 0 || i == n - 1){
                    C = 1.0;
                } else if(i % 2 == 0){
                    C = 2.0;
                } else {
                    C = 4.0;
                }
                p += C*fx(x0 + i*dX);
            }
            p *= dX/3.0;
            return 1 - (0.5 + p);
        }

        matrix m;

    public:

        // Calculates percent change in a dataframe
        std::vector<std::vector<double>> RateOfReturn(std::vector<std::vector<double>> x, std::string side){
            std::vector<std::vector<double>> result;
            std::vector<double> temp;

            for(int i = 1; i < x.size(); ++i){
                temp.clear();
                for(int j = 0; j < x[0].size(); ++j){
                    if(side == "new"){
                        temp.push_back(x[i-1][j]/x[i][j] - 1);
                    } else {
                        temp.push_back(x[i][j]/x[i-1][j] - 1);
                    }
                }
                result.push_back(temp);
            }

            return result;
        }

        // Calculates percent change for a vector
        std::vector<double> PRateOfReturn(std::vector<double> x, std::string side)
        {
            std::vector<double> res;
            for(int i = 1; i < x.size(); ++i){
                if(side == "new"){
                    res.push_back(x[i-1]/x[i] - 1);
                } else {
                    res.push_back(x[i]/x[i-1] - 1);
                }
            }
            return res;
        }

        // Reads a CSV file
        std::vector<std::vector<std::string>> DataFrame(std::string filename)
        {
            std::string temp;
            std::ifstream dataset(filename);

            std::vector<std::vector<std::string>> framedata;
            std::vector<std::string> tempdata;

            while(std::getline(dataset, temp)){
                std::string row;
                std::istringstream linear(temp);
                tempdata.clear();
                while(std::getline(linear, row, ',')){
                    tempdata.push_back(row);
                }
                framedata.push_back(tempdata);
            }

            return framedata;
        }

        // Extracts a double vector from the dataset
        std::vector<double> PullData(std::vector<std::vector<std::string>> df, std::string column)
        {
            int colindex = 0;
            for(auto & c : df[0]){
                if(c == column){
                    break;
                }
                colindex += 1;
            }
            std::vector<double> result;
            df.erase(df.begin());
            for(auto & k : df){
                result.push_back(atof(k[colindex].c_str()));
            }
            return result;
        }

        // Builds numeric dataset for Regression
        std::vector<std::vector<double>> BuildDF(std::vector<std::vector<double>> x, std::vector<double> b)
        {
            x.push_back(b);
            return x;
        }

        // Adds a '1' at the beginning of the matrix in order to 
        // include an intercept in your regression
        std::vector<std::vector<double>> RegIntercept(std::vector<std::vector<double>> x)
        {
            std::vector<std::vector<double>> reg;
            std::vector<double> temp;
            for(int i = 0; i < x.size(); ++i){
                temp.clear();
                temp.push_back(1.0);
                for(auto & t : x[i]){
                    temp.push_back(t);
                }
                reg.push_back(temp);
            }

            return reg;
        }
        
        // Calculates your beta coefficients
        std::vector<std::vector<double>> Regression(std::vector<std::vector<double>> x, std::vector<double> y)
        {
            std::vector<std::vector<double>> Y = m.LINEAR(y);
            std::vector<std::vector<double>> XTX, XTY;
            // Computes pieces of the regression formula B = (XTX)^-1 XTy
            XTX = m.MMULT(m.TRANSPOSE(x), x);
            XTX = m.INVERSE(XTX);
            XTY = m.MMULT(m.TRANSPOSE(x), Y);
            return m.MMULT(XTX, XTY);
        }

        // Calculates your mean matrix
        std::vector<std::vector<double>> Mean(std::vector<std::vector<double>> x)
        {
            int p = x.size(), n = x[0].size();
            double factor = 1.0 / (double) p;
            // Generates a vector where each value is 1
            std::vector<std::vector<double>> ones = m.LINEAR(m.nums(1.0, p));
            // Computes the mean of each column of the matrix
            ones = m.TRANSPOSE(ones);
            ones = m.MMULT(ones, x);
            ones = m.COEF(factor, ones);
            return ones;
        }

        // Calculates your covariance matrix or correlation matrix
        std::vector<std::vector<double>> Variance(std::vector<std::vector<double>> x, std::string vtype)
        {
            std::vector<std::vector<double>> res, Xu, px, sd, sdx;
            std::vector<std::vector<double>> mean;
            mean = m.TRANSPOSE(Mean(x));
            Xu = m.SUBMATRIX(x, m.SINGLE(mean)); // Subtracts the mean from the inputted matrix
            res = m.COEF(1.0/((double) x.size() - 1), m.MMULT(m.TRANSPOSE(Xu), Xu)); // Computes covariance
            if(vtype == "Correlation"){
                // Takes diagonal of covar matrix, square roots it, then multiplies the vector by its transpose and 
                // calculates the correlation matrix through the matrix divider function
                px = m.LINEAR(m.DIAGONAL(res));
                sd = m.POWER(px, 0.5);
                sdx = m.MMULT(sd, m.TRANSPOSE(sd));
                res = m.DIVISOR(res, sdx);
                return res;
            }

            return res;
        }

        // Calculates the weights of a portfolio with a set target rate
        std::vector<std::vector<double>> TargetRatePortfolio(std::vector<std::vector<double>> x, double r)
        {
            std::vector<std::vector<double>> cov = Variance(x, "Covariance");
            cov = m.COEF(2.0, cov); // Multiplies covariance matrix by 2
            std::vector<double> mu = m.SINGLE(m.TRANSPOSE(Mean(x)));
            int n = cov.size(); // Takes the length of the covariance matrix
            std::vector<std::vector<double>> res, C, weight;

            // Adds the mean and 1 vector to the covariance matrix
            for(int i = 0; i < n; ++i){
                cov[i].push_back(mu[i]);
                cov[i].push_back(1.0);
            }

            // Adds the mean and 1 vector transpose to the bottom of the matrix
            std::vector<double> tA, tB, B;
            for(int i = 0; i < n; ++i){
                tA.push_back(mu[i]);
                tB.push_back(1.0);
                B.push_back(0.0);
            }
            tA.push_back(0.0);
            tA.push_back(0.0);
            tB.push_back(0.0);
            tB.push_back(0.0);

            cov.push_back(tA);
            cov.push_back(tB);

            B.push_back(r);
            B.push_back(1.0);

            C = m.LINEAR(B); // Converts single vector into column vector

            // Calculates the weights using lagrange optimization
            res = m.MMULT(m.INVERSE(cov), C);
            for(int i = 0; i < n; ++i){
                weight.push_back({res[i]});
            }

            return weight;
        }

        // Calculates the weights for a min-variance portfolio
        std::vector<std::vector<double>> MinVarPortfolio(std::vector<std::vector<double>> x)
        {
            std::vector<std::vector<double>> res, cov;

            // Computes covariance matrix, multiplies it by 2.0 and takes the size of it
            cov = Variance(x, "Covariance");
            cov = m.COEF(2.0, cov);
            int n = cov.size();
            
            // Adds a 1 vector to the end of the covariance matrix
            for(int i = 0; i < n; ++i){
                cov[i].push_back(1.0);
            }

            // Computes bottom of covariance matrix with 1 vector and converts the bump vector to a column vector
            std::vector<double> temp, bump;
            for(auto & t : cov){
                temp.push_back(1.0);
                bump.push_back(0.0);
            }
            temp.push_back(0.0);
            bump.push_back(1.0);

            cov.push_back(temp);

            // Calculates the weights using lagrange optimization
            cov = m.MMULT(m.INVERSE(cov), m.LINEAR(bump));

            for(int i = 0; i < n; ++i){
                res.push_back({cov[i]});
            }

            return res;
        }

        // Generates an ANOVA table for regression
        void ANOVA(std::vector<std::vector<double>> x, std::vector<double> y)
        {
            std::vector<std::vector<double>> beta, yhat, e, mtx, sd;
            std::vector<double> tscore, pval;
            int n = y.size();

            // Calculates beta values
            beta = Regression(x, y);

            // Computes inputs
            yhat = m.MMULT(x, beta);
            double y_mu = m.mean(y);
            for(int i = 0; i < n; ++i){
                e.push_back({y[i] - yhat[i][0]});
            }

            // Calculates RSS, TSS, degrees of freedom, and ESS
            double RSS = m.MMULT(m.TRANSPOSE(e), e)[0][0];
            double TSS = 0;
            for(int i = 0; i < n; ++i){
                TSS += pow(y[i] - y_mu, 2);
            }
            double df = x.size() - x[0].size();
            double ESS = TSS - RSS;

            // Computes R-Squared and Adjusted R-Squared
            double rsq = 1.0 - RSS/TSS;
            double mo = x[0].size() - 1;
            double adj_rsq = 1 - (1 - rsq)*(n - 1)/(n - mo - 1);

            // Computes F statistic
            double F = (ESS/mo)/(RSS/(n - mo - 1));

            // Computes standard error
            double se = sqrt(RSS/(n - mo - 1));

            // Calculates standard deviation measures of each beta and is used to compute test statistic
            double factor = RSS / (df - 2);
            mtx = m.INVERSE(m.MMULT(m.TRANSPOSE(x), x));
            mtx = m.COEF(factor, mtx);
            sd = m.POWER(m.LINEAR(m.DIAGONAL(mtx)), 0.5);
            
            // Computation of test statistic and pvalue (using normal distribution rather than t-distribution)
            for(int i = 0; i < x[0].size(); ++i){
                tscore.push_back(beta[i][0]/sd[i][0]);
                pval.push_back(N(tscore[i]));
            }
            
            // Prints out ANOVA table
            std::cout << std::endl;
            std::cout << "ANOVA Regression Table\n" << std::endl;
            std::cout << "R-Squared: " << rsq << std::endl;
            std::cout << "Adj-R-Squared: " << adj_rsq << std::endl;
            std::cout << "F-Statistic: " << F << std::endl;
            std::cout << "Standard Error: " << se << std::endl;
            std::cout << std::endl;

            for(int i = 0; i < tscore.size(); ++i){
                std::cout << "Variable: " << i << std::setprecision(3) <<  "\t Beta: " << beta[i][0] << "\t StdErr: " << std::setprecision(3) << sd[i][0] << "\t TestStat: " << std::setprecision(3) << tscore[i] << "\t PValue: " << std::setprecision(3) << pval[i] << std::endl;
            }
        }
};

#endif
