#ifndef STATPACK_H
#define STATPACK_H
#define USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include "matrix.h"


class statpack {

    private:

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
            std::vector<std::vector<double>> ones = m.LINEAR(m.nums(1.0, p));
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
            Xu = m.SUBMATRIX(x, m.SINGLE(mean));
            res = m.MMULT(m.TRANSPOSE(Xu), Xu);
            if(vtype == "Correlation"){
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
            std::vector<double> mu = m.SINGLE(m.TRANSPOSE(Mean(x)));
            int n = cov.size();
            std::vector<std::vector<double>> res, C, weight;

            for(int i = 0; i < n; ++i){
                cov[i].push_back(mu[i]);
                cov[i].push_back(1.0);
            }

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

            C = m.LINEAR(B);

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
            cov = Variance(x, "Covariance");
            cov = m.COEF(2.0, cov);
            int n = cov.size();
            
            for(int i = 0; i < n; ++i){
                cov[i].push_back(1.0);
            }

            std::vector<double> temp, bump;
            for(auto & t : cov){
                temp.push_back(1.0);
                bump.push_back(0.0);
            }
            temp.push_back(0.0);
            bump.push_back(1.0);

            cov.push_back(temp);
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
            beta = Regression(x, y);

            yhat = m.MMULT(x, beta);
            double y_mu = m.mean(y);
            for(int i = 0; i < n; ++i){
                e.push_back({y[i] - yhat[i][0]});
            }
            double RSS = m.MMULT(m.TRANSPOSE(e), e)[0][0];
            double TSS = 0;
            for(int i = 0; i < n; ++i){
                TSS += pow(y[i] - y_mu, 2);
            }
            double df = x.size() - x[0].size();
            double ESS = TSS - RSS;
            double rsq = 1.0 - RSS/TSS;
            double mo = x[0].size() - 1;
            double adj_rsq = 1 - (1 - rsq)*(n - 1)/(n - mo - 1);
            double F = (ESS/mo)/(RSS/(n - mo - 1));
            double se = sqrt(RSS/(n - mo - 1));
            double factor = RSS / df;
            mtx = m.INVERSE(m.MMULT(m.TRANSPOSE(x), x));
            mtx = m.COEF(factor, mtx);
            sd = m.POWER(m.LINEAR(m.DIAGONAL(mtx)), 0.5);
            
            for(int i = 0; i < x[0].size(); ++i){
                tscore.push_back(beta[i][0]/sd[i][0]);
                pval.push_back(N(tscore[i]));
            }
            
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