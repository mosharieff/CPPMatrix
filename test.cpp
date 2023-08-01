#include <iostream>
#include "matrix.h"
#include "statpack.h"
#include <vector>

int main()
{
    matrix m;
    statpack s;

    std::vector<std::vector<double>> x = {{2, 3, -15, 5},
                                          {3, 5, -2, 6},
                                          {5, 100, -8, 3},
                                          {1, 2, 1, 1},
                                          {2, 666, -9, 2},
                                          {7, 6, 4, 5},
                                          {4, 2, -5, -7},
                                          {9, 3, -3, 6},
                                          {8, 9, 10, 1},
                                          {2, 15, 2, 0},
                                          {15, 9, 0, 8}}; 
                                          

    
    std::vector<double> y = {20, 30, 50, 18, 23, 45, 66, 92, 84, 21, 55};

    x = s.RegIntercept(x);
    s.ANOVA(x, y);


    return 0;
}