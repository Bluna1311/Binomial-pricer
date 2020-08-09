#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;



 // return the payoff of call
double payoff(double S, double k)
{
    return max(S-k, 0.);
}

// this will do the stock price tree and the go through a vanilla european option price at each step. will return the option price
double binomialTree(double S0, double k, double T, double r, double sigma, int n)
{

    double dt, u, d, q;
    dt = T / n;
    u = exp(sigma * sqrt(dt));
    d = exp(-sigma * sqrt(dt));
    q = (exp(r*dt)-d) / (u-d);
    // 2d vector of size n+1 for the time steps and other vector to have the size of n+1 for the same.
    vector<vector<double>> stockTree(n+1, vector<double>(n+1));
    
    for (int i = 0; i <= n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            stockTree[i][j] = S0 * pow(u, j) * pow(d, i-j);
        }
    }
    vector<vector<double>> valueTree(n+1, vector<double>(n+1));
    for (int j = 0; j <= n; j++)
    {
        valueTree[n][j] = payoff(stockTree[n][j],k);
    }
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = 0; j <= i; j++)
        {
            valueTree[i][j] = exp(-r*dt) * (q * valueTree[i+1][j+1] + (1-q) * valueTree[i+1][j]);
        }
    }
    return valueTree[0][0];
}

int main()
{



    
    double S0 = 190, k = 190., T = 1., r = 0.06, sigma = 0.2;
    // declare and initialise tree paramaters (steps in tree)
    
    // open output stream
    ofstream output("C://Users//Owner//Desktop//output C++//binomial-convergence.csv");
    // check it is open
    if (output.is_open())
    {
        cout << "File opened successfully" << endl;
        // output various n and V_n to file
        for (int n = 3; n <= 1000; n++)
        {
            output << n << " , " << binomialTree(S0, k, T, r, sigma, n) << endl;
        }
        cout << "File write complete" << endl;
    }
    else
    {
        cout << "File could not be opened for some reason.\n";
    }









}
   