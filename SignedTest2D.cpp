#include <iostream>
#include <string>
#include <sstream>
#include <iomanip> 
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <queue>
#include <stack>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <complex>
#include <numeric>
#include <fstream>

using namespace std;

long double variable3=0, variable4=0;
int variable1=0, variable2=0;
bool change = false;

string process(std::string const& s)
{
    string::size_type pos = s.find('/');
    if (pos != std::string::npos)
    {
        return s.substr(0, pos);
    }
    else
    {
        return s;
    }
}

void signtest(vector<int> v, int A) {

    ios::sync_with_stdio(false);
    cin.tie(0);
    // dimension of input matrix
    int n, i, j, k, s, zerocnt = 0;

    change = true;
    // pre-calculation needed for p values 

    n = A;

    vector<long double> f(n + 1);          // log of factorial: f[x] = log(x!)
    vector<long double> pow2(n + 1); // log power of 0.5:      pow2[x] = log(0.5^x)
    vector<vector<long double>> p(n + 1, vector<long double>(n + 1, 0));             // p[N][x] = 2 * P(i <= x) 

    f[0] = 0; pow2[0] = 0;
    for (i = 1; i <= n; i++) {
        f[i] = f[i - 1] + log(i);
        pow2[i] = pow2[i - 1] + log(.5);
    }

    for (i = 1; i <= n; i++) {
        for (j = 0; j <= i; j++) {
            p[i][j] = 2 * exp(f[i] - f[j] - f[i - j] + pow2[i]);
            if (j > 0) p[i][j] += p[i][j - 1];
        }
    }
    p[0][0] = 1;

    // prepare 2D prefix sum, sum[a][b] = # of 1s in v[i][j] where i <= a and j <= b (top left corner with a rows and b columns)

    // total # of 1 and -1 in the original matrix
    int total1 = 0, total_1 = 0;
    vector<int> sum(A + 1, 0);
    vector<int> sumneg(A + 1, 0);
    vector<int> sumzero(A + 1, 0);
    for (i = 1; i <= A; i++) {
        if (v[i] == 1) {
            total1++;
            sum[i] = sum[i - 1] + 1;
        }
        else if (v[i] == -1) {
            total_1++;
            sumneg[i] = sumneg[i - 1] + 1;
        }
        else {
            zerocnt++;
            sumzero[i] = sumzero[i - 1] + 1;
        }
    }

    long double pvalue = p[n - zerocnt][min(total_1, total1)];
    if (pvalue >= 0.05) {
        cout << "No change" << '\n';
        change = false;
        return;
    }

    // track the best answer
    long double bestp = 1E9;
    int minarea = n;

    long double pvalans = 0, coords1 = 0, coords2 = 0, newnumber = 0, newtotal = 0;

    // enumerate all window size a x b
    for (int a = 1; a <= A; a++) {
        if (a > minarea) continue;
        // enumerate all topleft cells
        for (int x = 1; x + a - 1 <= A; x++) {
            if (a == minarea) break;
            int xx = x + a - 1;

            // now calculate the # of 1s in the submatrix [x..xx]
            int c1, c_1, new1, new_1, nn, minn, c0;
            c1 = sum[xx] - sum[x - 1];
            c_1 = sumneg[xx] - sumneg[x - 1];
            c0 = sumzero[xx] - sumzero[x - 1];

            new1 = total1 - c1;
            new_1 = total_1 - c_1;
            minn = min(new1, new_1);
            nn = (n - zerocnt) - (a - c0);

            // calculate p
            long double pvalue = p[nn][minn];

            if ((pvalue >= 0.05) && (a < minarea)) {
                pvalans = pvalue;
                coords1 = x;
                coords2 = a;
                newnumber = minn;
                newtotal = nn;
                // log new best
                minarea = a;

                bestp = pvalue;
            }
        }
    }
    cout << "P-value " << setprecision(5) << fixed << pvalans << " found at by taking out " << coords2 << " cells in submatrix[" << coords1 << " .. " << coords1 + coords2 - 1 << "] with min(# of 1, # of -1) = " << newnumber << " out of " << newtotal << " cells\n";
  variable1 = coords1;
  variable2 = coords2;
  variable3 = pvalue;
  variable4 = pvalans;
}


int main() {
    clock_t start, end;
  
    ios_base::sync_with_stdio(false);
    cin.tie(0);

    int A, C;
    int location, age, gender;
    long long month = 1;
    int lastmonthcase[100] = {};
    string date, temp, lastmonth = "1";

    A = 9;
    C = 100;
  
    ofstream outfile;
    outfile.open("afile.txt");
  
    vector<int> v(A + 1, 0); //difference
    vector<vector<int>> v2(A + 1, vector<int>(C + 1, 0)); //original vectors

    ifstream fin("VancouverCovidData2.txt");

    fin >> temp >> temp >> temp >> temp >> temp;

    while (fin >> date >> location >> gender >> age >> temp) {

        if (date.rfind(lastmonth, 0) == 0) {
            v2[age][month]++;
            lastmonthcase[month]++;
        }
        else {
            lastmonth = process(date);
            month++;
        }
    }
    for (int i = 1; lastmonthcase[i] != 0; ++i) {
        for (int j = i + 1; lastmonthcase[j] != 0; ++j) {
            for (int a = 1; a <= A; ++a) {
                if ((long double)v2[a][i] / (long double)lastmonthcase[i] > (long double)v2[a][j] / (long double)lastmonthcase[j]) {
                    v[a] = 1;
                }
                else if ((long double)v2[a][i] / (long double)lastmonthcase[i] < (long double)v2[a][j] / (long double)lastmonthcase[j]) {
                    v[a] = -1;
                }
                else {
                    v[a] = 0;
                }
            }
            start = clock();
            signtest(v, A);
            end = clock();
            long double time_taken = (long double)(end - start) / double(CLOCKS_PER_SEC);
            if (change == true) outfile << i << "\t" << j << "\t" << variable3 << "\t" << variable4 << "\t" << variable2 << "\t" << variable1 << "\t" << variable1+variable2-1 << "\t" << time_taken << "\n";
          if (change == false) outfile << i << "\t" << j << "\t" << "---\t\t\t\t\t" << time_taken << '\n';
        }
    }
    outfile.close();
}

