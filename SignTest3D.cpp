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
#include <iostream>
#include <fstream>
#include <unistd.h>

using namespace std;

long double variable7=0, variable8=0;
int variable1=0, variable2=0, variable3=0, variable4=0, variable5=0, variable6=0;
bool change = false;
float maxn;

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

long double ReSign(vector<vector<vector<int>>> u, vector<vector<vector<int>>> w, int A, int B, int C, int totalmonthcaseu, int totalmonthcasew) {

  ios::sync_with_stdio(false);
  cin.tie(0);
  
  change = true;

  int n = A * B * C;
  int i, j, k, a, b, c, z, zerocnt = 0;
  float s, t;

  vector<vector<vector<float>>> v(A + 1, vector<vector<float> >(B + 1, vector <float>(C + 1, 0)));
  
  for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            if (totalmonthcaseu == 0 && w[i][j][k] == 0 || totalmonthcasew == 0 && u[i][j][k] == 0){
              v[i][j][k] = 0;
              zerocnt++;
            }
            else if (totalmonthcaseu == 0 && w[i][j][k] != 0){
              v[i][j][k] = -(double)(w[i][j][k]) / (double)totalmonthcasew;
            }
            else if (totalmonthcasew == 0 && u[i][j][k] != 0){
              v[i][j][k] = (double)(u[i][j][k]) / (double)totalmonthcaseu;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu < (long double)w[i][j][k] / (long double)totalmonthcasew) {
              v[i][j][k] = -(double)(w[i][j][k]) / totalmonthcasew - (double)(u[i][j][k]) / (double)totalmonthcaseu;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu > (long double)w[i][j][k] / (long double)totalmonthcasew) {
              v[i][j][k] = (double)(u[i][j][k]) / (double)totalmonthcaseu - (double)(w[i][j][k]) / totalmonthcasew;
            }
            else {
              v[i][j][k] = 0;
              zerocnt++;
            }
          }
        }
    }
  
    // pre-calculation needed for p values 

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

    // total # of 1 and -1 in the original matrix
    int total1 = 0, total_1 = 0;
    for (i = 1; i <= A; i++) {
      for (j = 1; j <= B; j++) {
        for (k = 1; k <= C; k++){
          if (v[i][j][k] == 1) total1++;
          else if (v[i][j][k] == -1) total_1++;
        }
      }
    }

    long double pvalue = p[n - zerocnt][min(total_1, total1)];
    if (pvalue >= 0.05) {
        change = false;
        return pvalue;
    }else{
      return 0;
    }
  }

void signtest(vector<vector<vector<int>>> u, vector<vector<vector<int>>> w, int A, int B, int C, int totalmonthcaseu, int totalmonthcasew) {

  ios::sync_with_stdio(false);
  cin.tie(0);
  
  change = true;

  int n = A * B * C;
  int i, j, k, a, b, c, z, zerocnt = 0;
  float s, t;

  vector<vector<vector<float>>> v(A + 1, vector<vector<float> >(B + 1, vector <float>(C + 1, 0)));
  
  for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            if (totalmonthcaseu == 0 && w[i][j][k] == 0 || totalmonthcasew == 0 && u[i][j][k] == 0){
              v[i][j][k] = 0;
              zerocnt++;
            }
            else if (totalmonthcaseu == 0 && w[i][j][k] != 0){
              v[i][j][k] = -(double)(w[i][j][k]) / (double)totalmonthcasew;
            }
            else if (totalmonthcasew == 0 && u[i][j][k] != 0){
              v[i][j][k] = (double)(u[i][j][k]) / (double)totalmonthcaseu;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu < (long double)w[i][j][k] / (long double)totalmonthcasew) {
              v[i][j][k] = -(double)(w[i][j][k]) / totalmonthcasew - (double)(u[i][j][k]) / (double)totalmonthcaseu;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu > (long double)w[i][j][k] / (long double)totalmonthcasew) {
              v[i][j][k] = (double)(u[i][j][k]) / (double)totalmonthcaseu - (double)(w[i][j][k]) / totalmonthcasew;
            }
            else {
              v[i][j][k] = 0;
              zerocnt++;
            }
          }
        }
    }
  
    // pre-calculation needed for p values 

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

    //number of units without including 0
    vector<vector<vector<int>>> nn (A + 1, vector<vector<int> >(B + 1, vector <int>(C + 1, 0)));
  
    for (i = 1; i <= A; i++) {
        for (j = 1; j <= B; j++) {
            z=0;
            for (k = 1; k <= C; k++){
              if (v[i][j][k] != 0) s++;
              nn[i][j][k] = nn[i - 1][j][k] + nn[i][j-1][k] - nn[i-1][j-1][k] + s;
            }
        }
    }

    // total # of 1 and -1 in the original matrix
    int total1 = 0, total_1 = 0;
    for (i = 1; i <= A; i++) {
      for (j = 1; j <= B; j++) {
        for (k = 1; k <= C; k++){
          if (v[i][j][k] == 1) total1++;
          else if (v[i][j][k] == -1) total_1++;
        }
      }
    }

    long double pvalue = p[n - zerocnt][min(total_1, total1)];
    if (pvalue >= 0.05) {
        change = false;
        cout << "No change" << '\n';
        return;
    }

    // track the best answer
    long double bestp = 1E9;
    int minarea = n;
    long double pvalans = 0, coords1 = 0, coords2 = 0, coords3 = 0, coords4 = 0, coords5 = 0, coords6 = 0, newnumber = 0, newtotal = 0;

    vector<vector<vector<int>>> u1(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
    vector<vector<vector<int>>> w1(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
  
    // enumerate all window size a x b
    for (a = 1; a <= A; a++) {
      if (pvalue >= 0.05) continue;
      for (b = 1; b <= B; b++) {
        for (c = 1; c <= C; c++){
          if (a * b * c > minarea) continue;

          // enumerate all topleft cells
          for (int x = 1; x + a - 1 <= A; x++) {
            for (int y = 1; y + b - 1 <= B; y++) {
              for (int z = 1; z + c - 1 <= C; z++){
                
                int xx = x + a - 1, yy = y + b - 1, zz = z + c - 1;

                int var1 = 0, var2 = 0;

                //REDO EVERYTHING WITH 0
                  for (int i = x; i <= xx; ++i){
                    for (int j = y; j <= yy; ++j){
                      for (int k = z; k <= zz; ++k){
                        var1 += u[i][j][k];
                        var2 += w[i][j][k];
                        u1[i][j][k] = u[i][j][k];
                        w1[i][j][k] = w[i][j][k];
                        u[i][j][k] = 0;
                        w[i][j][k] = 0;
                      }
                    }
                  }
                  long double pvalue2 = 0;
                  long double pvalue3 = ReSign(u, w, A, B, C, totalmonthcaseu - var1, totalmonthcasew - var2);
                  if (pvalue3 != 0) pvalue2 = pvalue3;

                  for (int i = x; i <= xx; ++i){
                    for (int j = y; j <= yy; ++j){
                      for (int k = z; k <= zz; ++k){
                        u[i][j][k] = u1[i][j][k];
                        w[i][j][k] = w1[i][j][k];
                      }
                    }
                  }
                
                int newn = nn[xx][yy][zz] - nn[x - 1][yy][zz] - nn[xx][y - 1][zz] - nn[xx][yy][z - 1] + nn[x - 1][y - 1][zz] + nn[xx][y - 1][z - 1] + nn[x - 1][yy][z - 1] - nn[x - 1][y - 1][z - 1];

                if ((pvalue2 >= 0.05) && (a * b * c < minarea)) {
                    pvalans = pvalue;
                    coords1 = x;
                    coords2 = y;
                    coords3 = a;
                    coords4 = b;
                    coords5 = z;
                    coords6 = c;
                    newnumber = maxn;
                    newtotal = (n - zerocnt) - newn;
                    // log new best
                    minarea = a * b * c;

                    bestp = pvalue;
                }
              }
            }
          }
        }
      }
    }
    cout << "P-value " << setprecision(5) << fixed << pvalans << " found at by taking out " << coords3 * coords4 * coords6 << " cells in submatrix[" << coords1 << " .. " << coords1 + coords3 - 1 << "][" << coords2 << " .. " << coords2 + coords4 - 1 << "][" << coords5 << " .. " << coords5 + coords6 - 1 << "] with min(# of 1, # of -1) = " << newnumber << " out of " << newtotal << " cells\n";
  variable1 = coords1;
  variable2 = coords2;
  variable3 = coords3;
  variable4 = coords4;
  variable5 = coords5;
  variable6 = coords6;
  variable7 = pvalue;
  variable8 = pvalans;
}


int main() {
    clock_t start, end;
  
    ios_base::sync_with_stdio(false);
    cin.tie(0);
  
    int A = 6, B = 9, C = 2, D = 100, month = 1;
  
    int location, age, gender;
    int totalmonthcase[D];
    string date, useless, lastmonth = to_string(month);

   ofstream outfile;
   outfile.open("afile.txt");
  
    ifstream fin("VancouverCovidData2.txt");
    fin >> useless >> useless >> useless >> useless >> useless;
  
    vector<vector<vector<int>>> u(A + 1, vector<vector<int> >(B + 1, vector <int>(C + 1, 0)));
      vector<vector<vector<int>>> w(A + 1, vector<vector<int> >(B + 1, vector <int>(C + 1, 0)));
    vector<vector<vector<vector<int>>>> v (A + 1,vector<vector<vector<int>>>(B + 1,vector<vector<int> >(C + 1,vector<int>(D + 1, 0))));

    while (fin >> date >> location >> gender >> age >> useless) {

        if (date.rfind(lastmonth, 0) == 0) {
            v[location][age][gender][month]++;
            totalmonthcase[month]++;
        }
        else {
            lastmonth = process(date);
            month++;
            v[location][age][gender][month]++;
            totalmonthcase[month]++;
        }
    }
  
    // for (int i = 1; totalmonthcase[i] != 0; ++i) {
    //     for (int j = i + 1; totalmonthcase[j] != 0; ++j) {
    ifstream fin2("edited3D.txt");
    int aa, bb, cc, dd, ee, ff;
    for (int i = 1; i == 1; ++i) {
        for (int j = 2; j <= D; ++j) {
          fin2 >> aa >> bb >> cc >> dd >> ee >> ff;
          
            for (int a = 1; a <= A; ++a) {
              for (int b = 1; b <= B; ++b) {
                for (int c = 1; c <= C; ++c){
                  u[a][b][c] = v[a][b][c][i];
                  w[a][b][c] = v[a][b][c][j];
                }
              }
            }
          
            start = clock();
            cout << totalmonthcase[i] << ' ' << totalmonthcase[j] << '\n';
            signtest(u, w, A, B, C, totalmonthcase[i], totalmonthcase[j]);
            end = clock();
          
            long double time_taken = (long double)(end - start) / double(CLOCKS_PER_SEC);
          
            if (change == true){
              outfile << i << "\t" << j << "\t" << variable7 << "\t" << variable8 << "\t" << variable3*variable4*variable6 << "\t" << variable1 << "\t" << variable1+variable3-1 << "\t" << variable2 << "\t" << variable2+variable4-1 << "\t" << variable5 << "\t" << variable5+variable6-1 << "\t" << time_taken << '\t';
              if (variable1 >= aa && variable1+variable3-1 <= bb && variable2 >= cc && variable2+variable4-1 <= dd && variable5 >= ee && variable5+variable6-1 <= ff){
                outfile << "yes\t" << (double)(variable1*(variable1+variable3-1)*variable2*(variable2+variable4-1)*variable5*(variable5+variable6-1))/(double)(aa*bb*cc*dd*ee*ff)*100 << "%\n";
              }
              else outfile << "no\n";
              
            }
            if (change == false) outfile << i << "\t" << j << "\t" << "---\t\t\t\t\t\t\t\t\t" << time_taken << '\n';
          
        }
    }
    
    outfile.close();
}
