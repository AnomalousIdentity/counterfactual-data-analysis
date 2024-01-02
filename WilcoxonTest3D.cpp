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

long double variable7=0, variable8=0;
int variable1=0, variable2=0, variable3=0, variable4=0, variable5=0, variable6=0;
double variable9 = 0;
bool change = false;
float maxn;

double R[1000] = { 0 };

double rmse(vector<vector<vector<int>>> predicted, vector<vector<vector<int>>> actual, int A, int B, int C)
{
	double sum = 0;

	for (int i = 1; i <= A; i++){
    for (int j = 1; j <= B; j++){
      for (int k = 1; k <= C; k++){
		    sum += pow((double)predicted[i][j][k]-(double)actual[i][j][k], 2);
      }
    }
  }
	return sqrt(sum / A*B*C);
}

long double normalCDF(long double value)
{
   return 0.5*erfc(-value * M_SQRT1_2);
}

bool sortbysec(const pair<long double, long double>& a, const pair<long double, long double>& b)
{
    return (a.first < b.first);
}

string process(string const& s)
{
    string::size_type pos = s.find('/');
    if (pos != string::npos)
    {
        return s.substr(0, pos);
    }
    else
    {
        return s;
    }
}

void rankit(double A[], int n) {
    vector <pair<double, double>> T(n);

    for (int i = 0; i < n; i++)
    {
        T[i].first = A[i];
        T[i].second = i;
    }

    sort(T.begin(), T.end(), sortbysec);

    double rank = 1, m = 1, i = 0;

    while (i < n) {
        double j = i;
        while (j < n - 1 && T[j].first == T[j + 1].first)
            j += 1;

        m = j - i + 1;

        for (int k = 0; k < m; k++) {
            int idx = T[i + k].second;
            R[idx] = (double)(rank + (m - 1) * 0.5);
        }
        rank += m;
        i += m;
    }
}

long double RerankVector(int A, int B, int C, int totalmonthcaseu, int totalmonthcasew, vector<vector<vector<int>>> u, vector<vector<vector<int>>> w){
  
    ios::sync_with_stdio(false);
    cin.tie(0);
  
    int n = A * B * C;
    int i, j, k, a, b, z;
    float s, t;

    double list[108];
    int negs[108] = {0};
    int zerocnt = 0;

    vector<vector<vector<float>>> v(A + 1, vector<vector<float> >(B + 1, vector <float>(C + 1, 0)));
  
    for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            if ((long double)u[i][j][k] / (long double)totalmonthcaseu < (long double)w[i][j][k] / (long double)totalmonthcasew) {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = (double)(w[i][j][k]) / totalmonthcasew - (double)(u[i][j][k]) / (double)totalmonthcaseu;
                negs[(i - 1) * B * C + (j - 1) * C + (k - 1)] = 1;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu > (long double)w[i][j][k] / (long double)totalmonthcasew) {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = (double)(u[i][j][k]) / (double)totalmonthcaseu - (double)(w[i][j][k]) / totalmonthcasew;
            }
            else {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = 0;
                zerocnt++;
            }
          }
        }
    }
  
    // rank list
    rankit(list, sizeof(list) / sizeof(list[0]));

    // adjustments to rank, nullifying zeros
    double temp = n;
    for (int i = 0; i < n; i++) if (R[i] < temp) temp = R[i];
    for (int i = 0; i < n; i++) {

        if (zerocnt > 0 && R[i] != temp) {
            R[i] -= zerocnt;
        }
        else if (zerocnt > 0) {
            R[i] = 0;
        }
        else break;
    }

    //reinsert back into vector, add negatives

    for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            (negs[(i - 1) * B * C + (j - 1) * C + (k - 1)] == 1) ? (v[i][j][k] = -R[(i - 1) * B * C + (j - 1) * C + (k - 1)]) : (v[i][j][k] = R[(i - 1) * B * C + (j - 1) * C + (k - 1)]);
          }
        }
    }

    // total # of positives and negatives in the original matrix
    float totalpos = 0, totalnegs = 0;
    for (i = 1; i <= A; i++) {
        for (j = 1; j <= B; j++) {
          for (k = 1; k <= C; k++){
            if (v[i][j][k] > 0) {
                totalpos += v[i][j][k];
            }
            else if (v[i][j][k] < 0) {
                totalnegs -= v[i][j][k];
            }
          }
        }
    }
  
    // get z scores
    int nnozero = n - zerocnt;
    long double zscore = (max(totalpos, totalnegs) - nnozero * (nnozero + 1) / 4) / sqrt(nnozero * (nnozero + 1) * (2 * nnozero + 1) / 24);
    long double pvalue = 2*(1-normalCDF(zscore));
    if (pvalue >= 0.05) {
      maxn = max(totalpos, totalnegs);
      return pvalue;
    }else{
      return 0;
    }
}

void WilcoxonTest(vector<vector<vector<int>>> u, vector<vector<vector<int>>> w, int A, int B, int C, int totalmonthcaseu, int totalmonthcasew)
{
    ios::sync_with_stdio(false);
    cin.tie(0);
  
    change = true;
  
    int n = A * B * C;
    int i, j, k, a, b, z;
    float s, t;

    double list[108];
    int negs[108] = {0};
    int zerocnt = 0;
  
    for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            if ((long double)u[i][j][k] / (long double)totalmonthcaseu < (long double)w[i][j][k] / (long double)totalmonthcasew) {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = (double)(w[i][j][k]) / totalmonthcasew - (double)(u[i][j][k]) / (double)totalmonthcaseu;
                negs[(i - 1) * B * C + (j - 1) * C + (k - 1)] = 1;
            }
            else if ((long double)u[i][j][k] / (long double)totalmonthcaseu > (long double)w[i][j][k] / (long double)totalmonthcasew) {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = (double)(u[i][j][k]) / (double)totalmonthcaseu - (double)(w[i][j][k]) / totalmonthcasew;
            }
            else {
                list[(i - 1) * B * C + (j - 1) * C + (k - 1)] = 0;
                zerocnt++;
            }
          }
        }
    }
  
    // rank list
    rankit(list, sizeof(list) / sizeof(list[0]));

    // adjustments to rank, nullifying zeros

    double temp = n;
    for (int i = 0; i < n; i++) if (R[i] < temp) temp = R[i];
    for (int i = 0; i < n; i++) {

        if (zerocnt > 0 && R[i] != temp) {
            R[i] -= zerocnt;
        }
        else if (zerocnt > 0) {
            R[i] = 0;
        }
        else break;
    }

    //reinsert back into vector, add negatives


    vector<vector<vector<float>>> v(A + 1, vector<vector<float> >(B + 1, vector <float>(C + 1, 0)));
    for (int i = 1; i <= A; i++) {
        for (int j = 1; j <= B; j++) {
          for (int k = 1; k <= C; k++){
            (negs[(i - 1) * B * C + (j - 1) * C + (k - 1)] == 1) ? (v[i][j][k] = -R[(i - 1) * B * C + (j - 1) * C + (k - 1)]) : (v[i][j][k] = R[(i - 1) * B * C + (j - 1) * C + (k - 1)]);
          }
        }
    }

    //number of units without including 0
    
    vector<vector<vector<float>>> nn(A + 1, vector<vector<float> >(B + 1, vector <float>(C + 1, 0)));

    for (i = 1; i <= A; i++) {
        for (j = 1; j <= B; j++) {
            z=0;
            if (v[i][j][k] != 0) z++;
            nn[i][j][k] = nn[i - 1][j][k] + nn[i][j-1][k] - nn[i-1][j-1][k] + z;
        }
    }

    // total # of positives and negatives in the original matrix
    float totalpos = 0, totalnegs = 0;
    for (i = 1; i <= A; i++) {
        for (j = 1; j <= B; j++) {
          for (k = 1; k <= C; k++){
            if (v[i][j][k] > 0) {
                totalpos += v[i][j][k];
            }
            else if (v[i][j][k] < 0) {
                totalnegs -= v[i][j][k];
            }
          }
        }
    }

    // get z scores
  
    int nnozero = n - zerocnt;
    long double zscore = (max(totalpos, totalnegs) - nnozero * (nnozero + 1) / 4) / sqrt(nnozero * (nnozero + 1) * (2 * nnozero + 1) / 24);
    long double pvalue = 2*(1-normalCDF(zscore));
    if (pvalue >= 0.05) {
        // for (k = 1; k <= C; k++){
        //   for (i = 1; i <= A; i++){
        //     for (j = 1; j <= B; j++){
        //       cout << v[i][j][k] << ' ';
        //     }
        //     cout << '\n';
        //   }
        //   cout << '\n';
        // }
        variable9 = rmse(u, w, A, B, C);
        cout << "no change needed, p value = " << pvalue << ", "<< max(totalpos, totalnegs) << " total out of " << nnozero << " cells." << '\n';
        change = false;
        return;
    }

    cout << "--------------------------------------" << '\n';
    for (k = 1; k <= C; k++){
      for (i = 1; i <= A; i++){
        for (j = 1; j <= B; j++){
          cout << v[i][j][k] << ' ';
        }
        cout << '\n';
      }
      cout << '\n';
    }
  
    // track the best answer

    float minarea = n+1;

    long double zvalans = 0, rmsevalans = 0, pvalans = 0, coords1 = 0, coords2 = 0, coords3 = 0, coords4 = 0, coords5 = 0, coords6 = 0, newnumber = 0, newtotal = 0;

    vector<vector<vector<int>>> u1(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
    vector<vector<vector<int>>> w1(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
    // enumerate all window size a x b

    for (a = 1; a <= A; a++) {
        for (b = 1; b <= B; b++) {
          for (int c = 1; c <= C; c++){
            if (a * b * c > minarea || a > minarea || b > minarea || c > minarea) continue;

            // enumerate all topleft cells
            for (int x = 1; x + a - 1 <= A; x++) {
              for (int y = 1; y + b - 1 <= B; y++) {
                for (int z = 1; z + c - 1 <= C; z++){
                  if (a * b == minarea) break;
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
                  
                  if (RerankVector(A, B, C, totalmonthcaseu - var1, totalmonthcasew - var2, u, w) != 0) pvalue2 = RerankVector(A, B, C, totalmonthcaseu - var1, totalmonthcasew - var2, u, w);
                  double rmseans = rmse(u, w, A, B, C);
                  for (int i = x; i <= xx; ++i){
                    for (int j = y; j <= yy; ++j){
                      for (int k = z; k <= zz; ++k){
                        u[i][j][k] = u1[i][j][k];
                        w[i][j][k] = w1[i][j][k];
                      }
                    }
                  }
                  
                  float newn = nn[xx][yy][zz] - nn[x - 1][yy][zz] - nn[xx][y - 1][zz] - nn[xx][yy][z - 1] + nn[x - 1][y - 1][zz] + nn[xx][y - 1][z - 1] + nn[x - 1][yy][z - 1] - nn[x - 1][y - 1][z - 1];
                  
                    float nnn = nnozero - newn;
                  
                    if (nnn == 0 || ((pvalue2 >= 0.05) && (a * b < minarea))) {
                      rmsevalans = rmseans;
                      pvalans = pvalue2;
                      coords1 = x;
                      coords2 = y;
                      coords3 = a;
                      coords4 = b;
                      coords5 = z;
                      coords6 = c;
                      newnumber = maxn;
                      newtotal = nnn - a*b*c;
                      // log new best
                      minarea = a * b * c;
                    }
                }
              }
            }
          }
        }
    }

    cout << "P-value " << setprecision(5) << fixed << pvalans << " found at by taking out " << setprecision(0) << fixed << coords3 * coords4 * coords6 << " cells in submatrix[" << coords1 << " .. " << coords1 + coords3 - 1 << "][" << coords2 << " .. " << coords2 + coords4 - 1 << "][" << coords5 << " .. " << coords5 + coords6 - 1 << "] with min(# of 1, # of -1) = " << newnumber << " out of " << newtotal << " cells\n" << setprecision(1);
  variable1 = coords1;
  variable2 = coords2;
  variable3 = coords3;
  variable4 = coords4;
  variable5 = coords5;
  variable6 = coords6;
  variable7 = pvalue;
  variable8 = pvalans;
  variable9 = rmsevalans;

    return;
}

int main() {
    clock_t start, end;
  
    ios_base::sync_with_stdio(false);
    cin.tie(0);
    
    int A = 6, B = 9, C = 2, D = 100, month = 1;
  
    ofstream outfile;
    outfile.open("afile.txt");
    int totalmonthcase[100] = {};

    int location, age, gender;
    string date, useless, lastmonth = "1";
    ifstream fin("VancouverCovidData2.txt");
    fin >> useless >> useless >> useless >> useless >> useless;
  
    // prepare test data
    vector<vector<vector<int>>> u(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
    vector<vector<vector<int>>> w(A + 1, vector<vector<int>>(B + 1, vector <int>(C + 1, 0)));
    vector<vector<vector<vector<int>>>> v(A + 1,vector<vector<vector<int>>>(B + 1,vector<vector<int> >(C + 1,vector<int>(D + 1, 0))));
  
    while (fin >> date >> location >> gender >> age >> useless) {
        if (process(date) == lastmonth) {
            v[location][age][gender][month]++;
            totalmonthcase[month]++;
        }
        else {
            lastmonth = process(date);
            month++;
        }
    }
    for (int i = 1; totalmonthcase[i] != 0; ++i) {
        for (int j = i + 1; totalmonthcase[j] != 0; ++j) {
            for (int a = 1; a <= A; ++a) {
                for (int b = 1; b <= B; ++b) {
                  for (int c = 1; c <= C; ++c){
                    u[a][b][c] = v[a][b][c][i];
                    w[a][b][c] = v[a][b][c][j];
                  }
                }
            }
            
            start = clock();
            WilcoxonTest(u, w, A, B, C, totalmonthcase[i], totalmonthcase[j]);
            end = clock();
          
            long double time_taken = (long double)(end - start) / double(CLOCKS_PER_SEC);
            
            if (change == true) outfile << i << "\t" << j << "\t" << variable7 << "\t" << variable8 << "\t" << variable3*variable4*variable6 << "\t" << variable1 << "\t" << variable1+variable3-1 << "\t" << variable2 << "\t" << variable2+variable4-1 << "\t" << variable5 << "\t" << variable5+variable6-1 << "\t" << time_taken << '\t' << variable9 << "\n";
            if (change == false) outfile << i << "\t" << j << "\t" << "---\t\t\t\t\t\t\t\t\t" << time_taken << '\t' << variable9 << '\n';
        }
    }
  
    outfile.close();
}