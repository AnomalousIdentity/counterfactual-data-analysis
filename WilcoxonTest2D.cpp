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

double R[1000] = { 0 };

bool sortbysec(const pair<long double, long double>& a, const pair<long double, long double>& b)
{
    return (a.first < b.first);
}

long double normalCDF(long double value)
{
   return 0.5*erfc(-value * M_SQRT1_2);
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

void WilcoxonTest(vector<int> u, vector<int> w, int A, int totalmonthcaseu, int totalmonthcasew)
{
    ios::sync_with_stdio(false);
    cin.tie(0);
  
    change = true;
  
    int n = A;
    int i, j, k, a, b, z;
    float s, t;

    double list[n];
    int negs[200] = {0};
    int zerocnt = 0;

    for (int i = 1; i <= A; i++) {
      if ((long double)u[i] / (long double)totalmonthcaseu < (long double)w[i] / (long double)totalmonthcasew) {
        list[i - 1] = (double)(w[i]) / totalmonthcasew - (double)(u[i]) / (double)totalmonthcaseu;
        negs[i - 1] = 1;
      }
      else if ((long double)u[i] / (long double)totalmonthcaseu > (long double)w[i] / (long double)totalmonthcasew) {
        list[i - 1] = (double)(u[i]) / (double)totalmonthcaseu - (double)(w[i]) / totalmonthcasew;
      }
      else {
        list[i - 1] = 0;
        zerocnt++;
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

    vector<float> v(A + 1, 0);
    for (int i = 1; i <= A; i++) {
        (negs[i - 1] == 1) ? (v[i] = -R[i - 1]) : (v[i] = R[i - 1]);
    }
    
    // prepare 2D prefix sum, sumpos[a][b] = total of positives in v[i][j] where i <= a and j <= b (top left corner with a rows and b columns)
    int ss, tt;
    vector<float> sumpos(A + 1, 0);
    vector<float> sumneg(A + 1, 0);
    vector<float> numpos(A + 1, 0);
    vector<float> numneg(A + 1, 0);

    //number of units without including 0
    vector<float> nn(A + 1, 0);
    for (i = 1; i <= A; i++) {
      if (v[i] > 0) sumpos[i] = sumpos[i-1] + v[i];
      if (v[i] > 0) numpos[i] = numpos[i-1] + 1;
      if (v[i] < 0) sumneg[i] = sumneg[i-1] - v[i];
      if (v[i] < 0) numneg[i] = numneg[i-1] + 1;
      if (v[i] != 0) nn[i] = nn[i-1] + 1;
    }

    // total # of positives and negatives in the original matrix
    float totalpos = 0, totalnegs = 0;
    queue<float> posabove, negsabove;
    for (i = 1; i <= A; i++) {
      if (v[i] > 0) {
        totalpos += v[i];
        posabove.push(v[i]);
      }
      else if (v[i] < 0) {
        totalnegs -= v[i];
        negsabove.push(-v[i]);
      }
    }

    //change in rank
    vector<float> rankch(A + 1, 0);
    for (i = 1; i <= A; i++) {
      if (v[i] > 0 && max(totalpos, totalnegs) < 0) {
        for (k = 0; k < negsabove.size(); k++) {
          if (negsabove.front() > v[i]) {
            rankch[i] ++;
          }
          else if (negsabove.front() == v[i]) {
            rankch[i] += 0.5;
          }
          negsabove.push(negsabove.front());
          negsabove.pop();
        }
      }
      else if (v[i] < 0 && max(totalpos, totalnegs) > 0) {
        for (k = 0; k < posabove.size(); k++) {
          if (posabove.front() > abs(v[i])) {
            rankch[i] ++;
          }
          else if (posabove.front() == abs(v[i])) {
            rankch[i] += 0.5;
          }
      posabove.push(posabove.front());
      posabove.pop();
        }
      }
    }

    //rankchange for top left
    vector<float> rankc(A + 1, 0);
    
    for (i = 1; i <= A; i++) {
      if (rankch[i] > 0) rankc[i] = rankc[i-1] + rankch[i];
    }
  
    // get z scores
    int nnozero = n - zerocnt;
    long double zscore = (max(totalpos, totalnegs) - nnozero * (nnozero + 1) / 4) / sqrt(nnozero * (nnozero + 1) * (2 * nnozero + 1) / 24);
    long double pvalue = 2*(1-normalCDF(zscore));

    // cout << "--------------------------------------" << '\n';
    // for (i = 1; i <= A; i++){
    //   cout << v[i] << ' ';
    // }
    // cout << '\n';

    if (pvalue >= 0.05) {
      // cout << "no change needed, " << max(totalpos, totalnegs) << " total out of " << nnozero << " cells." << '\n';
      cout << "no change needed " << ' '<< zscore << ' ' << pvalue << "\n";
      change = false;
      return;
    }
  
    // track the best answer

    float minarea = n+1;

    long double zvalans = 0, pvalans = 0, coords1 = 0, coords2 = 0, newnumber = 0, newtotal = 0;

    // enumerate all window size a x b

    for (a = 1; a <= A; a++) {
      if (a > minarea) continue;
        // enumerate all topleft cells
        for (int x = 1; x + a - 1 <= A; x++) {
          if (a == minarea) break;
          int xx = x + a - 1;
        // now calculate the # of 1s in the submatrix [x..xx], [y..yy]
          long double c1, c_1, new1, new_1, newn, maxn, rankchange;
          int nnn;
          c1 = numpos[xx] - numpos[x - 1];
          c_1 = numneg[xx] - numneg[x - 1];
          newn = nn[xx] - nn[x - 1];
          rankchange = rankc[xx] - rankc[x - 1];

          new1 = totalpos;
          new_1 = totalnegs;
          maxn = max(new_1, new1);
          nnn = nnozero - newn;
          (new_1 > new1) ? newn -= c1 : newn -= c_1;

          for (i = 0; i < newn; ++i) {
            maxn -= (nnozero - i);
          }
          maxn += rankchange;

          // calculate z
          long double zvalue = (maxn - (nnn * (nnn + 1) / 4)) / sqrt(nnn * (nnn + 1) * (2 * nnn + 1) / 24);
          long double pvalue = 2*(1-normalCDF(zvalue));

          if (nnn == 0 || (zvalue <= 1.96) && (a * b < minarea)) {
            if (nnn == 0) zvalue = 0;
            pvalans = pvalue;
            coords1 = x;
            coords2 = a;
            newnumber = maxn;
            newtotal = nnn;
            // log new best
            minarea = a;
          }
        }
    }

    cout << "P-value " << setprecision(5) << fixed << pvalans << " found at by taking out " << setprecision(0) << fixed << coords2 << " cells in submatrix[" << coords1 << " .. " << coords1 + coords2 - 1 << "] with min(# of 1, # of -1) = " << newnumber << " out of " << newtotal << " cells\n" << setprecision(1);
  variable1 = coords1;
  variable2 = coords2;
  variable3 = pvalue;
  variable4 = pvalans;

    return;
}

int main() {
    clock_t start, end;
    start = clock();
  
    ios_base::sync_with_stdio(false);
    cin.tie(0);

    int A = 9, B = 100, month = 1;

   ofstream outfile;
   outfile.open("afile.txt");
  
    int totalmonthcase[100] = {};
    int location, age, gender;
    string date, useless, lastmonth = "1";
    ifstream fin("VancouverCovidData2.txt");
    fin >> useless >> useless >> useless >> useless >> useless;
    // prepare test data
    
    vector<int> u(A + 1, 0);
    vector<int> w(A + 1, 0);
    vector<vector<int>> v(A + 1, vector<int>(B + 1, 0));

    while (fin >> date >> location >> gender >> age >> useless) {
        if (process(date) == lastmonth) {
            v[age][month]++;
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
                  u[a] = v[a][i];
                  w[a] = v[a][j];
                }
            }
            start = clock();
            WilcoxonTest(u, w, A, totalmonthcase[i], totalmonthcase[j]);
            end = clock();
            long double time_taken = (long double)(end - start) / double(CLOCKS_PER_SEC);
            
            if (change == true) outfile << i << "\t" << j << "\t" << variable3 << "\t" << variable4 << "\t" << variable2 << "\t" << variable1 << "\t" << variable1+variable2-1 << "\t" << time_taken << "\n";
            if (change == false) outfile << i << "\t" << j << "\t" << "---\t\t\t\t\t" << time_taken << '\n';
        }
    }
    outfile.close();
}