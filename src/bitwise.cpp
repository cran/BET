///*
// * bitwise.cpp
// *
// *  Created on: Jul 8, 2020
// *      Author: wanzh
// */
//

#include <bitset>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <limits>
#include <fstream>
#include <sstream>
#include <map>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace std;

class BETfunction
{
public:
  int d;
  int getStats();
  double getPvalue();
  std::vector<std::string> getInteraction();
  std::vector<std::vector<std::string>> getSymmInteraction();
  std::vector<int> getSymmStats();
  BETfunction(NumericMatrix& X, int depth);

private:
  bool debug = 0;
  int numThread = 4;
  std::vector<std::vector<double>> X;
  size_t p_size;
  std::vector<std::vector<std::vector<std::string>>> symm_interaction;
  std::vector<std::vector<std::string>> out_symminter;
  std::vector<std::vector<int>> symm_stats;
  std::vector<int> out_symmstats;
  std::vector<long long> b;
  int Stats = 0;
  double pvalue;
  std::vector<std::string> interaction_str;
  std::vector<std::vector<double>> imp(NumericMatrix& X);
  std::vector<std::vector<int>> interactions(int d);
  std::vector<std::string> interaction_index(int d);
  std::vector<long long> BIDs(int d);
  std::vector<std::vector<int>> ecdf_loc(std::vector<std::vector<double>>& X, int d);
  std::map<std::vector<int>, int> groupC(std::vector<std::vector<int>>& c);
  int locate(int c, long long b, int d);
  std::vector<std::vector<std::vector<int>>> CBIDs(std::map<std::vector<int>, int>& count, size_t p, std::vector<long long>& BIDs);
  double binomial(int n, int k, double p);
  double pbinom(int n, int k, double p);
  double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);
  double fisherExact(int a, int b, int c, int d);
  void printRes(std::vector<size_t>& var ,size_t numCount, size_t n, std::vector<int>& countValues, std::vector<std::vector<std::vector<int>>>& matrix, size_t j);
};

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

vector<vector<double>> BETfunction::imp(NumericMatrix& X)
{
  int r = X.nrow(), c = X.ncol();
  vector<vector<double>> v(r, vector<double>(c));
  for (int i = 0; i < r; i++){
    for (int j = 0; j < c; j++)
      v[i][j] = X(i, j);
  }
  return v;
}

vector<vector<int>> BETfunction::interactions(int d)
{
	// all interactions on depth = d
	vector<vector<int>> res;
	vector<int> nums(d);
	for (int i = 0; i < d; i++){
		nums[i] = i + 1;
	}
	for(int num: nums){
	  int n = res.size();
	  res.push_back({num});

	  for (int i = 0; i < n; ++i){
		  auto temp = res[i];
		  temp.push_back(num);
		  res.push_back(temp);
	  }
	}
	return res;
}

vector<string> BETfunction::interaction_index(int d)
{
	vector<vector<int>> itr = interactions(d);
	vector<string> res;
	string s, str;
	string str0(d, '0');
	for (size_t i = 0; i < itr.size(); i++)
	{
		s = accumulate(itr[i].begin()+1, itr[i].end(), to_string(itr[i][0]),
				[](const string& a, int b){return a + ':' + to_string(b);});
		str = str0;
		for (size_t j = 0; j < itr[i].size(); j++){
			str.replace(itr[i][j]-1, 1, "1");
		}
		res.push_back(s + " " + str);
	}
	return res;
}

vector<long long> BETfunction::BIDs(int d)
{
	int n = (int)round(pow(2, d-1));
	vector<long long> b;
	string str1(n, '1');
	string str2(n, '0');
	//	basic configurations
	b.push_back(stoll(str2 + str1, nullptr, 2)); // A0
	for (int i = 1; i < d; i++)
	{
	//	recursively: A_i = A_{i-1} XOR (A_{i-1} << (1 << (d-i)))
		b.push_back(b[i-1] ^ (b[i-1] << (1 << (d-i-1))));
	}

	vector<vector<long long>> subs;
	//	all subsets of basic configurations
	for (long long num: b){
		int n = subs.size();
		subs.push_back({num});

		for (int i = 0; i < n; ++i){
			auto temp = subs[i];
			temp.push_back(num);
			subs.push_back(temp);
		}
	}

	double len = (int)round(pow(2, d));
	string str((int)len, '1');
	unsigned long long a0 = stoull(str, nullptr, 2);
	vector<long long> res;
	long long x;
	for (size_t i = 0; i < subs.size(); i++)
	{
	//  NOT XOR
		x = a0;
		for (size_t j = 0; j < subs[i].size(); j++){
			x = ~(x ^ subs[i][j]);
		}
		x = x & a0; // keep only 2^d digits
		res.push_back(x);
	}
	return res; // BIDs on one dimension
}

vector<vector<int>> BETfunction::ecdf_loc(vector<vector<double>>& X, int d)
{
  int n = X.size(), p = X[0].size();
  vector<vector<int>> C(n, vector<int>(p));
  double u = 0;
  vector<double> sorted(n);
  if (p == 1)
  {
    for (int i = 0; i < n; i++){
      C[i][0] = (int) ceil(X[i][0] * (int)round(pow(2, d)));
    }
  }
  else{
    for (int j = 0; j < p; j++){
      for (int i = 0; i < n; i++){
        sorted[i] = X[i][j];
      }
      // store sorted data
      sort(sorted.begin(), sorted.end(), greater<double>());
      for (int i = 0; i < n; i++){
        //	empirical cdf
        vector<double>::iterator ans = find(sorted.begin(), sorted.end(), X[i][j]);
        int index = distance(sorted.begin(), ans);
        u = n - index;
        C[i][j] = (int) ceil(u/n * (int)round(pow(2, d)));
      }
    }
  }
  return C; // Cij
}

map<vector<int>, int> BETfunction::groupC(vector<vector<int>>& c)
{
	map<vector<int>, int> count;
	for (auto& v: c){
		count[v]++;
	}
	return count;
}

int BETfunction::locate(int c, long long b, int d)
{
	// the c th bits of b (BID)
	int move = (int)round(pow(2, d)) - c;
	int res = (b >> move) & 1;
	return res;
}

vector<vector<vector<int>>> BETfunction::CBIDs(map<vector<int>, int>& count, size_t p, vector<long long>& BIDs)
{
  size_t numCount = count.size(), numBIDs = BIDs.size();
  // store all Cij * BIDs in to a matrix: p * #BIDs * #Cijs
  vector<vector<vector<int>>> res(p, vector<vector<int>>(numBIDs, vector<int>(numCount)));

  for (int i = 0; i < (int)p; i++){
    size_t j = 0;
    vector<vector<int>> temp(numBIDs, vector<int>(numCount));
    // go over all BIDs:
    for (map<vector<int>,int>::iterator it=count.begin(); it!=count.end()&&(j < numCount); ++it, j++){
      // go over all p dimensions:
      for (size_t k = 0; k < numBIDs; k++){
        temp[k][j] = locate(it->first[i], BIDs[k], d);
        //				res[i][k][j] = locate(it->first[i], BIDs[k], d);
      }
    }
    res[i] = temp;
    //		cout << "thread " + to_string(omp_get_thread_num()) + "\n";
  }

  return res;
}


double BETfunction::binomial(int n, int k, double p)
{
	// binomial distribution with dp
	if (n < 0 || k < 0)
		return 0;
	// store previous results
	vector<vector<double> > dp(n + 1, vector<double>(k+1));
	// marginal condition
	dp[0][0] = 1;
	for (int i = 1; i < n+1; ++i){
		dp[i][0] = (1 - p) * dp[i-1][0];
	}
	for (int j = 1; j < k+1; ++j){
		dp[0][j] = 0;
	}
	for (int i = 1; i < n+1; ++i){
		for (int j = 1; j < k+1; ++j){
			// recursion
			dp[i][j] = (1 - p) * dp[i-1][j] + p * dp[i-1][j-1];
		}
	}
	return dp[n][k];
}

double BETfunction::pbinom(int n, int k, double p)
{
	// cumulative: P(Binom(n, p) <= k)
	if (n < 0 || k < 0)
		return 0;
	double pb = 0;
	for (int i = 0; i <= k; i++){
		pb += binomial(n, i, p);
	}
	return pb;
}

double BETfunction::logHypergeometricProb(double* logFacs, int a, int b, int c, int d)
{
  return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double BETfunction::fisherExact(int a, int b, int c, int d)
{
  int n = a + b + c + d;
  double* logFacs = new double[n+1];

  for(int i=1; i < n+1; ++i) {
    logFacs[i] = logFacs[i-1] + log((double)i); //store factorization
  }
  double logpCutoff = logHypergeometricProb(logFacs,a,b,c,d); // *** logFacs added
  double pFraction = 0;
  for(int x=0; x <= n; ++x) {
    if ( a+b-x >= 0 && a+c-x >= 0 && d-a+x >=0 ) {
      double l = logHypergeometricProb(logFacs, x,a+b-x,a+c-x,d-a+x);
      if ( l <= logpCutoff )
        pFraction += exp(l - logpCutoff);
    }
  }
  double logpValue = logpCutoff + log(pFraction);
  double hp = exp(logpCutoff)/2;
  delete [] logFacs;
  return exp(logpValue) - hp;
}

void BETfunction::printRes(vector<size_t>& var, size_t numCount, size_t n, vector<int>& countValues, vector<vector<vector<int>>>& matrix, size_t j)
{
  // numCount = #different locations
  // j: j'th combination of Xi

  // #rows of current interaction instead of total # of rows
  size_t p = matrix.size();
  vector<size_t> indices(p);

  vector<vector<int>> choose(p, vector<int>(numCount));
  size_t i = 0;
  size_t k = 0;

  size_t total = pow(pow(2, d) - 1, p);
  vector<int> symmStats(total);
  string str0(d, '0');
  //initial interaction: 0..0
  vector<vector<string>> symmInter(p_size, vector<string>(total, str0));
  vector<string> inter = interaction_index(d);

  while (i < p){
    // all combination of BIDs (p-dimensions)
    for (size_t j = 0; j < p; j++){
      symmInter[var[j]][k] = inter[indices[j]];

      // j'th variable:
      choose[j] = matrix[j][indices[j]];
    }

    int loc0 = 1, sum = 0;

    for (size_t k = 0; k < numCount; k++){
      for (size_t l = 0; l < p; l++){
        // NOT(XOR) p digits
        loc0 = (~(loc0 ^ choose[l][k])) & 1;
      }
      sum += loc0 * countValues[k];
      loc0 = 1;
    }

    //store symmetry statistics
    symmStats[k] = sum - (int)n/2;

    i = 0;
    while (i < matrix.size() && ++indices[i] == matrix[i].size())
      indices[i++] = 0;

    k++;
  }

  // fill in symmetry statistics and interaction
  symm_interaction[j] = symmInter;
  symm_stats[j] = symmStats;
}


int BETfunction::getStats()
{
	return Stats;
}

double BETfunction::getPvalue()
{
  return pvalue;
}

vector<string> BETfunction::getInteraction()
{
	return interaction_str;
}

vector<int> BETfunction::getSymmStats()
{
  return out_symmstats;
}

vector<vector<string>> BETfunction::getSymmInteraction()
{
  return out_symminter;
}

BETfunction:: BETfunction(NumericMatrix& X_R, int depth)
{
	d = depth;
  X = imp(X_R);
	size_t n = X.size();
	size_t p = X[0].size();
	p_size = X[0].size();

	interaction_str.resize(p);

	vector<vector<int>> c;
	// the empirical cdf of X and multiplxied by 2^d

	  c = ecdf_loc(X, d);
	  // create a hash table that count for same column
	  map<vector<int>, int> count = groupC(c);
	  size_t numCount = count.size();

	  // number of each location
	  vector<int> countValues;
	  for (map<vector<int>,int>::iterator it=count.begin(); it!=count.end(); ++it){
	    countValues.push_back(it->second);
	  }

	  b = BIDs(d); // all BIDs on one-dimension

	  // all locations operate with all BIDs
	  vector<vector<vector<int>>> matrix0 = CBIDs(count, p, b);

	  vector<vector<size_t>> res;

	  vector<size_t> nums(p);
	  for (size_t i = 0; i < p; i++){
	    nums[i] = i;
	  }
	  // all combinations of variables
	  for(size_t num: nums){
	    size_t n = res.size();
	    res.push_back({num});

	    for (size_t i = 0; i < n; ++i){
	      auto temp = res[i];
	      temp.push_back(num);
	      res.push_back(temp);
	    }
	  }

	  // all combinations of variables
	  if (p > 2){
	    vector<vector<size_t>>::iterator it;
	    for(it = res.begin(); it != res.end();){
	      if ((*it).size() == 1)
	        it = res.erase(it);
	      else
	        ++it;
	    }

	    for (size_t i = 0; i < p; i++){
	      interaction_str[i].resize(res.size());
	    }


	    vector<string> inter = interaction_index(d);

	    // steps of multi-processing
	    size_t total = res.size();
	    size_t remainder = total % numThread;
	    size_t step = total / numThread;
	    vector<size_t> interval(numThread, step);
	    for (size_t i = 0; i < remainder; i++){
	      interval[i] += 1;
	    }

	    vector<size_t> thread(numThread + 1);
	    for (int i = 1; i < numThread + 1; i++){
	      thread[i] = thread[i - 1] + interval[i - 1];
	    }

	    size_t comb = pow(2, p) - p -1;

	    symm_interaction.resize(comb);
	    symm_stats.resize(comb);

#ifdef _OPENMP
omp_set_num_threads(numThread);
#pragma omp parallel for
#endif

	    for (int thr = 0; thr < numThread; thr++){
	      //				double wtime = omp_get_wtime();

	      vector<vector<vector<int>>> matrix;

	      for (size_t j = thread[thr]; j < thread[thr+1]; j++){
	        for (size_t i = 0; i < res[j].size(); i++){
	          matrix.push_back(matrix0[res[j][i]]);
	          // if (debug)
	          //   cout << "X" << res[j][i]+1 << " ";
	        }

	        // if (debug)
	        //   cout << endl;

	        printRes(res[j], numCount, n, countValues, matrix, j);

	        // if (debug)
	        //   cout << "===========================" << endl;
	        matrix.clear();
	      }
	      //				wtime = omp_get_wtime() - wtime;
	      //				cout << "thread " + to_string(omp_get_thread_num()) + " time: " + to_string(wtime) +"\n";
	    }
	    for (size_t i = 0; i < comb; i++){
	      out_symminter.resize(p_size);
	      for (size_t j = 0; j < p_size; j++){
	        out_symminter[j].insert(out_symminter[j].end(), symm_interaction[i][j].begin(), symm_interaction[i][j].end());
	      }

	      out_symmstats.insert(out_symmstats.end(), symm_stats[i].begin(), symm_stats[i].end());
	    }

	  }
	  else{
	    // p < 2: do not adopt multi-processing
	    // if (debug)
	    //   cout << "X" << 1 << endl;
	    if (p == 2)
	      res = {{0, 1}};

	    interaction_str.resize(p);
	    for (size_t i = 0; i < p; i++){
	      interaction_str[i].resize(1);
	    }

	    symm_interaction.resize(1);
	    symm_stats.resize(1);

	    vector<vector<vector<int>>> matrix0 = CBIDs(count, p, b);

	    vector<string> inter = interaction_index(d);

	    printRes(res[0], numCount, c.size(), countValues, matrix0, 0);

	    out_symminter.resize(p_size);
	    for (size_t j = 0; j < p_size; j++){
	      out_symminter[j].insert(out_symminter[j].end(), symm_interaction[0][j].begin(), symm_interaction[0][j].end());
	    }

	    out_symmstats.insert(out_symmstats.end(), symm_stats[0].begin(), symm_stats[0].end());

	  }
	  // find the extreme stats
	  vector<int>::iterator findMax = max_element(out_symmstats.begin(), out_symmstats.end(), abs_compare);
	  size_t max_idx = distance(out_symmstats.begin(), findMax);
	  Stats = out_symmstats[max_idx];
	  pvalue = ((int)round(pow(2, (int)p_size*d)) - (int)p_size*((int)round(pow(2, d)) - 1) - 1) * fisherExact((int)n/4 + abs(Stats)/2, (int)n/4 - abs(Stats)/2, (int)n/4 - abs(Stats)/2, (int)n/4 + abs(Stats)/2);
	  for (size_t i = 0; i < p_size; i++){
	    interaction_str[i] = out_symminter[i][max_idx];
	  }

}

//[[Rcpp::export]]
List symmCpp(NumericMatrix& X, int d)
{
  BETfunction bet(X, d);
  size_t p = X.ncol();

  DataFrame symm = DataFrame::create(
    Named("SymmetryStatistics") = bet.getSymmStats()
  );

  for (size_t i = p; i > 0; i--){
    symm.push_front(bet.getSymmInteraction()[i-1], "X" + to_string(i));
  }

  List L = List::create(Named("SymmetryStatistics") = DataFrame(symm));

  return L;

}

//[[Rcpp::export]]
List BETCpp(NumericMatrix& X, int d)
{
  BETfunction bet(X, d);
  int n = X.nrow();
  int white = bet.getStats() + n/2;
  double zstats = (double)((1.0*white)/(1.0*n) - 0.5)/sqrt(0.25/n);

  DataFrame df = DataFrame::create();

  size_t p = X.ncol();
  for (size_t i = p; i > 0; i--){
    CharacterVector vi = {bet.getInteraction()[i-1]};
    df.push_front(bet.getInteraction()[i-1], "X" + to_string(i));
  }

  List L = List::create(Named("Interaction") = DataFrame(df), Named("Extreme.Asymmetry") = bet.getStats(), Named("p.value") = bet.getPvalue(), Named("z.statistic") = zstats);

  return L;

}
