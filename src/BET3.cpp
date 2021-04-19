//============================================================================
// Name        : BET2nd.cpp
// Author      : wanzhang
// Version     :
// Copyright   : Your copyright notice
// Description : BET 2nd version in C++, Ansi-style
//============================================================================

#include "BET.h"
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
#include <cassert>
#include <iterator>
#include <random>
#include <chrono>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace N;

inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// vector<vector<double> > BETfunction::imp(NumericMatrix& X)
// {
//   int r = X.nrow(), c = X.ncol();
//   vector<vector<double>> v(r, vector<double>(c));
//   for (int i = 0; i < r; i++){
//     for (int j = 0; j < c; j++)
//       v[i][j] = X(i, j);
//   }
//   return v;
// }

vector<vector<int> > BETfunction::interactions()
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

vector<vector<int>> BETfunction::interaction_mat()
{
  // matrix interaction for BEAST
  vector<vector<int>> itr = interactions();
  size_t len = itr.size();
  vector<vector<int>> res(len + 1, vector<int>(d));

  for(size_t i = 0; i < len; i++){
    for(size_t j = 0; j < itr[i].size(); j++){
      res[i+1][itr[i][j]-1] = 1;
    }
  }

  return res;
}

vector<string> BETfunction::interaction_index()
{
	// interactions in string: 1-d and binary
	vector<vector<int>> itr = interactions();
	vector<string> res;
	string str;
	string str0(d, '0');
	res.push_back(str0);
	for (size_t i = 0; i < itr.size(); i++)
	{
		// binary
		str = str0;
		for (size_t j = 0; j < itr[i].size(); j++){
			str.replace(itr[i][j]-1, 1, "1");
		}
		res.push_back(str);
	}
	return res;
}

vector<long long> BETfunction::BIDs()
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

vector<vector<int>> BETfunction::ecdf_loc(vector<vector<double>>& X)
{
  int n = X.size(), p = X[0].size();
  vector<vector<int>> C(n, vector<int>(p));
  double u = 0;
  vector<double> sorted(n);
  if (p == 1)
  {
    for (int i = 0; i < n; i++){
      // throw error in R
      //			if (X[i][0] > 1 or X[i][0] < 0){
      //				// test whether unif[0,1]
      //				throw "Data out of range [0, 1]!";
      //			}
      C[i][0] = (int) ceil(X[i][0] * (int)round(pow(2, d)));
    }
  }
  else{
    if (unifMargin)
      // already uniformed
      for (int j = 0; j < p; j++){
        for (int i = 0; i < n; i++){
          C[i][j] = (int)X[i][j];
        }
      }
      else
      {
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

int BETfunction::locate(int c, long long b)
{
  // the c th bits of b (BID)
  int move = (int)round(pow(2, d)) - c;
  int res = (b >> move) & 1;
  return res;
}

vector<vector<vector<int>>> BETfunction::CBIDs(map<vector<int>, int>& count)
{
  size_t numCount = count.size(), numBIDs = (int)round(pow(2, d));
  // store all Cij * BIDs in to a matrix: p * #BIDs * #Cijs
  vector<vector<vector<int>>> res(p, vector<vector<int>>(numBIDs, vector<int>(numCount)));

  for (int i = 0; i < (int)p; i++){
    size_t j = 0;
    vector<vector<int>> temp(numBIDs, vector<int>(numCount, 1));
    // go over all BIDs:
    for (map<vector<int>, int>::iterator it=count.begin(); it!=count.end()&&(j < numCount); ++it, j++){
      // go over all p dimensions:
      for (size_t k = 0; k < numBIDs - 1; k++){
        temp[k+1][j] = locate(it->first[i], bid[k]);
      }
    }
    res[i] = temp;
  }
  return res;
}

vector<vector<size_t>> BETfunction::allComb()
{
  // # of BIDs
  size_t b = (int)round(pow(2, d));

  // # of combinations = 2^d^p
  size_t rows = (int)round(pow(b, p));

  // store all combinations
  vector<vector<size_t>> res(rows, vector<size_t>(p));

  vector<size_t> indices(p);
  vector<size_t> choose(p);

  // p columns, each column from 0 to 2^d-1
  vector<vector<size_t>> comb(p, vector<size_t>(b));
  for (size_t i = 0; i < p; i++){
    for (size_t j = 0; j < b; j++){
      comb[i][j] = j;
    }
  }

  size_t i = 0;
  size_t k = 0;
  while (i < p){
    for (size_t j = 0; j < p; j++){
      choose[j] = comb[j][indices[j]];
    }

    // store into res
    res[k] = choose;

    i = 0;
    while (i < comb.size() && ++indices[i] == comb[i].size())
      indices[i++] = 0;

    k++;
  }
  return res;
}

static bool abs_compare(int a, int b)
{
  return (std::abs(a) < std::abs(b));
}

bool BETfunction::isIndex(vector<size_t>& idx)
{
  bool res = 1;
  if(!testUnif){
    for(size_t i = 0; i < indepIndex.size(); i++){
      bool nonZero = 0;
      for(size_t j = 0; j < indepIndex[i].size(); j++){
        nonZero = nonZero || (idx[indepIndex[i][j]-1] != 0);
      }
      res = res && nonZero;
    }
  }

  return res;
}

vector<int> BETfunction::symmstats(vector<int>& countValues, vector<vector<vector<int>>>& matrix, size_t sampleSize)
{
  size_t total = (int)round(pow(2, p*d));
  vector<int> symmstats(total, 0);
  size_t numCount = countValues.size();
  // vector<vector<int>> inter_idx(total);

#ifdef _OPENMP
  omp_set_num_threads(numThread);
#pragma omp parallel for
#endif

  for (size_t th = 1; th < numThread + 1; th++){
    //		double wtime = omp_get_wtime();
    for (size_t i = thread[th - 1]; i < thread[th]; i++){
      string bi = "";
      // vector<int> iidx(p);
      for (size_t v = 0; v < p; v++){
        // store interaction
        out_symminter[v][i] = inter[allidx[i][v]];
        bi += inter[allidx[i][v]];
        // iidx[v] = allidx[i][v];
      }
      binary_inter[i] = bi;
      // inter_idx[i] = iidx;
      if ( ( (count(allidx[i].begin(), allidx[i].end(), 0) <= (int)(p - 2)) && (testUnif ||isIndex(allidx[i]) ) ) || (p == 1 && i > 0) ){
        int loc0 = 1, sum = 0;
        for (size_t j = 0; j < numCount; j++){
          for (size_t k = 0; k < p; k++){
            loc0 = (~(loc0 ^ matrix[k][allidx[i][k]][j])) & 1;
          }
          sum += loc0 * countValues[j];
          loc0 = 1;
        }
        //store symmetry statistics
        symmstats[i] = 2 * sum - (int)sampleSize;
      }
    }
  }

  // find the extreme stat
  vector<int>::iterator findMax = max_element(symmstats.begin(), symmstats.end(), abs_compare);
  size_t max_idx = distance(symmstats.begin(), findMax);

  // max interaction index
  // vector<int> max_interaction = inter_idx[max_idx];
  string max_interaction = binary_inter[max_idx];

  // count most frequent max interaction
  countInteraction[max_interaction]++;

  return symmstats;
}

vector<size_t> BETfunction::unreplaceShuffle(size_t size, size_t max_size)
{
  assert(size <= max_size);
  vector<size_t> b(size);
  vector<size_t> t(max_size);

  for(size_t i = 0; i < max_size; i++){
    t[i] = i;
  }

  // shuffle 0~max_size
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // shuffle(t.begin(), t.end(), gen);
  random_shuffle(t.begin(), t.end(), randWrapper);

  // get first size entry
  for(size_t i = 0; i < size; i++){
    b[i] = t[i];
  }

  return b;
}

vector<double> BETfunction::subsample(size_t m, size_t B)
{
  size_t total = (int)round(pow(2, p*d));
  vector<double> res(total);

  // bool debug = 1;

  // subsample B times
  for(size_t b = 0; b < B; b++){
    vector<size_t> idx = unreplaceShuffle(m, n);

    // if(debug){
    //   for(size_t i = 0; i < m; i++){
    //     cout << idx[i] << " ";
    //   }
    //   cout << " The " << b+1 << " time index " << endl;
    // }

    //	new data
    vector<vector<double> > X_m(m, vector<double>(p));
    for(size_t i = 0; i < m; i++){
      X_m[i] = X[idx[i]];
    }

    vector<vector<int>> c = ecdf_loc(X_m);

    // create a map that count for observations in the same location
    map<vector<int>, int> mapC = groupC(c);

    // number of each location
    vector<int> countValues;
    for (map<vector<int>, int>::iterator it=mapC.begin(); it!=mapC.end(); ++it){
      countValues.push_back(it->second);
    }

    // 3-dimension matrix: Xij ~ BIDs
    vector<vector<vector<int>>> matrix = CBIDs(mapC);

    // binary expansion on X_m
    vector<int> S = symmstats(countValues, matrix, m);

    // if(debug){
    //   for(size_t i = 0; i < S.size(); i++){
    //     cout << S[i] << " ";
    //   }
    //   cout << " The " << b+1 << " time sample " << endl;
    // }


    // sum S1 to SB
    for(size_t i = 0; i < total; i++){
      res[i] += (double)S[i]/(double)m;
    }
  }

  // take mean
  for(size_t i = 0; i < total; i++){
    res[i] = res[i]/(double)B;
  }

  return res;
}

vector<double> BETfunction::softthreshold(vector<double>& S, double lambda)
{
  size_t length = S.size();
  vector<double> res(length);
  for(size_t i = 0; i < length; i++){
    if(S[i] >= 0){
      if((double)S[i] - lambda > 0){
        res[i] = (double)S[i] - lambda;
      }
    }else{
      if((double)S[i] + lambda < 0){
        res[i] = (double)S[i] + lambda;
      }
    }
  }

  return res;
}

double BETfunction::logHypergeometricProb(double* logFacs, int a, int b, int c, int d)
{
  return logFacs[a+b] + logFacs[c+d] + logFacs[a+c] + logFacs[b+d] - logFacs[a] - logFacs[b] - logFacs[c] - logFacs[d] - logFacs[a+b+c+d];
}

double BETfunction::fisherExact(int a, int b, int c, int d)
{
  int n = a + b + c + d;
  //	vector<double> logFacs(n+1);
  double* logFacs = new double[n+1] ();

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

double BETfunction::binomial(int n, int k, double p)
{
  // binomial distribution with dp
  // store previous results
  double **ret = new double*[n + 1];
  for (int i = 0; i < n + 1; ++i)
    ret[i] = new double[k + 1] ();

  ret[0][0] = 1.0;
  for (int i = 1; i < n + 1; ++i)
    ret[i][0] = (1.0 - p) * ret[i - 1][0];
  for (int j = 1; j < k + 1; ++j)
    ret[0][j] = 0.0;
  for (int i = 1; i < n + 1; ++i){
    for (int j = 1; j < k + 1; ++j){
      ret[i][j] = (1.0 - p) * ret[i - 1][j] + p * ret[i - 1][j - 1];
    }
  }

  for (int i = 0; i < n + 1; ++i)
    delete [] ret[i];
  delete [] ret;

  return ret[n][k];
}

double BETfunction::pbinom(int n, int k, double p)
{
  // cumulative: P(Binom(n, p) >= k)
  if (n < 0 || k < 0)
    return 0;
  // mid p value correction
  double pb = binomial(n, k, p) / 2;
  for (int i = k + 1; i <= n; i++){
    pb += binomial(n, i, p);
  }
  return pb;
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

vector<string> BETfunction::getBinary()
{
  return binary_inter;
}

double BETfunction::getBeastStat()
{
  return BeastStat;
}

// vector<vector<int>> BETfunction::getBeastInteraction()
// {
//   return BeastInter;
// }

string BETfunction::getBeastInteraction()
{
  return freqInter;
}


void BETfunction::Beast(size_t m, size_t B, double lambda, bool test_uniformity, bool test_independence, vector<vector<size_t>>& independence_index)
{
  // bool debug = 1;

  testUnif = test_uniformity;
  testIndep = test_independence;

  if(testIndep){
    indepIndex = independence_index;
  }else{
    indepIndex = {{}};
    for(size_t i = 0; i < p; i++){
      indepIndex[0].push_back(i);
    }
  }

  // first subsample
  vector<double> S1 = subsample(m, B);

  // if(debug){
  //   for(size_t i = 0; i < S1.size(); i++){
  //     cout << S1[i] << " ";
  //   }
  //   cout << "  S1" <<endl;
  // }

  // soft-thresholding lambda = sqrt(2log(2^pd)/n)
  //		lambda = sqrt(2 * log(pow(2, p*d)) / n) / 2;

  //soft-threshold
  vector<double> T1 = softthreshold(S1, lambda);
  // if(debug){
  //   for(size_t i = 0; i < S1.size(); i++){
  //     cout << T1[i] << " ";
  //   }
  //   cout << "  T1" <<endl;
  // }

  // second subsample
  vector<double> S2 = subsample(m, B);

  // if(debug){
  //   for(size_t i = 0; i < S2.size(); i++){
  //     cout << S2[i] << " ";
  //   }
  //   cout << "  S2" <<endl;
  // }

  //soft-threshold
  vector<double> T2 = softthreshold(S2, lambda);
  // if(debug){
  //   for(size_t i = 0; i < S1.size(); i++){
  //     cout << T2[i] << " ";
  //   }
  //   cout << "  T2" <<endl;
  // }

  for(size_t i = 0; i < S1.size(); i++){
    BeastStat += T1[i]*T2[i];
  }

  // BeastStat = pow(BeastStat, 2);

  using pair_type = decltype(countInteraction)::value_type;
  auto mostFreq = max_element(countInteraction.begin(), countInteraction.end(), [] (const pair_type & p1, const pair_type & p2) {
    return p1.second < p2.second;
  });



  // if(debug){
  //   cout <<"interaction : " << mostFreq->first << endl;
  // }

  freqInter = mostFreq->first;

  // for(size_t i = 0; i < p; i++){
  //   BeastInter.push_back(inter_mat[freqInter[i]]);
  // }

}


BETfunction:: BETfunction(vector<vector<double>>& X_R, int depth, bool unif, bool asymptotic, bool test_uniformity, bool test_independence, vector<vector<size_t>>& independence_index)
{
	unifMargin = unif;
  asympt = asymptotic;
  X = X_R;
	d = depth;
	n = X.size();
	p = X[0].size();

	testUnif = test_uniformity;
	testIndep = test_independence;

	if(testIndep){
	  indepIndex = independence_index;
	}else{
	  indepIndex = {{}};
	  for(size_t i = 0; i < p; i++){
	    indepIndex[0].push_back(i);
	  }
	}

	bid = BIDs();

	inter = interaction_index();
	inter_mat = interaction_mat();

	// all variables go over all bids:
	allidx = allComb();
	size_t total = (int)round(pow(2, p*d));


	// store symmetry statistics
	out_symmstats.resize(total, 0);

	// store interactions
	vector<string> temp(total);
	out_symminter.resize(p, temp);

	// store binary interactions
	binary_inter.resize(total);

	// start multiprocessing
#ifdef _OPENMP
	if(p >= 2)
	  numThread = 4;
#endif

	// loops: 2^pd
	// steps of multi-processing
	size_t remainder = total % numThread;
	size_t step = total / numThread;

	// # of assignments each thread
	vector<size_t> interval(numThread, step);
	for (size_t i = 0; i < remainder; i++){
	  interval[i] += 1;
	}

	thread.resize(numThread + 1);
	for (size_t i = 1; i < numThread + 1; i++){
	  thread[i] = thread[i - 1] + interval[i - 1];
	}


	// empirical cdf
	vector<vector<int>> c = ecdf_loc(X);

	// create a map that count for observations in the same location
	map<vector<int>, int> mapC = groupC(c);
	size_t numCount = mapC.size();

	// number of each location
	vector<int> countValues;
	for (map<vector<int>, int>::iterator it=mapC.begin(); it!=mapC.end(); ++it){
	  countValues.push_back(it->second);
	}

	// 3-dimension matrix: Xij ~ BIDs
	vector<vector<vector<int>>> matrix = CBIDs(mapC);

	// binary expansion on X
	out_symmstats = symmstats(countValues, matrix, n);

	// find the extreme stat
	vector<int>::iterator findMax = max_element(out_symmstats.begin(), out_symmstats.end(), abs_compare);
	size_t max_idx = distance(out_symmstats.begin(), findMax);

	// maxStat
	Stats = out_symmstats[max_idx];
	// max interaction
	interaction_str.resize(p);

	// store marginal
	int Sa = 0, Sb = 0;

	for (size_t i = 0; i < p; i++){
	  vector<string>::iterator ans = find(inter.begin(), inter.end(), out_symminter[i][max_idx]);
	  int index = distance(inter.begin(), ans);
	  interaction_str[i] = inter[index];
	  if (p == 2 && i == 0){
	    // count marginal
	    for (size_t j = 0; j < numCount; j++){
	      Sa += matrix[i][index][j] * countValues[j];
	    }
	  }
	  if (p == 2 && i == 1){
	    // count marginal
	    for (size_t j = 0; j < numCount; j++){
	      Sb += matrix[i][index][j] * countValues[j];
	    }
	  }
	}

	if (p > 2 || asympt){
	  // asymptotic p value: normal N (0, n)
	  double z = abs(Stats)/sqrt((double)n);
	  if (p == 1){
	    pvalue = (1-0.5 * erfc(-z * sqrt(0.5))) * ((int)round(pow(2, (int)d)) - 1) * 2;
	  }else{
	    pvalue = (1-0.5 * erfc(-z * sqrt(0.5))) * ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * 2;
	  }
	}else if (unifMargin){
	  // p value: if unif >> binomial
	  // bonferroni
	  pvalue = ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * pbinom(n, (abs(Stats) + n) / 2, 0.5);
	}else if (p == 1){
	  // p value: p = 1 >> binomial
	  // bonferroni
	    pvalue = ((int)round(pow(2, (int)d)) - 1) * pbinom(n, (abs(Stats) + n) / 2, 0.5);
	}else if (p == 2){
	  // p = 2: fisher exact
	  int n11 = ((n + Stats) / 2 + Sb - (n - Sa))/2;
	  int n00 = (n + Stats) / 2 - n11;
	  int n10 = Sa - n11;
	  int n01 = Sb - n11;
	  pvalue = ((int)round(pow(2, (int)p*d)) - (int)p*((int)round(pow(2, d)) - 1) - 1) * fisherExact(n11, n01, n10, n00);
	}

	if (pvalue > 1){
	  pvalue = 1;
	}
	if (pvalue < 0){
	  pvalue = 0;
	}

}

vector<vector<double>> imp(NumericMatrix& X)
{
  int r = X.nrow(), c = X.ncol();
  vector<vector<double>> v(r, vector<double>(c));
  for (int i = 0; i < r; i++){
    for (int j = 0; j < c; j++)
      v[i][j] = X(i, j);
  }
  return v;
}

//[[Rcpp::export]]
List symmCpp(NumericMatrix& X_R, int d, bool unif, bool test_uniformity, bool test_independence, List& independence_index)
{
  vector<vector<double>> X = imp(X_R);
  vector<vector<size_t>> idx;
  if(test_independence){
    for(auto& v: independence_index){
      idx.push_back(v);
    }
  }else{
    idx = {{}};
  }
  // vector<vector<size_t>> indp = {{}};

  BETfunction bet(X, d, unif, 1, test_uniformity, test_independence, idx);
  size_t p = X_R.ncol();
  // size_t length = (int)round(pow(2, (int)p*d));

  DataFrame symm = DataFrame::create(
    Named("SymmetryStatistics") = bet.getSymmStats()
  );

  // for (size_t i = p; i > 1; i--){
  //   vector<string>::const_iterator first = bet.getSymmInteraction().end() - length*(p-i+1);
  //   vector<string>::const_iterator last = bet.getSymmInteraction().end() - length*(p-i);
  //   vector<string> temp(first, last);
  //   symm.push_front(temp, "X" + to_string(i));
  // }
  //
  // vector<string> temp1(length);
  // for (size_t i = 0; i < length; i++){
  //   temp1[i] = bet.getSymmInteraction()[i];
  // }
  // symm.push_front(temp1, "X1");

  for (size_t i = p; i > 0; i--){
    symm.push_front(bet.getSymmInteraction()[i-1], "X" + to_string(i));
  }

  symm.push_front(bet.getBinary(), "Binary Index");

  List L = List::create(Named("SymmetryStatistics") = DataFrame(symm));

  return L;

}

//[[Rcpp::export]]
List BETCpp(NumericMatrix& X_R, int d, bool unif, bool asymptotic, bool test_uniformity, bool test_independence, List& independence_index)
{
  vector<vector<double>> X = imp(X_R);
  vector<vector<size_t>> idx;
  if(test_independence){
    for(auto& v: independence_index){
      idx.push_back(v);
    }
  }else{
    idx = {{}};
  }

  BETfunction bet(X, d, unif, asymptotic, test_uniformity, test_independence, idx);
  size_t n = X.size();
  double zstats = abs(bet.getStats())/sqrt(n);

  DataFrame df = DataFrame::create();

  size_t p = X[0].size();
  for (size_t i = p; i > 0; i--){
    // CharacterVector vi = {bet.getInteraction()[i-1]};
    df.push_front(bet.getInteraction()[i-1], "X" + to_string(i));
  }

  List L = List::create(Named("Interaction") = DataFrame(df), Named("Extreme.Asymmetry") = bet.getStats(), Named("p.value.bonf") = bet.getPvalue(), Named("z.statistic") = zstats);

  return L;

}

//[[Rcpp::export]]
List BeastCpp(NumericMatrix& X_R, int d, size_t m, size_t B, bool unif, double lambda, bool test_uniformity, bool test_independence, List& independence_index, String method, int numPerm)
{
  vector<vector<double>> X = imp(X_R);
  // independence index
  vector<vector<size_t>> idx;
  if(test_independence){
    for(auto& v: independence_index){
      idx.push_back(v);
    }
  }else{
    idx = {{}};
  }

  size_t n = X.size();
  size_t p = X[0].size();

  BETfunction bet(X, d, unif, 1, 1, 0, idx);

  bet.Beast(m, B, lambda, test_uniformity, test_independence, idx);

  double beastStat = bet.getBeastStat();

  //null distribution simulation
  vector<double> nullD(numPerm);

  vector<vector<double>> indp;
  NumericMatrix indp_R(n, p);

  // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // std::default_random_engine gen(seed);
  // std::normal_distribution<double> dis(0,1);

  vector<size_t> t(n);
  for(size_t i = 0; i < n; i++){
    t[i] = i;
  }

  double pv;

  // IntegerMatrix beastinter(p, d);
  // for(size_t i = 0; i < p; i++){
  //   for(size_t j = 0; j < d; j++)
  //     beastinter(i, j) = bet.getBeastInteraction()[i][j];
  // }

  // else if(method == "Y"){
  //   // user provide null distribution
  //   size_t l = null_simu.size();
  //   vector<double> nullsim(l);
  //   for(size_t i = 0; i < l; i++){
  //     nullsim[i] = null_simu[i];
  //   }
  //   pv = (double)count_if(nullsim.begin(), nullsim.end(), [&beastStat](double x) { return (x >= beastStat); }) / (double)l;
  //   List L = List::create(Named("Interaction") = bet.getBeastInteraction(), Named("BEAST.Statistic") = beastStat, Named("p.value") = pv);
  //   return L;
  // }

  if(method == "NA"){
    List L = List::create(Named("Interaction") = bet.getBeastInteraction(), Named("BEAST.Statistic") = beastStat);
    return L;
  }else{
    // 1-dim: only simulation
    if(p == 1){
      method = "s";
    }

    if(method == "p"){
      // permutation
      // for(size_t j = 0; j < p; j++){
      //   for(size_t i = 0; i < n; i++){
      //     indp[i][j] = runif(1);
      //   }
      // }
      for(size_t j = 0; j < p; j++){
        indp_R(_, j) = rnorm(n);
      }
      indp  = imp(indp_R);
      vector<vector<double>> newdata = indp;
      for(int sample = 0; sample < numPerm; sample++){
        // uniformity:
        if(test_uniformity){
          for(size_t col = 1; col < p; col++){

            // shuffle 0~n
            // random_device rd;
            // mt19937 gen(rd());
            // shuffle(t.begin(), t.end(), gen);
            random_shuffle(t.begin(), t.end(), randWrapper);

            for(size_t i = 0; i < n; i++){
              newdata[i][col] = indp[t[i]][col];
            }
          }
          BETfunction bet0(newdata, d, unif, 1, 1, 0, idx);
          bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
          nullD[sample] = bet0.getBeastStat();
        }else{
          // independence:
          for(size_t g = 1; g < idx.size(); g++){

            // shuffle 0~n
            // random_device rd;
            // mt19937 gen(rd());
            // shuffle(t.begin(), t.end(), gen);
            random_shuffle(t.begin(), t.end(), randWrapper);

            for(size_t v = 0; v < idx[g].size(); v++){
              for(size_t i = 0; i < n; i++){
                newdata[i][idx[g][v]-1] = indp[t[i]][idx[g][v]-1];
              }
            }
          }
        }
        // null distribution:
        BETfunction bet0(newdata, d, unif, 1, 1, 0, idx);
        bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
        nullD[sample] = bet0.getBeastStat();
      }
    }else if(method == "s"){
      // simulation
      for(int sample = 0; sample < numPerm; sample++){
        // for(size_t j = 0; j < p; j++){
        //   for(size_t i = 0; i < n; i++){
        //     indp[i][j] = dis(gen);
        //   }
        // }
        for(size_t j = 0; j < p; j++){
          indp_R(_, j) = rnorm(n);
        }
        indp = imp(indp_R);
        BETfunction bet0(indp, d, unif, 1, 1, 0, idx);
        bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
        nullD[sample] = bet0.getBeastStat();
      }
    }

    pv = (double)count_if(nullD.begin(), nullD.end(), [&beastStat](double x) { return (x >= beastStat); }) / (double)numPerm;

    List L = List::create(Named("Interaction") = bet.getBeastInteraction(), Named("BEAST.Statistic") = beastStat, Named("Null.Distribution") = nullD, Named("p.value") = pv);

    return L;
  }


}

//[[Rcpp::export]]
NumericVector nullCpp(size_t n, size_t p, int d, size_t m, size_t B, double lambda, bool test_uniformity, bool test_independence, List& independence_index, String method, int numPerm)
{
  // independence index
  vector<vector<size_t>> idx;
  if(test_independence){
    for(auto& v: independence_index){
      idx.push_back(v);
    }
  }else{
    idx = {{}};
  }

  //null distribution simulation
  NumericVector nullD (numPerm);

  vector<vector<double>> indp;
  NumericMatrix indp_R(n, p);

  vector<size_t> t(n);
  for(size_t i = 0; i < n; i++){
    t[i] = i;
  }

  // 1-dim: only simulation
  if(p == 1){
    method = "s";
  }

  if(method == "p"){
    // permutation
    for(size_t j = 0; j < p; j++){
      indp_R(_, j) = rnorm(n);
    }
    indp  = imp(indp_R);
    vector<vector<double>> newdata = indp;
    for(int sample = 0; sample < numPerm; sample++){
      // uniformity:
      if(test_uniformity){
        for(size_t col = 1; col < p; col++){

          // shuffle 0~n
          random_shuffle(t.begin(), t.end(), randWrapper);

          for(size_t i = 0; i < n; i++){
            newdata[i][col] = indp[t[i]][col];
          }
        }
        BETfunction bet0(newdata, d, 0, 1, 1, 0, idx);
        bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
        nullD[sample] = bet0.getBeastStat();
      }else{
        // independence:
        for(size_t g = 1; g < idx.size(); g++){

          // shuffle 0~n
          // random_device rd;
          // mt19937 gen(rd());
          // shuffle(t.begin(), t.end(), gen);
          random_shuffle(t.begin(), t.end(), randWrapper);

          for(size_t v = 0; v < idx[g].size(); v++){
            for(size_t i = 0; i < n; i++){
              newdata[i][idx[g][v]-1] = indp[t[i]][idx[g][v]-1];
            }
          }
        }
      }
      // null distribution:
      BETfunction bet0(newdata, d, 0, 1, 1, 0, idx);
      bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
      nullD[sample] = bet0.getBeastStat();
    }
  }else if(method == "s"){
    // simulation
    for(int sample = 0; sample < numPerm; sample++){
      // for(size_t j = 0; j < p; j++){
      //   for(size_t i = 0; i < n; i++){
      //     indp[i][j] = dis(gen);
      //   }
      // }
      for(size_t j = 0; j < p; j++){
        indp_R(_, j) = rnorm(n);
      }
      indp = imp(indp_R);
      BETfunction bet0(indp, d, 0, 1, 1, 0, idx);
      bet0.Beast(m, B, lambda, test_uniformity, test_independence, idx);
      nullD[sample] = bet0.getBeastStat();
    }
  }

  return nullD;

}
