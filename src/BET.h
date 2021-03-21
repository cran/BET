/*
 * BET.h
 *
 *  Created on: Oct 20, 2020
 *      Author: wanzh
 */

#ifndef BET_H_
#define BET_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <random>
#include <Rcpp.h>

namespace N
{
	class BETfunction
	{
	   public:
		int d;
		int getStats();
		double getPvalue();
		std::vector<std::string> getInteraction();
		std::vector<std::vector<std::string>> getSymmInteraction();
		std::vector<int> getSymmStats();
		std::vector<std::string> getBinary();
		void Beast(std::size_t m, std::size_t B, double lambda, bool test_uniformity, bool test_independence, std::vector<std::vector<std::size_t>>& independence_index);
		double getBeastStat();
		std::vector<std::vector<int>> getBeastInteraction();
		BETfunction(std::vector<std::vector<double>>& X_R, int depth, bool unif, bool asymptotic, bool test_uniformity, bool test_independence, std::vector<std::vector<std::size_t>>& independence_index);

	   private:
		bool unifMargin;
    bool asympt;
    bool testUnif;
    bool testIndep;
    std::vector<std::vector<std::size_t>> indepIndex;
		std::size_t numThread = 1;
		std::size_t p;
		std::size_t n;
		std::vector<std::vector<double>> X;
		std::vector<std::string> binary_inter;
		// std::vector<std::string> symminter;
		std::vector<std::vector<std::string>> out_symminter;
		std::vector<int> out_symmstats;
		std::vector<long long> bid;
		std::vector<std::string> inter;
		std::vector<std::vector<int>> inter_mat;
		std::vector<std::vector<std::size_t>> allidx;
		std::vector<std::size_t> thread;
		std::map<std::vector<int>, int> countInteraction;
		std::vector<int> freqInter;
		std::vector<std::vector<int>> BeastInter;
		int Stats = 0;
		double BeastStat = 0;
		double pvalue = 0;
		// std::vector<std::vector<double>> imp(Rcpp::NumericMatrix& X);
		std::vector<std::string> interaction_str;
		std::vector<std::vector<int>> interactions();
		std::vector<std::vector<int>> interaction_mat();
		std::vector<std::string> interaction_index();
		std::string binaryToNum(std::string binary);
		std::vector<long long> BIDs();
		std::vector<std::vector<int>> ecdf_loc(std::vector<std::vector<double>>& X);
		std::map<std::vector<int>, int> groupC(std::vector<std::vector<int>>& c);
		int locate(int c, long long b);
		std::vector<std::vector<std::vector<int>>> CBIDs(std::map<std::vector<int>, int>& count);
		std::vector<std::vector<std::size_t>> allComb();
		bool isIndex(std::vector<std::size_t>& idx);
		std::vector<int> symmstats(std::vector<int>& countValues, std::vector<std::vector<std::vector<int>>>& matrix, size_t sampleSize);
		std::vector<std::size_t> unreplaceShuffle(std::size_t size, std::size_t max_size);
		std::vector<double> subsample(size_t m, size_t B);
		std::vector<double> softthreshold(std::vector<double>& S, double lambda);
		double logHypergeometricProb(double* logFacs, int a, int b, int c, int d);
		double fisherExact(int a, int b, int c, int d);
		double binomial(int n, int k, double p);
		double pbinom(int n, int k, double p);
	};
}



#endif /* BET_H_ */
