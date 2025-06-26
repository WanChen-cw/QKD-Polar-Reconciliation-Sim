/*!
 * \file
 * \brief Class tools::Frozenbits_generator.
 */
#ifndef FROZENBITS_GENERATOR_HPP_
#define FROZENBITS_GENERATOR_HPP_

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include <bitset>

class Frozenbits_generator
{
private:
	protected:
		const int K; /*!< Number of information bits in the frame. */
		const int N; /*!< Codeword size (or frame size). */
        int N2;
		std::string path;
		std::vector<uint32_t> best_channels; /*!< The best channels in a codeword sorted by descending order. */
        std::vector<uint32_t> error_sort; /*!< 表示每个子信道错误，越大表示错误概率越高. */

public:
	Frozenbits_generator(const int K, const int N, const std::string& path = " ");	
	~Frozenbits_generator();

	void generate(std::vector<bool>& frozen_bits);//set frozenbit site
    void generate(std::vector<bool>& frozen_bits,int Nf);//set frozenbit site
    void generate(std::vector<bool>& frozen_bits, int & Nf_real,int & Nf_curidx, int blength , const std::vector<int> &seg_s, int n_seg);//set frozenbit site
    //void generate(std::vector<bool>& frozen_bits, int K1);//segment set frozenbit site
    //Pre-select pc location
    //Based on reliability estimation at the position of the minimum code word spacing
	void generatepc(int p,int number_pc,std::vector<bool>& pc_bits);
    void generatepcdescend(int p, int number_pc, std::vector<bool>& pc_bits);

    void encodegenerate(const std::vector<bool>& pc_bits,std::vector<std::vector<int>>&  pcindex);

		
protected:
	void get_best_channels();
    void _get_best_channels();
};

Frozenbits_generator
::Frozenbits_generator(const int K, const int N, const std::string& path)
	:K(K),N(N),path(path), error_sort(N)
{
    N2 = N - K;
	this->get_best_channels();
}

Frozenbits_generator
::~Frozenbits_generator()
{
}

void Frozenbits_generator
::generate(std::vector<bool>& frozen_bits)
{
	
	std::fill(frozen_bits.begin(), frozen_bits.end(), true);
	for (auto i = 0; i < K; i++)
		frozen_bits[best_channels[i]] = false;
}

void Frozenbits_generator
::generate(std::vector<bool>& frozen_bits,int Nf)
{
    int K1 = N - Nf;
    std::fill(frozen_bits.begin(), frozen_bits.end(), true);
    for (auto i = 0; i < K1; i++)
        frozen_bits[best_channels[i]] = false;
}

void Frozenbits_generator
::generate(std::vector<bool>& frozen_bits, int &Nf_real ,int &Nf_curidx, int blength ,const std::vector<int> &seg_s, int n_seg)
{
    int c = N / std::pow(2, n_seg-1);
    int idx = N;
    if (seg_s[0] == 1) {
        idx = 0;
    }
    for (int i = 1; i < n_seg; i++) {
        if (seg_s[i] == 1) {
            idx=c* std::pow(2, i-1);
            break;
        }
    }

    for (auto i = 0; i < blength; i++) {
        if (best_channels[N-Nf_curidx-i]>=idx) {
            frozen_bits[best_channels[N - Nf_curidx - i]] = true;
            Nf_real++;
        }
    }
    Nf_curidx += blength;
}

void Frozenbits_generator
::generatepc(int p, int number_pc, std::vector<bool>& pc_bits)
{
    int f = number_pc;

    std::vector<uint32_t> w(N, 1);
    w[1] = 2;
    for (int i = 1; i < static_cast<int>(log2(N)); ++i) {
        int start = pow(2, i);
        int end = pow(2, i + 1)-1;
        for (int j = start; j <= end; ++j) {
            w[j] = w[j-start] * 2;
        }
    }
    std::vector<uint32_t> buff1(K+f, 0);
    for (int i = 0; i < K + f; i++)
    {
        buff1[i] = w[best_channels[i]];
    }
    int wmin = *std::min_element(buff1.begin(), buff1.end());
    std::vector<uint32_t> site_wmin1;
    std::vector<uint32_t> site_wmin2;
    for (int i = 0; i < K + f; ++i) {
        if (w[best_channels[i]] == wmin) {
            site_wmin1.push_back(best_channels[i]);
        }
        else if (w[best_channels[i]] == 2 * wmin) {
            site_wmin2.push_back(best_channels[i]);
        }
    }
    int number_wmin = site_wmin1.size();
    int f1;
    int f2;
    if (f <= number_wmin) {
        f1 = f;
        f2 = 0;
    }
    else {
        f1 = number_wmin;
        f2 = f - number_wmin;
    }
    std::fill(pc_bits.begin(), pc_bits.end(), true);
    for (auto i = 0; i < K+f; i++)
        pc_bits[best_channels[i]] = false;
    for (int i = 0; i < f1; i++)
    {
        pc_bits[site_wmin1[i]]=true;
    }
    for (int i = 0; i < f2; i++)
    {
        pc_bits[site_wmin2[i]] = true;
    }
}

void Frozenbits_generator
::generatepcdescend(int p, int number_pc, std::vector<bool>& pc_bits)
{
    int f = number_pc;

    std::vector<uint32_t> w(N, 1);
    w[1] = 2;
    for (int i = 1; i < static_cast<int>(log2(N)); ++i) {
        int start = pow(2, i);
        int end = pow(2, i + 1) - 1;
        for (int j = start; j <= end; ++j) {
            w[j] = w[j - start] * 2;
        }
    }
    std::vector<uint32_t> buff1(K + f, 0);
    for (int i = 0; i < K + f; i++)
    {
        buff1[i] = w[best_channels[i]];
    }
    int wmin = *std::min_element(buff1.begin(), buff1.end());
    std::vector<uint32_t> site_wmin1;
    std::vector<uint32_t> site_wmin2;
    for (int i = 0; i < K + f; ++i) {
        if (w[best_channels[i]] == wmin) {
            site_wmin1.push_back(best_channels[i]);
        }
        else if (w[best_channels[i]] == 2 * wmin) {
            site_wmin2.push_back(best_channels[i]);
        }
    }
    int number_wmin = site_wmin1.size();
    /*int f1 = site_wmin1.size();
    int f2 = site_wmin2.size();*/
    int f1;
    int f2;
    if (f <= number_wmin) {
        f1 = f;
        f2 = 0;
    }
    else {
        f1 = number_wmin;
        f2 = f - number_wmin;
    }
    std::fill(pc_bits.begin(), pc_bits.end(), true);
    for (auto i = 0; i < K + f1+f2; i++)
        pc_bits[best_channels[i]] = false;


    //std::sort(site_wmin1.begin(), site_wmin1.end());
    //int x = site_wmin1.size() / f;
    //for (int i = 0; i < f1; i++)
    //{
    //    pc_bits[site_wmin1[i*x]] = true;
    //}
    for (int i = 0; i < f1; i++)
    {
        pc_bits[site_wmin1[site_wmin1.size()-1-i]] = true;
    }
    for (int i = 0; i < f2; i++)
    {
        pc_bits[site_wmin2[site_wmin2.size() - 1 - i]] = true;
    }
}

void Frozenbits_generator
::encodegenerate(const std::vector<bool>& pc_bits, std::vector<std::vector<int>>& pcindex)
{
    int idx1 = 0;
    for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
        if (pc_bits[encodepc_loop] == true) {
            pcindex[idx1].push_back(encodepc_loop);
            idx1++;
        }
    }
    for (int idx = 1; idx < N2; idx++) {
        if ((pcindex[idx][0] - pcindex[idx - 1][0] - 1) <= N2 - idx) {
            int insertidx = 0;
            for (int infoidx = pcindex[idx - 1][0] + 1; infoidx < pcindex[idx][0]; infoidx++) {
                pcindex[idx + insertidx].push_back(infoidx);
                insertidx++;
            }
        }
        else {
            std::vector<int> path_ordered(pcindex[idx][0] - pcindex[idx - 1][0] - 1);
            for (int i = 0; i < path_ordered.size(); ++i) {
                path_ordered[i] = pcindex[idx - 1][0] + i + 1;
            }
            std::sort(path_ordered.begin(), path_ordered.end(),
                [this](int a, int b) { return this->error_sort[a] > this->error_sort[b]; });
            std::sort(path_ordered.begin(), path_ordered.begin()+ N2 - idx);
            for (int i=0 ; i< N2 - idx; i++) {
                pcindex[idx+ i].push_back(path_ordered[i]);
            }
        }
    }
}




void Frozenbits_generator
::get_best_channels()
{
	std::ifstream inputFile(this->path);
	if (!inputFile.is_open()) {
		std::cout << "Failed to open channel parameter file: " << this->path << std::endl;
		return;
	}
	uint32_t value;
	while (inputFile >> value) {
		this->best_channels.push_back(value);
	}
	inputFile.close();
    for (int idx = 0; idx < best_channels.size(); idx++)
    {
        error_sort[best_channels[idx]] = idx;
    }
}

void Frozenbits_generator
::_get_best_channels()
{
    std::vector<double> weight(N,0);
    for (int i = 0; i < N; i++) {
        std::bitset<64> bitset2(i);
    
        for (int j = 0; j < 64; j++) {
            weight[i] += bitset2[j] * pow(1.1892,j);
        }
    }
    for (int i = 0; i <N; ++i) {
        this->best_channels.push_back(i) ;
    }
    std::sort(this->best_channels.begin(), this->best_channels.end(),
        [&weight](int a, int b) { return weight[a] > weight[b]; });
    for (int idx = 0; idx < best_channels.size(); idx++)
    {
        error_sort[best_channels[idx]] = idx;
    }
}
#endif /* FROZENBITS_GENERATOR_HPP_ */
