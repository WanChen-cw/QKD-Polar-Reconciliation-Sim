/*!
 * \file
 * \brief Class module::Decoder_polar_fast.
 */
#ifndef DECODER_POLAR_FAST_HPP_
#define DECODER_POLAR_FAST_HPP_
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Mycrc.hpp"


template <typename B = int, typename R = float>
class Decoder_polar_fast
{
private:
    int N;
    std::vector<int> lambda_offset;
    std::vector<int> llr_layer_vec;
    std::vector<int> bit_layer_vec;

    std::vector<std::vector<int>> node_type_structure;
    std::vector<std::vector<int>> code_structure;
    int cnt_structure;
    std::vector<int> psi_vec;
public:
    Decoder_polar_fast(const int& N);
    ~Decoder_polar_fast() {}

    bool fastscl(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        std::vector<B>& dec_bit,
        uint32_t crc, int L);
   
    void get_node_structure(const std::vector<bool>& frozen_bits);

    void get_psi_for_advanced_sc_decoder();
protected:
    void node_identifier(const std::vector<bool>& f, const std::vector<int>& z);
    void get_lambdaoffset();
    void get_bit_layer();
    void get_llr_layer();
    std::vector<int> generateIndices(int size);
    void encode(std::vector<B>& bits,int M){
        for (auto k = (M >> 1); k > 0; k >>= 1)
            for (auto j = 0; j < M; j += 2 * k)
                for (auto i = 0; i < k; i++)
                    bits[j + i] = bits[j + i] ^ bits[k + j + i];
    }
    R f(R L1, R L2);

    void Rate1(int M, int L,
        const std::vector<std::vector<R>>& LLR_array,
        std::vector<int>& activepath,
        std::vector<R>& PM,
        std::vector<std::vector<int>>& lazy_copy,
        std::vector<std::vector<B>>& candidate_codeword);

    void SPC(int M, int L,
        const std::vector<std::vector<R>>& LLR_array,
        std::vector<int>& activepath,
        std::vector<R>& PM,
        std::vector<std::vector<int>>& lazy_copy,
        std::vector<std::vector<B>>& candidate_codeword);

    void REP(int M, int L,
        const std::vector<std::vector<R>>& LLR_array,
        std::vector<int>& activepath,
        std::vector<R>& PM,
        std::vector<std::vector<int>>& lazy_copy,
        std::vector<std::vector<B>>& candidate_codeword);
};

template<typename B, typename R>
Decoder_polar_fast<B, R>
::Decoder_polar_fast(const int& N)
    :N(N), bit_layer_vec(N, 0), llr_layer_vec(N, 0)
{
    this->get_lambdaoffset();
    this->get_bit_layer();
    this->get_llr_layer();
}

template<typename B, typename R>
bool Decoder_polar_fast<B, R>
::fastscl(const std::vector<R>& llr, const std::vector<B>& frozen_bits_value,std::vector<B>& dec_bit, uint32_t crc, int L)
{
    int node_size = node_type_structure.size();
    bool right = false;
    int m = log2(N);
    std::vector<std::vector<R>>         P(L, std::vector<R>(N-1, 0.0));
    std::vector<std::vector<B>>         C(2*L, std::vector<B>(2*N-1, 0));
    std::vector<R>                      PM(L, 0.0);
    std::vector<int>                    activepath(L, 0);
    activepath[0] = 1;
    std::vector<std::vector<int>>       lazy_copy(L, std::vector<int>(m, 0));
    for (int i = 0; i < m; ++i) {
        lazy_copy[0][i] = 0;
    }

    //BUFF
    int index_1;
    int index_2;
    int llr_layer;
    int current_index;
    int M ;
    int reduced_layer;
    int psi_mod_2 ;
    std::vector<std::vector<R>> sub_P(L);
    std::vector<std::vector<B>> candidate_codeword(L);
    std::vector<B> ux;
    for (volatile int i_node = 0; i_node < node_type_structure.size(); ++i_node) {
        //llr update
        M = node_type_structure[i_node][1];
        reduced_layer = log2(M);
        psi_mod_2 = psi_vec[i_node] % 2;
        current_index = node_type_structure[i_node][0];
        llr_layer = llr_layer_vec[current_index];
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if (current_index == 0) {
                index_1 = lambda_offset[m-1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[l_index][beta + index_1 - 1] = f(llr[beta], llr[beta + index_1]);
                }
                for (int i_layer = m - 2; i_layer >= reduced_layer; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < index_2; ++beta) {
                        P[l_index][beta- 1] = f(P[l_index][beta+index_1-1], P[l_index][beta+index_2-1]);
                    }
                }
            }
            else if (current_index == N / 2) {
                index_1 = lambda_offset[m-1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[l_index][beta + index_1-1] = (1-2*C[2*l_index][beta+index_1-1])*llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= reduced_layer; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < index_2; ++beta) {
                        P[l_index][beta- 1] = f(P[l_index][beta + index_1 - 1], P[l_index][beta + + index_2 - 1]);
                    }
                }
            }
            else {
                index_1 = lambda_offset[llr_layer];
                index_2 = lambda_offset[llr_layer + 1];
                for (int beta = index_1; beta < index_2; ++beta) {
                    P[l_index][beta-1] = (1 - 2 * C[2 * l_index][beta-1]) *P[lazy_copy[l_index][llr_layer + 1]][beta + index_1 - 1]
                        + P[lazy_copy[l_index][llr_layer + 1]][beta + index_2- 1];
                }
                for (int i_layer = llr_layer - 1; i_layer >= reduced_layer; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < index_2; ++beta) {
                        P[l_index][beta-1] = f(P[l_index][beta + index_1- 1], P[l_index][beta+index_2 - 1]);
                    }
                }
            }
        }

        //node process
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            sub_P[l_index].resize(M);
            std::copy(P[l_index].begin() + M - 1, P[l_index].begin() + M - 1 + M, sub_P[l_index].begin());
        }
        for (int l_index = 0; l_index < L; ++l_index) {
            candidate_codeword[l_index].resize(M);
        }
        std::copy(frozen_bits_value.begin() + current_index, frozen_bits_value.begin() + current_index + M, candidate_codeword[0].begin());
        switch (node_type_structure[i_node][2]) {
        case -1: // RATE 0
            ux.resize(M);
            std::copy(frozen_bits_value.begin() + current_index, frozen_bits_value.begin() + current_index + M,ux.begin());
            encode(ux, M);
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                std::copy(ux.begin(), ux.end(), C[psi_mod_2 + 2 * l_index].begin() + M -1);
                for (int bete = 0; bete < M; bete++) {
                    PM[l_index] += log(1+exp(-(1-2*C[psi_mod_2+2*l_index][M+bete-1]) * P[l_index][M + bete - 1]));
                }
            }
            break;
        case 1: // RATE 1
            Rate1(M,L, sub_P, activepath, PM, lazy_copy, candidate_codeword);
            for(int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                std::copy(candidate_codeword[l_index].begin(), candidate_codeword[l_index].end(), C[psi_mod_2 + 2 * l_index].begin() + M - 1);
            }
            break;
        case 2: // REP
            REP(M, L, sub_P, activepath, PM, lazy_copy, candidate_codeword);
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                std::copy(candidate_codeword[l_index].begin(), candidate_codeword[l_index].end(), C[psi_mod_2 + 2 * l_index].begin() + M - 1);
            }
            break;
        case 3: // SPC------------
            SPC(M, L, sub_P, activepath, PM, lazy_copy, candidate_codeword);
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                std::copy(candidate_codeword[l_index].begin(), candidate_codeword[l_index].end(), C[psi_mod_2 + 2 * l_index].begin() + M - 1);
            }
            break;
        }

        //partial-sum return------------------------
        int bit_layer = bit_layer_vec[current_index + M - 1];
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if (psi_mod_2 == 1) {
                for (int i_layer = reduced_layer; i_layer < bit_layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[2 * l_index + 1][beta + index_1 - 1] = (C[2 * lazy_copy[l_index][i_layer]][beta - 1] + C[2 * l_index + 1][beta - 1]) % 2;
                        C[2 * l_index + 1][beta + index_2 - 1] = C[2 * l_index + 1][beta - 1];
                    }
                }
                index_1 = lambda_offset[bit_layer];
                index_2 = lambda_offset[bit_layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[2 * l_index][beta + index_1 - 1] = (C[2 * lazy_copy[l_index][bit_layer]][beta - 1] + C[2 * l_index + 1][beta - 1]) % 2; // 左列延迟复制
                    C[2 * l_index][beta + index_2 - 1] = C[2 * l_index + 1][beta - 1];
                }
            }
        }

        // Lazy copy------------------
        if (i_node < node_size-1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[current_index + M] ; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[l_index][i_layer] = l_index;
                }
            }
        }
    }

    //path selection
    std::vector<int> path_ordered(L);
    for (int i = 0; i < L; ++i) {
        path_ordered[i] = i; // Initialize with path indices
    }
    std::sort(path_ordered.begin(), path_ordered.end(),
        [&PM](int a, int b) { return PM[a] < PM[b]; });
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::copy(C[2*path_num].begin()+N-1, C[2*path_num].end(), dec_bit.begin());
        encode(dec_bit,N);
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}


template<typename B, typename R>
void Decoder_polar_fast<B, R>
::get_psi_for_advanced_sc_decoder()
{
    psi_vec.resize(node_type_structure.size(), 0);
    for (size_t i = 0; i < psi_vec.size(); ++i) {
        int psi = node_type_structure[i][0];
        int M = node_type_structure[i][1];
        int reduced_layer = static_cast<int>(log2(M));

        for (int j = 0; j < reduced_layer; ++j) {
            psi = psi / 2;
        }
        psi_vec[i] = psi;
    }
}

template <typename B, typename R>
void  Decoder_polar_fast<B, R>
::get_node_structure(const std::vector<bool>& frozen_bits) {
    node_type_structure.clear();  // 清空节点结构
    code_structure.resize(frozen_bits.size(), std::vector<int>(3, 0));
    cnt_structure = 1;
    //node_identifier(frozen_bits, std::vector<int>(frozen_bits.begin(), frozen_bits.end()));
    node_identifier(frozen_bits, generateIndices(frozen_bits.size()));


    // Reshape the code_structure
    for (const auto& row : code_structure) {
        if (row[2] != 0) {
            node_type_structure.push_back(row);
        }
    }
}

template <typename B, typename R>
void Decoder_polar_fast<B, R>
::get_lambdaoffset()
{
    int numBits = static_cast<int>(std::log2(this->N)) + 1;
    int powerOfTwo = 1;
    for (int i = 0; i < numBits; ++i) {
        this->lambda_offset.push_back(powerOfTwo);
        powerOfTwo *= 2;
    }
}

template <typename B, typename R>
void Decoder_polar_fast<B, R>
::get_bit_layer()
{
    for (int phi = 0; phi < this->N; ++phi) {
        int psi = phi / 2;
        int layer = 0;
        while (psi % 2 == 1) {
            psi = psi / 2;
            layer++;
        }
        this->bit_layer_vec[phi] = layer;
    }
}

template <typename B, typename R>
void Decoder_polar_fast<B, R>
::get_llr_layer()
{
    for (int phi = 1; phi < this->N; ++phi) {
        int psi = phi;
        int layer = 0;
        while (psi % 2 == 0) {
            psi = psi / 2;
            layer++;
        }
        this->llr_layer_vec[phi] = layer;
    }
}

template<typename B, typename R>
R Decoder_polar_fast<B, R>
::f(R L1, R L2)
{
    return std::copysign(1.0, L1) * std::copysign(1.0, L2) * std::min(std::fabs(L1), std::fabs(L2));
}

template<typename B, typename R>
std::vector<int> Decoder_polar_fast<B, R>
::generateIndices(int size) {
    std::vector<int> indices(size);
    for (int i = 0; i < size; ++i) {
        indices[i] = i;
    }
    return indices;
}

template <typename B, typename R>
void Decoder_polar_fast<B, R>
::node_identifier(const std::vector<bool>& f, const std::vector<int>& z) {
    int N = f.size();

    if (std::all_of(f.begin(), f.end() - 1, [](int x) { return x == 1; }) && (f.back() == 0)) {
        code_structure[cnt_structure][0] = z[0];
        code_structure[cnt_structure][1] = N;
        code_structure[cnt_structure][2] = 2;
        cnt_structure++;
    }
    else {
        if ((f[0] == 1) && std::all_of(f.begin() + 1, f.end(), [](int x) { return x == 0; })) {
            code_structure[cnt_structure][0] = z[0];
            code_structure[cnt_structure][1] = N;
            code_structure[cnt_structure][2] = 3;
            cnt_structure++;
        }
        else {
            if (std::all_of(f.begin(), f.end(), [](int x) { return x == 0; })) {
                code_structure[cnt_structure][0] = z[0];
                code_structure[cnt_structure][1] = N;
                code_structure[cnt_structure][2] = 1;
                cnt_structure++;
            }
            else {
                if (std::all_of(f.begin(), f.end(), [](int x) { return x == 1; })) {
                    code_structure[cnt_structure][0] = z[0];
                    code_structure[cnt_structure][1] = N;
                    code_structure[cnt_structure][2] = -1;
                    cnt_structure++;
                }
                else {
                    node_identifier(std::vector<bool>(f.begin(), f.begin() + N / 2), std::vector<int>(z.begin(), z.begin() + N / 2));
                    node_identifier(std::vector<bool>(f.begin() + N / 2, f.end()), std::vector<int>(z.begin() + N / 2, z.end()));
                }
            }
        }
    }
}

template<typename B, typename R>
void Decoder_polar_fast<B, R>
::Rate1(int M, int L, const std::vector<std::vector<R>>& LLR_array, std::vector<int>& activepath, std::vector<R>& PM, std::vector<std::vector<int>>& lazy_copy, std::vector<std::vector<B>>& candidate_codeword)
{
    std::vector<std::vector<int>> llr_ordered(L, std::vector<int>(M, 0));
    std::vector<int> sub_path_ordered(M);
    std::vector<int> sub_lazy_copy(L, 0);
    for (int i = 0; i < L; ++i) {
        sub_lazy_copy[i] = i;
    }
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;

        for (int i = 0; i < M; ++i) {
            sub_path_ordered[i] = i; 
        }
        std::sort(sub_path_ordered.begin(), sub_path_ordered.end(),[&](int a, int b) { return std::abs(LLR_array[l_index][a]) < std::abs(LLR_array[l_index][b]); });
        llr_ordered[l_index] = sub_path_ordered;
        for (int i = 0; i < M; ++i) {
            candidate_codeword[l_index][i] = LLR_array[l_index][i] < 0 ? 1 : 0;
        }
    }

    for (int i_split = 0; i_split < std::min(M, L - 1); ++i_split) {
        std::vector<std::vector<R>> PM_rate_1(2, std::vector<R>(L, std::numeric_limits<R>::max()));
        std::vector<std::vector<int>> compare_rate_1(2, std::vector<int>(L, 0));

        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            PM_rate_1[0][l_index] = PM[l_index];
            PM_rate_1[1][l_index] = PM[l_index] + std::abs(LLR_array[sub_lazy_copy[l_index]][llr_ordered[sub_lazy_copy[l_index]][i_split]]);
        }

        int num_surviving_path = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
        std::vector<int> PM_inx(2 * L);
        for (int i = 0; i < PM_inx.size(); ++i) {
            PM_inx[i] = i;
        }
        std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
            return PM_rate_1[a / L][a % L] < PM_rate_1[b / L][b % L];
            });
        for (int i = 0; i < num_surviving_path; ++i) {
            int row = PM_inx[i] / L;
            int col = PM_inx[i] % L;
            compare_rate_1[row][col] = 1;
        }
        int kill_cnt = 0;
        std::vector<int> kill_index;
        for (int i = 0; i < L; ++i) {
            if (compare_rate_1[0][i] == 0 && compare_rate_1[1][i] == 0) {
                activepath[i] = 0;
                kill_cnt++;
                kill_index.push_back(i);
            }
        }

        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((compare_rate_1[0][l_index]  + compare_rate_1[1][l_index])==2)
            {
                int new_index = kill_index[kill_cnt - 1];
                --kill_cnt;
                activepath[new_index] = 1;
                sub_lazy_copy[new_index] = sub_lazy_copy[l_index];
                candidate_codeword[new_index] = candidate_codeword[l_index];
                candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][i_split]] = candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][i_split]] == 1 ? 0 : 1;
                PM[new_index] = PM_rate_1[1][l_index];
            }
        }
    }

    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;
        lazy_copy[l_index] = lazy_copy[sub_lazy_copy[l_index]];
    }
}

template<typename B, typename R>
void Decoder_polar_fast<B, R>
::SPC(int M, int L, const std::vector<std::vector<R>>& LLR_array, std::vector<int>& activepath, std::vector<R>& PM, std::vector<std::vector<int>>& lazy_copy, std::vector<std::vector<B>>& candidate_codeword)
{
    B ux = candidate_codeword[0][0];
    std::vector<std::vector<int>> llr_ordered(L, std::vector<int>(M, 0));
    std::vector<int> sub_path_ordered(M,0);

    std::vector<int> parity_check_track(L, 0);
    std::vector<int> sub_lazy_copy(L, 0);
    for (int i = 0; i < L; ++i) {
        sub_lazy_copy[i] = i;
    }

    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;
        std::vector<R> abs_LLR(LLR_array[l_index].size());
        std::transform(LLR_array[l_index].begin(), LLR_array[l_index].end(), abs_LLR.begin(),[](R v) {return std::abs(v);});
        for (int i = 0; i < M; ++i) {
            sub_path_ordered[i] = i;
        }
        std::sort(sub_path_ordered.begin(), sub_path_ordered.end(), [&](int a, int b) { return abs_LLR[a] < abs_LLR[b]; });
        llr_ordered[l_index] = sub_path_ordered;
        for (int i = 0; i < M; ++i) {
            candidate_codeword[l_index][i] = LLR_array[l_index][i] < 0 ? 1 : 0;
        }
        if ((std::accumulate(candidate_codeword[l_index].begin(), candidate_codeword[l_index].end(), 0)  % 2)  == (ux != 1)) 
        {
            candidate_codeword[l_index][llr_ordered[l_index][0]] = candidate_codeword[l_index][llr_ordered[l_index][0]]==1?0:1;
            parity_check_track[l_index] = 1;
            PM[l_index] += abs_LLR[llr_ordered[l_index][0]];
        }
    }
    std::vector<R> LLR(LLR_array[0].size());
    for (int t = 1; t < std::min(M, L); ++t) {
        std::vector<std::vector<R>> PM_rate_1(2, std::vector<R>(L, std::numeric_limits<R>::max()));
        std::vector<std::vector<int>> compare_rate_1(2, std::vector<int>(L, 0));
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            LLR = LLR_array[sub_lazy_copy[l_index]];
            PM_rate_1[0][l_index] = PM[l_index];
            PM_rate_1[1][l_index] = PM[l_index] + std::abs(LLR[llr_ordered[sub_lazy_copy[l_index]][t]]) + (1 - 2 * parity_check_track[l_index]) * std::abs(LLR[llr_ordered[sub_lazy_copy[l_index]][0]]);
        }
        int num_surviving_path = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
        std::vector<int> PM_inx(2 * L);
        for (int i = 0; i < PM_inx.size(); ++i) {
            PM_inx[i] = i;
        }
        std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
            return PM_rate_1[a / L][a % L] < PM_rate_1[b / L][b % L];
            });
        for (int i = 0; i < num_surviving_path; ++i) {
            int row = PM_inx[i] / L;
            int col = PM_inx[i] % L;
            compare_rate_1[row][col] = 1;
        }
        int kill_cnt = 0;
        std::vector<int> kill_index;
        for (int i = 0; i < L; ++i) {
            if (compare_rate_1[0][i] == 0 && compare_rate_1[1][i] == 0) {
                activepath[i] = 0;
                kill_cnt++;
                kill_index.push_back(i);
            }
        }
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((compare_rate_1[0][l_index] + compare_rate_1[1][l_index]) == 2)
            {
                int new_index = kill_index[kill_cnt - 1];
                --kill_cnt;
                activepath[new_index] = 1;
                sub_lazy_copy[new_index] = sub_lazy_copy[l_index];
                candidate_codeword[new_index] = candidate_codeword[l_index];
                candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][t]] = candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][t]] == 1 ? 0 : 1;
                candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][0]] = candidate_codeword[new_index][llr_ordered[sub_lazy_copy[l_index]][0]] == 1 ? 0 : 1;
                parity_check_track[new_index] = parity_check_track[new_index] == 1 ? 0 : 1;
                PM[new_index] = PM_rate_1[1][l_index];
            }
        }
    }
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;
        lazy_copy[l_index] = lazy_copy[sub_lazy_copy[l_index]];
    }
}

template<typename B, typename R>
void Decoder_polar_fast<B, R>
::REP(int M, int L, const std::vector<std::vector<R>>& LLR_array, std::vector<int>& activepath, std::vector<R>& PM, std::vector<std::vector<int>>& lazy_copy, std::vector<std::vector<B>>& candidate_codeword)
{
    int num_surviving_path = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
    std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
    std::vector<B> candidate_codeword1 = candidate_codeword[0];
    candidate_codeword1[M - 1] = 0;
    encode(candidate_codeword1, M);
    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0)
            continue;
        R Delta0 = 0;
        R Delta1 = 0;
        for (int i_llr = 0; i_llr < M; ++i_llr) {
            if (candidate_codeword1[i_llr] == 0) {
                if (LLR_array[l_index][i_llr] < 0) {
                    Delta0 -= LLR_array[l_index][i_llr];
                }
                else {
                    Delta1 += LLR_array[l_index][i_llr];
                }
            }
            else {
                if (LLR_array[l_index][i_llr] > 0) {
                    Delta0 += LLR_array[l_index][i_llr];
                }
                else {
                    Delta1 -= LLR_array[l_index][i_llr];
                }
            }
        }
        PM_pair[0][l_index] = PM[l_index] + Delta0;
        PM_pair[1][l_index] = PM[l_index] + Delta1;
    }
    std::vector<int> PM_inx(2 * L);
    for (int i = 0; i < PM_inx.size(); ++i) {
        PM_inx[i] = i;
    }
    std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
        return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
        });
    std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
    for (int i = 0; i < num_surviving_path; ++i) {
        int row = PM_inx[i] / L;
        int col = PM_inx[i] % L;
        compare[row][col] = 1;
    }
    std::vector<int> kill_index;
    int kill_cnt = 0;
    for (int i = 0; i < L; ++i) {
        if (compare[0][i] == 0 && compare[1][i] == 0) {
            activepath[i] = 0;
            kill_cnt++;
            kill_index.push_back(i);
        }
    }

    for (int l_index = 0; l_index < L; ++l_index) {
        if (activepath[l_index] == 0) {
            continue;
        }
        int path_state = compare[0][l_index] * 2 + compare[1][l_index];
        switch (path_state) {
        case 1:
            std::transform(candidate_codeword1.begin(), candidate_codeword1.end(),
                candidate_codeword[l_index].begin(),
                [](B value) { return 1 - value; });
            PM[l_index] = PM_pair[1][l_index];
            break;
        case 2:
            candidate_codeword[l_index] = candidate_codeword1;
            PM[l_index] = PM_pair[0][l_index];
            break;
        case 3: {
            int new_index = kill_index[kill_cnt - 1];
            --kill_cnt;
            activepath[new_index] = 1;
            lazy_copy[new_index] = lazy_copy[l_index];
            candidate_codeword[l_index] = candidate_codeword1;
            std::transform(candidate_codeword1.begin(), candidate_codeword1.end(),
                candidate_codeword[new_index].begin(),
                [](B value) { return 1 - value; });
            PM[l_index] = PM_pair[0][l_index];
            PM[new_index] = PM_pair[1][l_index];
            break;
        }
        }
    }
}
#endif // DECODER_POLAR_FAST_HPP_


