/*!
 * \file
 * \brief Class module::Decoder_polar.
 */
#ifndef DECODER_POLAR2_HPP_
#define DECODER_POLAR2_HPP_

#include <map>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Mycrc.hpp"
//#include <Eigen/Core>



template <typename B = int>
std::vector<B> CSR(std::vector<B> reg) {
    B feedback = reg[reg.size()-1] ; // Get the LSB (bit 0) as feedback
    for (size_t i = reg.size() - 1; i > 0; --i) {
        reg[i] = reg[i - 1]; // Right-shift all bits
    }
    reg[0] = feedback; // Set the shifted-in bit
    return reg;
}

template <typename B = int, typename R = float>
class Decoder_polar
{
private:
    int N;
    int K;
    bool Q;
    std::vector<bool>  frozen_bits; // true means frozen, false means set to 0/1
    std::vector<int> lambda_offset;
    std::vector<int> llr_layer_vec;
    std::vector<int> bit_layer_vec;
public:
    Decoder_polar(const int& K, const int& N, const std::vector<bool>& frozen_bits, const bool Q = 0);
    ~Decoder_polar();

    const std::vector<bool>& get_frozen_bits() const;
    void set_frozen_bits(const std::vector<bool>& frozen_bits);

    //优化冻结位置下的译码
    bool decodepc2(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit, 
        std::vector<std::vector<int>>& pcindex,
        int p,uint32_t crc, int L, const size_t frame_id = 0);

    bool decodepc(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p,uint32_t crc, int L, const size_t frame_id = 0);


    bool decodepc8(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<B>& frozen_bits_value2,
        std::vector<int>& column1,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, const size_t frame_id = 0);

    //超出容量位置下一个循环检测
    bool decodepc7(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        std::vector<double> &z,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, const size_t frame_id = 0);

    //超出信道容量不保护
    bool decodepc6(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<B>& u,
        const std::vector<bool>& pcsite,
        const std::vector<bool>& pcsite_u,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, const size_t frame_id = 0);

    //编码优化挑空
    bool decodepc5(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, const size_t frame_id = 0);

    //不同长度校验位
    bool decodepc4(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<B>& frozen_bits_value2,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p,int p2, uint32_t crc, int L);

    //冻结位置先pm比较后删除路径
    bool decodepc3(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, const size_t frame_id = 0);

    bool decodepcsegment(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p,
        uint32_t crc,
        const std::vector<uint32_t> &crc_v,
        std::vector<int>& seg_s,
        int L, const size_t frame_id = 0);

    bool decode(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        std::vector<B>& dec_bit,
        uint32_t crc, int L, const size_t frame_id = 0);
    
    //
    bool decodescl(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& fsite,
        std::vector<B>& dec_bit,
        uint32_t crc, int L, const size_t frame_id = 0);

    bool decodescltest(const std::vector<R>& llr,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& fsite,
        std::vector<B>& dec_bit,
        uint32_t crc, int L,
        std::vector<std::vector<int>>& error_site);

    bool decodepctest(const std::vector<R>& llr,
        const std::vector<B>& pc,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L, 
        std::vector<std::vector<int>>& error_site);

    bool decodepctest2(const std::vector<R>& llr,
        const std::vector<B>& pc,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L,
        std::vector<std::vector<int>>& error_site);

    bool decodepctest3(const std::vector<R>& llr,
        const std::vector<B>& pc,
        const std::vector<B>& frozen_bits_value,
        const std::vector<bool>& pcsite,
        std::vector<B>& dec_bit,
        int p, uint32_t crc, int L,
        std::vector<std::vector<int>>& error_site);

protected:
    void get_lambdaoffset();
    void get_bit_layer();
    void get_llr_layer();

    R f(R L1, R L2);
};

template <typename B, typename R>
Decoder_polar<B, R>
::Decoder_polar(const int& K, const int& N, const std::vector<bool>& frozen_bits, const bool Q)
    : K(K), N(N), Q(Q), frozen_bits(N, 0), bit_layer_vec(N, 0), llr_layer_vec(N, 0)
{
    this->get_lambdaoffset();
    this->get_bit_layer();
    this->get_llr_layer();
    this->set_frozen_bits(frozen_bits);
}

template <typename B, typename R>
Decoder_polar<B, R>
::~Decoder_polar()
{
}

template <typename B, typename R>
const std::vector<bool>& Decoder_polar<B, R>
::get_frozen_bits() const
{
    return this->frozen_bits;
}

template <typename B, typename R>
void Decoder_polar<B, R>
::set_frozen_bits(const std::vector<bool>& frozen_bits)
{
    std::copy(frozen_bits.begin(), frozen_bits.end(), this->frozen_bits.begin());
}


template<typename B, typename R>
bool Decoder_polar<B, R>
::decode(const std::vector<R>& llr,
    const std::vector<B>& frozen_bits_value,
    std::vector<B>& dec_bit,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = frozen_bits_value;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L, 0.0));
    std::vector<std::vector<B>>       C(N - 1, std::vector<B>(2 * L, 0));
    std::vector<std::vector<B>>       u(L, std::vector<B>(K, 0));
    std::vector<R>                      PM(L, 0.0);
    std::vector<int>                    activepath(L, 0);
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][lazy_copy[layer + 1][l_index]] + P[beta + index_1 + index_2 - 1][lazy_copy[layer + 1][l_index]];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        //if now we decode an unfrozen bit
        if (frozen_bits[phi] == 0) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }

            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<size_t> PM_inx(2 * L);
            for (size_t i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](size_t a, size_t b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (size_t i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                        //u[index] = u[l_index];
                    
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    break;
                }
                }
            }
            cnt_u = cnt_u + 1;
        }
        //frozen bit operation
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                if (P[0][l_index] < 0 && frozen_bits_value[phi] == 0) {
                    PM[l_index] = PM[l_index] - P[0][l_index];
                }
                if (P[0][l_index] > 0 && frozen_bits_value[phi] == 1) {
                    PM[l_index] = PM[l_index] + P[0][l_index];
                }
                if (phi_mod_2 == 0) {
                    C[0][2 * l_index] = frozen_bits_value[phi];
                }
                else {
                    C[0][2 * l_index + 1] = frozen_bits_value[phi];
                }
            }
        }

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }


        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(K);
        for (int i = 0; i < K; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            if (frozen_bits[i] == 0) {
                dec_bit[i] = polar_info_esti[j];
                j++;
            }
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
        /*else {
            j = 0;
            for (int i = 0; i < N; ++i) {
                if (frozen_bits[i] == 0) {
                    polar_info_esti_all[i] = polar_info_esti[j];
                    j++;
                }
            }
        }*/
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));

    /*int regsite_cur = 0;
    int regsite_last = -1;*/

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                //if (regsite_cur != (regsite_last + 1) % p) {
                //    for (int x = 0; x < regsite_last + p + 1 - regsite_cur; x++)
                //        reg[l_index] = CSR(reg[l_index]);
                //}

                if (pc[phi] == reg[l_index][0]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    reg[l_index] = CSR(reg[l_index]);
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                }

            }
            //regsite_last = regsite_cur;
        }
        cnt_u = cnt_u + 1;
        //regsite_cur = (regsite_cur + 1) % p;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc8(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& PCE,
    std::vector<int>& column1,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));

    /*int regsite_cur = 0;
    int regsite_last = -1;*/

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                if (phi == column1[9]) {
                    int xxxx = 0;
                    for (int encodepc_loop = 0; encodepc_loop < 9; ++encodepc_loop) {
                        xxxx = xxxx ^ u[l_index][column1[encodepc_loop]]; // XOR operation
                    }
                    if (PCE[0] == xxxx) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                }
                else if (phi == column1[19]) {
                    int xxxx = 0;
                    for (int encodepc_loop = 0; encodepc_loop < 9; ++encodepc_loop) {
                        xxxx = xxxx ^ u[l_index][column1[encodepc_loop+10]]; // XOR operation
                    }
                    if (PCE[1] == xxxx) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                }
                else {
                if (pc[phi] == reg[l_index][0]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    reg[l_index] = CSR(reg[l_index]);
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                }
                }
            }
            //regsite_last = regsite_cur;
        }
        cnt_u = cnt_u + 1;
        //regsite_cur = (regsite_cur + 1) % p;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc7(const std::vector<R>& llr,
    const std::vector<B>& pc,
    std::vector<double>& z,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));
    std::vector<double>                 z_reg(p, 0);
    /*int regsite_cur = 0;
    int regsite_last = -1;*/

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }

            z_reg[0] = z_reg[0] + z[phi];
        }
        else {
            if (z_reg[0] + z[phi] <= 1) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    if (activepath[l_index] == 0)
                        continue;
                    if (pc[phi] == reg[l_index][0]) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    }
                    reg[l_index] = CSR(reg[l_index]);
                }
                z_reg[0] = 0;
            }
            else {
                for (int l_index = 0; l_index < L; ++l_index) {
                    if (activepath[l_index] == 0)
                        continue;
                    if (pc[phi] == 0) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                    }
                    reg[l_index] = CSR(reg[l_index]);
                }
            }
        }
        z_reg = CSR(z_reg);
        cnt_u = cnt_u + 1;
        //regsite_cur = (regsite_cur + 1) % p;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc6(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& u_pc,
    const std::vector<bool>& pcsite,
    const std::vector<bool>& pcsite_u,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));

    /*int regsite_cur = 0;
    int regsite_last = -1;*/

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                if (pcsite_u[phi] == false) {
                    if (pc[phi] == reg[l_index][0]) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                }
                else {
                    if (u_pc[phi] == 0) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg[l_index] = CSR(reg[l_index]);
                    }
                }

            }
            //regsite_last = regsite_cur;
        }
        cnt_u = cnt_u + 1;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc5(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<B>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<B>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<B>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<B>                      activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));
    std::vector<B>                      last_statu(p, 0);
    activepath[0] = 1;

    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }
    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;
    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                B pc_est= reg[l_index][0];
                if ((phi - last_statu[0]) >= p * 3) {
                    for (int i = phi - 2*p; i > last_statu[0]; i -= 2 * p)
                        {   pc_est = pc_est ^ u[l_index][i];}
                }

                if (pc[phi] == pc_est) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                }
                reg[l_index] = CSR(reg[l_index]);
            }    

            last_statu[0] = phi;
        }
        last_statu = CSR(last_statu);
        cnt_u = cnt_u + 1;
      
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc4(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& pc2,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    int p2,
    uint32_t crc,
    int L)
{
    std::vector<B>                      phi_status(p, 0);
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));
    std::vector<std::vector<B>>         reg2(L, std::vector<B>(p2, 0));

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    reg2[l_index][0] = reg2[l_index][0] == 1 ? 0 : 1;
                    reg2[l_index] = CSR(reg2[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg2[l_index] = CSR(reg2[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    reg2[index] = reg2[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg2[l_index] = CSR(reg2[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    reg2[index][0] = reg2[index][0] == 1 ? 0 : 1;
                    reg2[index] = CSR(reg2[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;
                if (phi_status[phi%p]==0) {
                    if (pc[phi] == reg[l_index][0]) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg2[l_index][0] = reg2[l_index][0] == 1 ? 0 : 1;
                    }
                }
                else {
                    if (pc2[phi] == reg2[l_index][0]) {
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                        reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                        reg2[l_index][0] = reg2[l_index][0] == 1 ? 0 : 1;
                    }
                }
                reg[l_index] = CSR(reg[l_index]);
                reg2[l_index] = CSR(reg2[l_index]);
            }
        }
        phi_status[phi % p] = 1 - phi_status[phi % p];
        cnt_u = cnt_u + 1;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }
            }
        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}




template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepc2(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    std::vector<std::vector<int>>  &  pcindex,
    int p,
    uint32_t crc,
    int L,
    const size_t frame_id)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<std::vector<B>>         PCest(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    int pcnumber = 0;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1= lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                            P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<size_t> PM_inx(2 * L);
            for (size_t i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](size_t a, size_t b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (size_t i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                for (int j = 1; j < pcindex[pcnumber].size(); j++) {
                     PCest[l_index][phi] = PCest[l_index][phi] ^ u[l_index][pcindex[pcnumber][j]];
                 }

                if (pc[phi] == PCest[l_index][phi]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                }
                else{
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                }

            }
            pcnumber++;
        }
        cnt_u = cnt_u + 1;

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
                dec_bit[i] = polar_info_esti[j];
                j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepcsegment(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    const std::vector<uint32_t>& crc_v,
    std::vector<int >& seg_s,
    int L,
    const size_t frame_id)
{
    bool right = false;
    int crcindex = N / std::pow(2, seg_s.size() - 1);
    int crc_v_dex = 0;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));

    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }
    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;
    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                //if (regsite_cur != (regsite_last + 1) % p) {
                //    for (int x = 0; x < regsite_last + p + 1 - regsite_cur; x++)
                //        reg[l_index] = CSR(reg[l_index]);
                //}

                if (pc[phi] == reg[l_index][0]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    reg[l_index] = CSR(reg[l_index]);
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                }

            }
            //regsite_last = regsite_cur;
        }
        cnt_u = cnt_u + 1;
        //regsite_cur = (regsite_cur + 1) % p;
        int   countL = 0;
        if (phi == crcindex-1) {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (crc_v[crc_v_dex] != calculateCRC32(u[l_index], crcindex)) {
                    activepath[l_index] = 0;
                    countL++;
                }
            }
            crcindex *= 2;
            crc_v_dex++;
        }
        if (countL == 16) {
            seg_s[crc_v_dex-1]=1;
            break;
        }
          
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }
        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodescl(const std::vector<R>& llr,
    const std::vector<B>& frozen_bits_value, 
    const std::vector<bool>& fsite, 
    std::vector<B>& dec_bit, 
    uint32_t crc, int L, const size_t frame_id)
{
    bool right = false;
    dec_bit = frozen_bits_value;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L, 0.0));
    std::vector<std::vector<B>>       C(N - 1, std::vector<B>(2 * L, 0));
    std::vector<std::vector<B>>       u(L, std::vector<B>(K, 0));
    std::vector<R>                      PM(L, 0.0);
    std::vector<int>                    activepath(L, 0);
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][lazy_copy[layer + 1][l_index]] + P[beta + index_1 + index_2 - 1][lazy_copy[layer + 1][l_index]];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        //if now we decode an unfrozen bit
        if (fsite[phi] == 0) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }

            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<size_t> PM_inx(2 * L);
            for (size_t i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](size_t a, size_t b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (size_t i = 0; i < middle; ++i) {
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    //u[index] = u[l_index];

                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    break;
                }
                }
            }
            cnt_u = cnt_u + 1;
        }
        //frozen bit operation
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                if (P[0][l_index] < 0 && frozen_bits_value[phi] == 0) {
                    PM[l_index] = PM[l_index] - P[0][l_index];
                }
                if (P[0][l_index] > 0 && frozen_bits_value[phi] == 1) {
                    PM[l_index] = PM[l_index] + P[0][l_index];
                }
                if (phi_mod_2 == 0) {
                    C[0][2 * l_index] = frozen_bits_value[phi];
                }
                else {
                    C[0][2 * l_index + 1] = frozen_bits_value[phi];
                }
            }
        }

        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }


        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(K);
        for (int i = 0; i < K; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            if (fsite[i] == 0) {
                dec_bit[i] = polar_info_esti[j];
                j++;
            }
        }
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            break;
        }
        /*else {
            j = 0;
            for (int i = 0; i < N; ++i) {
                if (frozen_bits[i] == 0) {
                    polar_info_esti_all[i] = polar_info_esti[j];
                    j++;
                }
            }
        }*/
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodescltest(const std::vector<R>& llr,
    const std::vector<B>& frozen_bits_value,
    const std::vector<bool>& fsite,
    std::vector<B>& dec_bit,
    uint32_t crc, int L, 
    std::vector<std::vector<int>>& error_site)
{
    bool right = false;
    dec_bit = frozen_bits_value;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L, 0.0));
    std::vector<std::vector<B>>       C(N - 1, std::vector<B>(2 * L, 0));
    std::vector<std::vector<B>>       u(L, std::vector<B>(K, 0));
    std::vector<R>                      PM(L, 0.0);
    std::vector<int>                    activepath(L, 0);
    std::vector<int>                    is_error(L, 0);
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;
    int cur_phi = 0;
    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][lazy_copy[layer + 1][l_index]] + P[beta + index_1 + index_2 - 1][lazy_copy[layer + 1][l_index]];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        //if now we decode an unfrozen bit
        if (fsite[phi] == 0) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }

            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<size_t> PM_inx(2 * L);
            for (size_t i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](size_t a, size_t b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (size_t i = 0; i < middle; ++i) {
                int row = PM_inx[i] / L;
                int col = PM_inx[i] % L;
                compare[row][col] = 1;
            }

            std::vector<int> kill_index;
            int kill_cnt = 0;
            for (int i = 0; i < L; ++i) {
                if (compare[0][i] == 0 && compare[1][i] == 0) {
                    activepath[i] = 0;
                    is_error[i] = 0;
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    if (frozen_bits_value[phi] != 1) { 
                        is_error[l_index] = 1; }
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    if (frozen_bits_value[phi] != 0) { 
                        is_error[l_index] = 1; }
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());

                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    if (frozen_bits_value[phi] != 0) { 
                        is_error[l_index] = 1; }
                    else { 
                        is_error[index] = 1; }
                    break;
                }
                }
            }
            cnt_u = cnt_u + 1;
        }
        //frozen bit operation
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                if (P[0][l_index] < 0 && frozen_bits_value[phi] == 0) {
                    PM[l_index] = PM[l_index] - P[0][l_index];
                }
                if (P[0][l_index] > 0 && frozen_bits_value[phi] == 1) {
                    PM[l_index] = PM[l_index] + P[0][l_index];
                }
                if (phi_mod_2 == 0) {
                    C[0][2 * l_index] = frozen_bits_value[phi];
                }
                else {
                    C[0][2 * l_index + 1] = frozen_bits_value[phi];
                }
            }
        }
        int sum_error = 0;
        for (int num_error : is_error) {
            sum_error += num_error;
        }
        cur_phi = phi;
        if (sum_error == 16) {
            break;
        }
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }


        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(K);
        for (int i = 0; i < K; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            if (fsite[i] == 0) {
                dec_bit[i] = polar_info_esti[j];
                j++;
            }
        }
        
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            for (int i = 0; i < l_index; i++) { error_site.pop_back(); }
            
            break;
        }
        //if(l_index==0){
        std::vector<int> error_site1;
        for (size_t i = 0; i <= cur_phi; ++i) {
            if (dec_bit[i] != frozen_bits_value[i]) {
                error_site1.push_back(i); 
            }
        }
        error_site.push_back(error_site1);
        //}
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepctest(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& frozen_bits_value,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    std::vector<std::vector<int>>& error_site)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));
    std::vector<int>                    is_error(L, 0);

    int cur_phi = 0;
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
                int row = PM_inx[i] / L;
                int col = PM_inx[i] % L;
                compare[row][col] = 1;
            }

            std::vector<int> kill_index;
            int kill_cnt = 0;
            for (int i = 0; i < L; ++i) {
                if (compare[0][i] == 0 && compare[1][i] == 0) {
                    activepath[i] = 0;
                    is_error[i] = 0;
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 1) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    else {
                        is_error[index] = 1;
                    }
                    break;
                }
                }
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                //if (regsite_cur != (regsite_last + 1) % p) {
                //    for (int x = 0; x < regsite_last + p + 1 - regsite_cur; x++)
                //        reg[l_index] = CSR(reg[l_index]);
                //}

                if (pc[phi] == reg[l_index][0]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    reg[l_index] = CSR(reg[l_index]);
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                }

            }
            //regsite_last = regsite_cur;
        }
        cnt_u = cnt_u + 1;
        //regsite_cur = (regsite_cur + 1) % p;
        int sum_error = 0;
        for (int num_error : is_error) {
            sum_error += num_error;
        }
        cur_phi = phi;
        if (sum_error == 16) {
            break;
        }
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }
        
        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            for (int i = 0; i < l_index; i++) { error_site.pop_back(); }
            break;
        }
        //if (l_index == 0) {
        std::vector<int> error_site1;
        for (size_t i = 0; i <= cur_phi; ++i) {
            if (dec_bit[i] != frozen_bits_value[i]) {
                error_site1.push_back(i) ;
            }
        }
        error_site.push_back(error_site1);
        //}
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepctest2(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& frozen_bits_value,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    std::vector<std::vector<int>>& error_site)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<std::vector<B>>         reg(L, std::vector<B>(p, 0));
    std::vector<int>                    is_error(L, 0);

    int cur_phi = 0;
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
                int row = PM_inx[i] / L;
                int col = PM_inx[i] % L;
                compare[row][col] = 1;
            }

            std::vector<int> kill_index;
            int kill_cnt = 0;
            for (int i = 0; i < L; ++i) {
                if (compare[0][i] == 0 && compare[1][i] == 0) {
                    activepath[i] = 0;
                    is_error[i] = 0;
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 1) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    else {
                        is_error[index] = 1;
                    }
                    break;
                }
                }
                
            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                if (pc[phi] == reg[l_index][0]) {
                    PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                }
                else {
                    PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                }
                reg[l_index] = CSR(reg[l_index]);
            }
        }
        cnt_u = cnt_u + 1;
        
        int sum_error = 0;
        for (int num_error : is_error) {
            sum_error += num_error;
        }
        cur_phi = phi;
        if (sum_error == 16) {
            break;
        }
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }

        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            for (int i = 0; i < l_index; i++) { error_site.pop_back(); }
            break;
        }
        //if (l_index == 0) {
        std::vector<int> error_site1;
        for (size_t i = 0; i <= cur_phi; ++i) {
            if (dec_bit[i] != frozen_bits_value[i]) {
                error_site1.push_back(i);
            }
        }
        error_site.push_back(error_site1);
        //}
    }
    return right;
}

template<typename B, typename R>
bool Decoder_polar<B, R>
::decodepctest3(const std::vector<R>& llr,
    const std::vector<B>& pc,
    const std::vector<B>& frozen_bits_value,
    const std::vector<bool>& pcsite,
    std::vector<B>& dec_bit,
    int p,
    uint32_t crc,
    int L,
    std::vector<std::vector<int>>& error_site)
{
    bool right = false;
    dec_bit = pc;
    int m = log2(N);
    int cnt_u = 0;
    std::vector<std::vector<int>>       lazy_copy(m, std::vector<int>(L, 0));
    std::vector<std::vector<R>>         P(N - 1, std::vector<R>(L));
    std::vector<std::vector<B>>         C(N - 1, std::vector<int>(2 * L));
    std::vector<std::vector<B>>         u(L, std::vector<int>(N, 0));
    std::vector<R>                      PM(L, 0.0);//初始值必须为0
    std::vector<int>                    activepath(L, 0);
    std::vector<int>                    is_error(L, 0);

    int cur_phi = 0;
    activepath[0] = 1;
    for (int i = 0; i < m; ++i) {
        lazy_copy[i][0] = 0;
    }

    int layer;
    int phi_mod_2;
    int index_1;
    int index_2;

    std::vector<std::vector<B>>         reg(L, std::vector<B>(p*2, 0));
    std::vector<std::vector<B>>         reglast(L, std::vector<B>(p*2, 0));
    std::vector<int> _pc(p, 0);

    for (int phi = 0; phi < N; ++phi) {
        layer = llr_layer_vec[phi];
        phi_mod_2 = phi % 2;
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;

            if (phi == 0) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = f(llr[beta], llr[beta + index_1]);
                }

                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
            else if (phi == N / 2) {
                index_1 = lambda_offset[m - 1];
                for (int beta = 0; beta < index_1; ++beta) {
                    int x_tmp = C[beta + index_1 - 1][2 * l_index];
                    P[beta + index_1 - 1][l_index] = (1 - 2 * x_tmp) * llr[beta] + llr[beta + index_1];
                }
                for (int i_layer = m - 2; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }

            }
            else {//----------------------
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                int L1 = lazy_copy[layer + 1][l_index];
                for (int beta = 0; beta < index_1; ++beta) {
                    P[beta + index_1 - 1][l_index] = (1 - 2 * C[beta + index_1 - 1][2 * l_index]) *
                        P[beta + index_2 - 1][L1] + P[beta + index_1 + index_2 - 1][L1];
                }
                for (int i_layer = layer - 1; i_layer >= 0; --i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = 0; beta < index_1; ++beta) {
                        P[beta + index_1 - 1][l_index] = f(P[beta + index_2 - 1][l_index], P[beta + index_1 + index_2 - 1][l_index]);
                    }
                }
            }
        }

        if (pcsite[phi] == false) {
            std::vector<std::vector<R>> PM_pair(2, std::vector<R>(L, std::numeric_limits<R>::max()));
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0) {
                    continue;
                }
                if (P[0][l_index] >= 0) {
                    PM_pair[0][l_index] = PM[l_index];
                    PM_pair[1][l_index] = PM[l_index] + P[0][l_index];
                }
                else {
                    PM_pair[0][l_index] = PM[l_index] - P[0][l_index];
                    PM_pair[1][l_index] = PM[l_index];
                }
            }
            int middle = std::min(2 * std::accumulate(activepath.begin(), activepath.end(), 0), L);
            std::vector<int> PM_inx(2 * L);
            for (int i = 0; i < PM_inx.size(); ++i) {
                PM_inx[i] = i;
            }
            std::sort(PM_inx.begin(), PM_inx.end(), [&](int a, int b) {
                return PM_pair[a / L][a % L] < PM_pair[b / L][b % L];
                });
            std::vector<std::vector<int>> compare(2, std::vector<int>(L, 0));
            for (int i = 0; i < middle; ++i) {
                int row = PM_inx[i] / L;
                int col = PM_inx[i] % L;
                compare[row][col] = 1;
            }

            std::vector<int> kill_index;
            int kill_cnt = 0;
            for (int i = 0; i < L; ++i) {
                if (compare[0][i] == 0 && compare[1][i] == 0) {
                    activepath[i] = 0;
                    is_error[i] = 0;
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
                    u[l_index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 1;
                    PM[l_index] = PM_pair[1][l_index];
                    reg[l_index][0] = reg[l_index][0] == 1 ? 0 : 1;
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 1) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 2:
                    u[l_index][cnt_u] = 0;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    PM[l_index] = PM_pair[0][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    break;
                case 3: {
                    int index = kill_index[kill_cnt - 1];
                    --kill_cnt; // pop stack
                    activepath[index] = 1;
                    // Lazy copy
                    for (int i = 0; i < lazy_copy.size(); ++i) {
                        lazy_copy[i][index] = lazy_copy[i][l_index];
                    }
                    //u[index] = u[l_index];
                    std::copy(u[l_index].begin(), u[l_index].begin() + cnt_u, u[index].begin());
                    u[l_index][cnt_u] = 0;
                    u[index][cnt_u] = 1;
                    C[0][2 * l_index + phi_mod_2] = 0;
                    C[0][2 * index + phi_mod_2] = 1;
                    reg[index] = reg[l_index];
                    reglast[index] = reglast[l_index];
                    PM[l_index] = PM_pair[0][l_index];
                    PM[index] = PM_pair[1][l_index];
                    reg[l_index] = CSR(reg[l_index]);
                    reg[index][0] = reg[index][0] == 1 ? 0 : 1;
                    reg[index] = CSR(reg[index]);
                    if (frozen_bits_value[phi] != 0) {
                        is_error[l_index] = 1;
                    }
                    else {
                        is_error[index] = 1;
                    }
                    break;
                }
                }

            }
        }
        else {
            for (int l_index = 0; l_index < L; ++l_index) {
                if (activepath[l_index] == 0)
                    continue;

                if (_pc[phi % (p)] == 0) {
                    if (pc[phi] == (reglast[l_index][phi % (p)] ^ reg[l_index][0] ^ reg[l_index][p])) {
                        //
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                    }
                    
                }
                else {
                    if (pc[phi] == (reglast[l_index][phi % (p)+p] ^ reg[l_index][0] ^ reg[l_index][p])) {
                        //
                        PM[l_index] = PM[l_index] - (P[0][l_index] < 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 0;
                        C[0][2 * l_index + phi_mod_2] = 0;
                    }
                    else {
                        PM[l_index] = PM[l_index] + (P[0][l_index] > 0 ? P[0][l_index] : 0);
                        u[l_index][cnt_u] = 1;
                        C[0][2 * l_index + phi_mod_2] = 1;
                    }
                    
                }
                reglast[l_index][phi % (p)] = reg[l_index][0];
                reglast[l_index][phi % (p)+p] = reg[l_index][p];
                reg[l_index] = CSR(reg[l_index]);
            }
            _pc[phi%p] = 1 - _pc[phi%p];
        }
        cnt_u = cnt_u + 1;

        int sum_error = 0;
        for (int num_error : is_error) {
            sum_error += num_error;
        }
        cur_phi = phi;
        if (sum_error == 16) {
            break;
        }
        //partial-sum return
        for (int l_index = 0; l_index < L; ++l_index) {
            if (activepath[l_index] == 0)
                continue;
            if ((phi_mod_2 == 1) && (phi != N - 1)) {
                layer = bit_layer_vec[phi];
                for (int i_layer = 0; i_layer < layer; ++i_layer) {
                    index_1 = lambda_offset[i_layer];
                    index_2 = lambda_offset[i_layer + 1];
                    for (int beta = index_1; beta < 2 * index_1; ++beta) {
                        C[beta + index_1 - 1][2 * l_index + 1] = (C[beta - 1][2 * lazy_copy[i_layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                        C[beta + index_2 - 1][2 * l_index + 1] = C[beta - 1][2 * l_index + 1];
                    }
                }
                index_1 = lambda_offset[layer];
                index_2 = lambda_offset[layer + 1];
                for (int beta = index_1; beta < 2 * index_1; ++beta) {
                    C[beta + index_1 - 1][2 * l_index] = (C[beta - 1][2 * lazy_copy[layer][l_index]] + C[beta - 1][2 * l_index + 1]) % 2; // 左列延迟复制
                    C[beta + index_2 - 1][2 * l_index] = C[beta - 1][2 * l_index + 1];
                }

            }

        }

        //lazy_copy----------------
        if (phi < N - 1) {
            for (int i_layer = 0; i_layer <= llr_layer_vec[phi + 1]; ++i_layer) {
                for (int l_index = 0; l_index < L; ++l_index) {
                    lazy_copy[i_layer][l_index] = l_index;
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
    // Perform path selection and decode
    for (int l_index = 0; l_index < L; ++l_index) {
        int path_num = path_ordered[l_index];
        std::vector<B> polar_info_esti(N);
        for (int i = 0; i < N; ++i) {
            polar_info_esti[i] = u[path_num][i]; // Extract the specific column
        }
        int j = 0;
        for (int i = 0; i < N; ++i) {
            dec_bit[i] = polar_info_esti[j];
            j++;
        }

        if (crc == calculateCRC32(dec_bit, dec_bit.size())) {
            right = true;
            for (int i = 0; i < l_index; i++) { error_site.pop_back(); }
            break;
        }
        //if (l_index == 0) {
        std::vector<int> error_site1;
        for (size_t i = 0; i <= cur_phi; ++i) {
            if (dec_bit[i] != frozen_bits_value[i]) {
                error_site1.push_back(i);
            }
        }
        error_site.push_back(error_site1);
        //}
    }
    return right;
}


int oat_hash( int h,int key)
{
    h += key;
    h += (h << 10);
    h ^= (h >> 6);
    h += (h << 3);
    h ^= (h >> 11);
    h += (h << 15);

    return h;
}


template <typename B, typename R>
void Decoder_polar<B, R>
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
void Decoder_polar<B, R>
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
void Decoder_polar<B, R>
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
R Decoder_polar<B, R>
::f(R L1, R L2)
{
    return std::copysign(1.0, L1) * std::copysign(1.0, L2) * std::min(std::fabs(L1), std::fabs(L2));
}

#endif // DECODER_POLAR2_HPP_