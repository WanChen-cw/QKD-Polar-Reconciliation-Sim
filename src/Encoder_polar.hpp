/*!
 * \file
 * \brief Class module::Encoder_polar.
 */
#ifndef ENCODER_POLAR_HPP_
#define ENCODER_POLAR_HPP_
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <sstream>
#include <algorithm>

template <typename B = int>
class Encoder_polar
{
private:
	std::vector<bool>  frozen_bits; // true means frozen, false means set to 0/1
	std::vector<B>     X_N_tmp;
	int N;
	int K;
	bool Q;//

public:
	Encoder_polar(const int& K, const int& N,const std::vector<bool>& frozen_bits,const bool Q=0 );
	~Encoder_polar();

	const std::vector<bool>& get_frozen_bits() const;
	void set_frozen_bits(const std::vector<bool>& frozen_bits);

	void encode( B* U, B* X);
	void encode(std::vector<B> &U, std::vector<B> &X);
	void encode(std::vector<B>& U);
	bool is_codeword(const B* X_N);

protected:
	void convert(const B* U_K, B* U_N);
	void convert(std::vector<B> &U_K, std::vector<B> &U_N);
	void deconvert(const B* U_K, B* U_N);
	void deconvert(std::vector<B> &U_K, std::vector<B> &U_N);
	void light_encode(B* bits);
	void light_encode(std::vector<B> &bits);
};
	




template <typename B>
Encoder_polar<B>
::Encoder_polar(const int& K, const int& N, const std::vector<bool>& frozen_bits, const bool Q)
	: K(K), N(N), Q(Q), frozen_bits(N)
{
	
	this->set_frozen_bits(frozen_bits);
}

template<typename B>
Encoder_polar<B>::~Encoder_polar()
{
}


template <typename B>
void Encoder_polar<B>
::encode(B* U, B* X)
{
	if (this->Q == 0)
	{
		this->convert(U, X);
		this->light_encode(X);
	}
	else
	{
		this->light_encode(U);

		this->deconvert(U, X);
	}
}

template<typename B>
void Encoder_polar<B>
::encode(std::vector<B>& U, std::vector<B>& X)
{
	if (this->Q == 0)
	{
		this->convert(U, X);
		this->light_encode(X);
	}
	else
	{
		this->light_encode(U);
		this->deconvert(U, X);
	}
}

template<typename B>
void Encoder_polar<B>
::encode(std::vector<B>& U){
	this->light_encode(U);
}

template <typename B>
void Encoder_polar<B>
::light_encode(B* bits)
{
	for (auto k = (this->N >> 1); k > 0; k >>= 1)
		for (auto j = 0; j < this->N; j += 2 * k)
			for (auto i = 0; i < k; i++)
				bits[j + i] = bits[j + i] ^ bits[k + j + i];
}

template<typename B>
void Encoder_polar<B>::
light_encode(std::vector<B>&  bits)
{
	for (auto k = (this->N >> 1); k > 0; k >>= 1)
		for (auto j = 0; j < this->N; j += 2 * k)
			for (auto i = 0; i < k; i++)
				bits[j + i] = bits[j + i] ^ bits[k + j + i];
}

template <typename B>
void Encoder_polar<B>
::convert(const B* U_K, B* U_N)
{
	if (U_K == U_N)
	{
		std::vector<B> U_K_tmp(this->K);
		std::copy(U_K, U_K + this->K, U_K_tmp.begin());

		auto j = 0;
		for (unsigned i = 0; i < frozen_bits.size(); i++)
			U_N[i] = (frozen_bits[i]) ? (B)0 : U_K_tmp[j++];
	}
	else
	{
		auto j = 0;
		for (unsigned i = 0; i < frozen_bits.size(); i++)
			U_N[i] = (frozen_bits[i]) ? (B)0 : U_K[j++];
	}
}

template<typename B>
void Encoder_polar<B>::
convert(std::vector<B> &U_K, std::vector<B> &U_N)
{
	if (U_K == U_N)
	{
		std::vector<B> U_K_tmp(this->K);
		std::copy(U_K.begin(), U_K.end(), U_K_tmp.begin());

		auto j = 0;
		for (unsigned i = 0; i < frozen_bits.size(); i++)
			U_N[i] = (frozen_bits[i]) ? (B)0 : U_K_tmp[j++];
	}
	else
	{
		auto j = 0;
		for (unsigned i = 0; i < frozen_bits.size(); i++)
			U_N[i] = (frozen_bits[i]) ? (B)0 : U_K[j++];
	}
}

template<typename B>
void Encoder_polar<B>
::deconvert(const B* U_K, B* U_N)
{
	for (unsigned i = 0; i < frozen_bits.size(); i++)
		U_N[i] = (frozen_bits[i]) ? (B)U_K[i] :(B) 2;
}

template<typename B>
void Encoder_polar<B>::
deconvert(std::vector<B> &U_K, std::vector<B> &U_N)
{
	for (unsigned i = 0; i < frozen_bits.size(); i++)
		U_N[i] = (frozen_bits[i]) ? (B)U_K[i] : (B)2;
}

template <typename B>
bool Encoder_polar<B>
::is_codeword(const B* X_N)
{
	std::copy(X_N, X_N + this->N, this->X_N_tmp.data());

	for (auto k = (this->N >> 1); k > 0; k >>= 1)
		for (auto j = 0; j < this->N; j += 2 * k)
		{
			for (auto i = 0; i < k; i++)
				this->X_N_tmp[j + i] = this->X_N_tmp[j + i] ^ this->X_N_tmp[k + j + i];

			if (this->frozen_bits[j + k - 1] && this->X_N_tmp[j + k - 1])
				return false;
		}

	return true;
}

template <typename B>
void Encoder_polar<B>
::set_frozen_bits(const std::vector<bool>& frozen_bits)
{
	std::copy(frozen_bits.begin(), frozen_bits.end(), this->frozen_bits.begin());
}


template <typename B>
const std::vector<bool>& Encoder_polar<B>
::get_frozen_bits() const
{
	return this->frozen_bits;
}
#endif // ENCODER_POLAR_HPP_