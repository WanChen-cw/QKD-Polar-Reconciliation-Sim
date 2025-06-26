/*!
 * \file
 * \brief Class module::CRC_polynomial.
 */
#ifndef MYCRC_HPP_
#define MYCRC_HPP_

#include <vector>
#include <string>
template <typename B = int>
class Mycrc
{
private:
	unsigned	polynomial_packed;
	int			size;
	int			K;
	std::vector<B> polynomial;
	std::vector<B> buff_crc;
public:
	Mycrc(const int& K, const std::string& poly_key="0x04C11DB7", const int& size = 32);
	~Mycrc();
	
	void build(const B* U_K1, B* U_K2, const size_t frame_id);
	void generate(const B* U_in,B* U_out,const int off_in,const int off_out,const int loop_size);
	void get_crc(const B* U_in,  B* Ucrc);
	void get_crc(const std::vector<B>& U_in, std::vector<B>& Ucrc);
};

template<typename B>
Mycrc<B>
::Mycrc(const int& K, const std::string& poly_key, const int& size)
	:K(K), size(size), polynomial_packed((unsigned)std::stoul(poly_key, 0, 16)), polynomial(0), buff_crc(K)
{
	polynomial.push_back(1);
	for (auto i = 0; i < this->size; i++)
	polynomial.push_back((polynomial_packed >> ((this->size - 1) - i)) & 1);
}

template<typename B>
Mycrc<B>
::~Mycrc()
{
}

template<typename B>
void Mycrc<B>
::build(const B* U_K1, B* U_K2, const size_t frame_id)
{
		std::copy(U_K1, U_K1 + this->K, U_K2);
		this->_generate(U_K1, U_K2, 0, this->K, this->K);
}

template<typename B>
void Mycrc<B>
::generate(const B* U_in, B* U_out, const int off_in, const int off_out, const int loop_size)
{
	std::copy(U_in + off_in, U_in + off_in + loop_size, buff_crc.begin());
	std::fill(buff_crc.begin() + loop_size, buff_crc.begin() + loop_size + this->size, (B)0);
	for (auto i = 0; i < loop_size; i++)
		if (buff_crc[i])
			for (auto j = 0; j <= this->size; j++)
				if (this->polynomial[j])
					buff_crc[i + j] = !buff_crc[i + j];

	if (U_out != buff_crc.data())
		std::copy(buff_crc.begin() + loop_size, buff_crc.begin() + loop_size + this->size, U_out + off_out);
}

template<typename B>
void Mycrc<B>
::get_crc(const B* U_in,  B* Ucrc)
{
	std::copy(U_in , U_in+ this->K, buff_crc.begin());
	std::fill(buff_crc.begin() + this->K, buff_crc.begin() + this->K + this->size, (B)0);
	for (auto i = 0; i < this->K; i++)
		if (buff_crc[i])
			for (auto j = 0; j <= this->size; j++)
				if (this->polynomial[j])
					buff_crc[i + j] = !buff_crc[i + j];
	std::copy(buff_crc.begin() + this->K, buff_crc.begin() + this->K + this->size, Ucrc);
}

template<typename B>
void Mycrc<B>
::get_crc(const std::vector<B>& U_in, std::vector<B>& Ucrc)
{
	this->buff_crc.resize(this->K+ this->size);
	std::copy(U_in.begin(), U_in.end(), buff_crc.begin());
	std::fill(buff_crc.begin() + this->K, buff_crc.begin() + this->K + this->size, (B)0);
	for (auto i = 0; i < this->K; i++)
		if (buff_crc[i])
			for (auto j = 0; j <= this->size; j++)
				if (this->polynomial[j])
					buff_crc[i + j] = !buff_crc[i + j];
	std::copy(buff_crc.begin() + this->K, buff_crc.begin() + this->K + this->size, Ucrc.begin());
}


uint32_t calculateCRC32(const std::vector<int>&  data, size_t size) {
	uint32_t crc = 0xFFFFFFFF;
	uint32_t crcPolynomial = 0x04C11DB7;//0x814141AB 0x04C11DB7 320x1EDC6F41 0x814141AB 0x32583499
	for (size_t i = 0; i < size; ++i) {
		crc ^= (data[i] << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
	}
	return crc;
}

void calculateCRC32(uint32_t &crc,const std::vector<int>& data, size_t size) {
	uint32_t crcPolynomial = 0x04C11DB7;//0x814141AB 0x04C11DB7 320x1EDC6F41 0x814141AB 0x32583499
	for (size_t i = 0; i < size; ++i) {
		crc ^= (data[i] << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
	}
}

void calculateCRC32(uint32_t& crc, const int  data, size_t size) {
	uint32_t crcPolynomial = 0x04C11DB7;//0x814141AB 0x04C11DB7 320x1EDC6F41 0x814141AB 0x32583499
	for (size_t i = 0; i < size; ++i) {
		crc ^= (data << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
	}
}

std::vector<uint32_t> calculateCRC32(const std::vector<int>& data, size_t size,int n_seg) {
	std::vector<uint32_t> crc_v;
	int c= size/std::pow(2, n_seg-1);
	uint32_t crc = 0xFFFFFFFF;
	uint32_t crcPolynomial = 0x04C11DB7;//0x814141AB 0x04C11DB7 320x1EDC6F41 0x814141AB 0x32583499
	for (size_t i = 0; i < size; ++i) {
		crc ^= (data[i] << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
		if (i == c - 1) {
			crc_v.push_back(crc);
			c *= 2;
		}
	}
	return crc_v;
}
#endif /*MYCRC_HPP_ */


