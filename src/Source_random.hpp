
#ifndef SOURCE_RANDOM_HPP_
#define SOURCE_RANDOM_HPP_

#include <random>

template <typename B = int>
class Source_random 
{
private:
	int K;
	std::random_device rd;
	std::mt19937 rd_engine; // Mersenne Twister 19937
	// std::minstd_rand rd_engine; // LCG
#ifdef _MSC_VER
	std::uniform_int_distribution<short> uniform_dist;
#else
	std::uniform_int_distribution<B> uniform_dist;
#endif

public:
	Source_random(const int K, const int seed = 17);

	~Source_random();

	void generate(std::vector<B>& U_K, const size_t frame_id = 0);

protected:
			
};


template<typename B>
Source_random<B>
::Source_random(const int K, const int seed)
	:K(K), rd(), rd_engine(rd()+ seed), uniform_dist(0, 1)
{
}

template<typename B>
Source_random<B>
::~Source_random()
{
}

template<typename B>
void Source_random<B>
::generate(std::vector<B>& U_K, const size_t frame_id)
{
	for (auto i = 0; i < this->K; i++)
		U_K[i] = (B)this->uniform_dist(this->rd_engine);
}

#endif /* SOURCE_RANDOM_HPP_ */