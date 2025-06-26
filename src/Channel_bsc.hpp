
#ifndef CHANNEL_BSC_HPP_
#define CHANNEL_BSC_HPP_

#include <random>
#include <cmath>

template <typename B = int, typename R = float>
class Channel_bsc
{
private:
	float llr0;
	float llr1;
	float ep;
	std::random_device rd;             // 随机种子生成器
	std::mt19937 gen;                  // 随机数生成器引擎
	std::uniform_int_distribution<int> distribution;

public:
	Channel_bsc(const float ep ,const int seed = 0);
	~Channel_bsc();
	void add_noise(const std::vector<B>& U, std::vector<R>& llr);
	void bittollr(const std::vector<B>& U, std::vector<R>& llr);

protected:
};
template<typename B, typename R>
Channel_bsc<B, R>
::Channel_bsc(const float ep, const int seed)
	:ep(ep), rd(), gen(rd()), distribution(1, 10000)
{
	llr0 =log2((1 - ep) / ep);
	llr1 =log2(ep /(1- ep));
}

template<typename B, typename R>
 Channel_bsc<B, R>
::~Channel_bsc()
{
}

 template<typename B, typename R>
 inline void Channel_bsc<B, R>
::add_noise(const std::vector<B>& U, std::vector<R>& llr)
 {
	  //B tmp= bern_dist(this->rd_engine);
	 for (unsigned i = 0; i < llr.size(); i++){
		 // 生成随机数
		 int random_number = distribution(gen);
		 B tmp = random_number > 10000*ep ? 0 : 1;
		 llr[i] = (U[i]==tmp)?llr0:llr1;
	 }
 }

 template<typename B, typename R>
 inline void Channel_bsc<B, R>
::bittollr(const std::vector<B>& U, std::vector<R>& llr)
 {
	 for (unsigned i = 0; i < llr.size(); i++) {
	 llr[i] = U[i]? llr1 : llr0;
	 }
 }
#endif /* CHANNEL_BSC_HPP_ */


