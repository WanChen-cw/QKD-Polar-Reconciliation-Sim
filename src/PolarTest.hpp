#ifndef TESTPOLAR_HPP_
#define TESTPOLAR_HPP_
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <fstream>

#include "Source_random.hpp"
#include "Frozenbits_generator.hpp"
#include "Encoder_polar.hpp"
#include "Channel_bsc.hpp"
#include "Decoder_polar2.hpp"      //对应encode2  3
#include "Decoder_polar_fast.hpp" 
 //#include "Decoder_polar.hpp"    //对应encode0  1
#include "Mycrc.hpp"

void saveToTxt(const std::vector<std::vector<int>>& data, const char* filename) {
	std::ofstream txtFile(filename);
	if (!txtFile.is_open()) {
		std::cerr << "Error creating text file" << std::endl;
		return;
	}

	for (const auto& row : data) {
		for (const auto& val : row) {
			txtFile << val << " ";
		}
		txtFile << std::endl;
	}

	txtFile.close();
}


void fastpolartest()
{
	int N = 65536;
	for (int N2 = 20350; N2 <= 22500; N2 += 50) {
		int K = N - N2;
		float ep = 0.05;
		std::string path = "conf/N_65536_005.pc";
		std::vector<bool>	fp(N);
		std::vector<int>	u(N);
		std::vector<float>  llr(N);
		std::vector<int>	decbuff(N);

		std::unique_ptr<Source_random<>>		source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
		std::unique_ptr<Frozenbits_generator>	fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
		fb_generate->generate(fp);
		std::unique_ptr<Encoder_polar<>>		encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
		std::unique_ptr<Channel_bsc<>>			channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 1));
		Decoder_polar_fast<> decoder(N);

		decoder.get_node_structure(fp);
		decoder.get_psi_for_advanced_sc_decoder();

		int t = 100;
		volatile int a = 0;
		volatile int count1 = 0;
		while (count1 < t) {
			a++;
			source->generate(u);
			channel->add_noise(u, llr);
			encoder->encode(u);
			uint32_t ccc = calculateCRC32(u, u.size());
			//if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, 16)) {
			//	count1++;
			//}
			if (!decoder.fastscl(llr, u, decbuff, ccc, 16)) {
				count1++;
			}
		}
		std::cout << std::endl << "number_fbit: " << N2
			<< "  frame_number: " << a
			<< "  errorpc_number: " << count1 
			<< "  fer1:  " << (float)count1 / a;
	}
}



void polartest1()
{
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 200;
	int K = N-N2;
	const int p =9;
	float ep = 0.02;
	int L = 16;
	std::string path = "conf/N_1024_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits
	
	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 1000; 
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	std::cout << std::endl << "start" << std::endl;
	while (count2 < t){
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		//std::vector<int> reg2(p, h);
		//for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
		//	reg2[0] = oat_hash(reg2[0] , u[encodepc_loop]); // XOR operation
		//	PC2[encodepc_loop] = reg2[0] >0;
		//	reg2 = CSR(reg2); // Perform cyclic shift right
		//}  //hash


		//std::vector<int> reg(p, 0);
		//for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
		//	PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
		//	reg[0] = PC[encodepc_loop];
		//	reg = CSR(reg); // Perform cyclic shift right
		//}  //huawei 5G
		//if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, L)) {
		//	count1++;
		//}
		if (!decoder->decodescl(llr, u, fp,decbuff2, ccc, L)){
			count2++;
		}
	}
	std::cout << std::endl
		<< "frame_number: " <<a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorscl: " << count2 << " fer2:  " 
		<< (float)count2 / a
		<< std::endl ;
	 
}

void repolartest()
{
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 5950;
	int blength = 350;
	int K = N - N2;
	
	const int p = 9;
	float ep = 0.01;
	std::string path = "conf/N_65536_001.pc";
	volatile float H2 = -ep * log2(ep) - (1 - ep) * log2(1 - ep);
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 100;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int nf = 0;
	volatile int nt = 0;
	std::cout << std::endl << "start" << std::endl;
	while (count1 < t ) {
		a++;
		nf += N2;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		std::vector<int> reg(p, 0);
		uint32_t ccc = calculateCRC32(u, u.size());
		fb_generate->generate(fp);
		int times = 1;
		
		while (!decoder->decodescl(llr, u, fp, decbuff2, ccc, 16))
		{

			if (times == 4) {
				count1++;
				break;
			}
			nf += blength;
			fb_generate->generate(fp,N2+ times* blength);
			times++;
		}
		nt += times;
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< " errorpc_number: " << count1 << " fer1: " << (float)count1 / a
		<<" f: "<<(float)nf/a/N/H2
		<< " t: " << (float)nt / a
		<< std::endl;

}

void repolartestpc()
{
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 18000;
	int blength = 2000;
	int K = N - N2;

	const int p = 9;
	float ep = 0.04;
	std::string path = "conf/N_65536_002.pc";
	volatile float H2 = -ep * log2(ep) - (1 - ep) * log2(1 - ep);
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 100;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int nf = 0;
	volatile int nt = 0;
	std::cout << std::endl << "start" << std::endl;
	while (count1 < t) {
		a++;
		nf += N2;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		std::vector<int> reg(p, 0);
		uint32_t ccc = calculateCRC32(u, u.size());
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		fb_generate->generate(fp);
		int times = 1;

		while (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, 16))
		{

			if (times == 10) {
				count1++;
				break;
			}
			nf += blength;
			fb_generate->generate(fp, N2 + times * blength);
			times++;
		}
		nt += times;
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< " errorpc_number: " << count1 << " fer1: " << (float)count1 / a
		<< " f: " << (float)(nf / a + 32) / N / H2
		<< " t: " << (float)nt / a
		<< std::endl;

}

void segrepolartest(){
	//分段
	int h = 0xdeadbeef;
	int N = 65536;
	int n_seg = 4;
	int N2 = 10100;
	int blength = 500;
	int K = N - N2;

	const int p = 9;
	float ep = 0.02;
	std::string path = "conf/N_65536_002.pc";
	volatile float H2 = -ep * log2(ep) - (1 - ep) * log2(1 - ep);
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 100;
	volatile int a = 0;
	volatile int count1 = 0;
	std::vector<int> countseg (4,0);
	volatile int nf = 0;
	volatile int nt = 0;
	std::cout << std::endl << "start" << std::endl;
	while (count1 < t) {
		a++;
		int NF = N2;
		int NF_ = N2;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		std::vector<int> reg(p, 0);
		uint32_t ccc = calculateCRC32(u, u.size());
		std::vector<uint32_t> crc_v = calculateCRC32(u, N, n_seg);
		std::vector<int> seg_s (n_seg,0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		fb_generate->generate(fp);
		int times = 1;

		/*if (!decoder->decodepcsegment(llr, PC, fp, decbuff, p, ccc, crc_v, seg_s, 16)) {
			count1++;
		}
		for (int j = 0; j < n_seg; j++) {
			if (seg_s[j] == 1)
				countseg[j]= countseg[j]+1;
		}*/
		
		while (!decoder->decodepcsegment(llr, PC, fp, decbuff, p, ccc,crc_v, seg_s, 16))
		{
			if (times == 4) {
				count1++;
				break;
			}
			//fb_generate->generate(fp, N2 + times * blength);
			fb_generate->generate(fp,NF ,NF_, blength, seg_s, n_seg);
			times++;
		}
		nf += NF;
		nt += times;
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< " errorpc_number: " << count1 << " fer1: " << (float)count1 / a
		<< " f: " << (float)(nf / a + 32*4) / N / H2
		<< " t: " << (float)nt / a
		<< std::endl;

}

void test_nc1()
{
	//QBER与R不匹配
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 5000;
	int blength = 1000;
	int K = N - N2;

	const int p = 9;
	float ep = 0.01;
	std::string path = "conf/N_65536_001.pc";
	volatile float H2 = -ep * log2(ep) - (1 - ep) * log2(1 - ep);
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 100;
	volatile int a = 0;
	volatile int nf = 0;
	volatile int nt = 0;
	std::cout << std::endl << "start" << std::endl;
	while (a/t < t) {
		a++;
		nf += N2;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		std::vector<int> reg(p, 0);
		uint32_t ccc = calculateCRC32(u, u.size());
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		fb_generate->generate(fp);
		int times = 1;

		while (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, 16))
		{
			nf += blength;
			fb_generate->generate(fp, N2 + times * blength);
			times++;
		}
		nt += times;
	}
	std::cout << std::endl
		<< path << "  "
		<< "frame_number: " << a
		<< " f: " << (float)(nf / a + 32) / N / H2
		<< " t: " << (float)nt / a
		<< std::endl;

}

void polar_ml()
{
	//两个极化码合并
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 10500;
	int K = N - N2;
	const int p = 9;
	float ep = 0.02;
	std::string path = "conf/N_65536_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 100;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;

	int N2N = N * 2;
	int K2K = K * 2;
	std::unique_ptr<Decoder_polar<>>decoder2N = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K2K, N2N, fp, 1));
	std::vector<bool>	fp2(N2N);
	std::fill(fp2.begin(), fp2.end(), true);
	for (auto i = 0; i < N; i++)
		fp2[i] = fp[i];
	for (auto i = 0; i < N; i++)
		fp2[i+N] = fp[i];
	std::vector<int>	u2(N);
	std::vector<float>  llr2(N);
	std::vector<int>	decbuffN2N(N2N);
	while (count1 < t) {
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		if (!decoder->decodescl(llr, u, fp, decbuff2, ccc, 16)) {
			count1++;
			source->generate(u2);
			channel->bittollr(u2, llr2);
			std::vector<float>  llrN2N;
			llrN2N.insert(llrN2N.end(), llr2.begin(), llr2.end());
			llrN2N.insert(llrN2N.end(), llr.begin(), llr.end());
			std::vector<int>	u2N(N2N);
			encoder->encode(u2);
			for (int i = 0; i < N; i++) {
				u2N[i] = u[i] ^ u2[i];
			}
			for (int i = 0; i < N; i++) {
				u2N[i+N] =u[i];
			}
			uint32_t ccc2 = calculateCRC32(u2, u2.size());
			calculateCRC32(ccc2,u, u.size());
			
			if (!decoder2N->decodescl(llrN2N, u2N, fp2, decbuffN2N, ccc2, 16)) {
				count2++;
			}


		
		}
	}
	std::cout << std::endl << " number_fbit" << N2
		<< "frame_number: " << a
		<< "  errorpc_number: " << count1 << "  fer1:  " << (float)count1 / a
		<< std::endl;

}

void polartest2()
{
	//不同长度奇偶校验
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 220;
	int K = N - N2;
	const int p = 9;
	const int p2 = 17;	
	float ep = 0.02;
	int L = 16;
	std::string path = "conf/N_1024_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits
	std::vector<int>	PC3(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 500;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	while ( count1 < t||count2<t) {
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		//std::vector<int> reg(p, 0);
		//std::vector<int> last_statu(p, 0);
		//for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
		//	reg[0] = reg[0] ^ u[encodepc_loop]; // XOR operation
		//	PC[encodepc_loop] = reg[0];
		//	if ((encodepc_loop - last_statu[0]) >= p * 3) {
		//		for(int i= encodepc_loop-2*p;i> last_statu[0];i-=2*p)
		//		PC[encodepc_loop] = PC[encodepc_loop] ^ u[i];
		//	}
		//	if (fp[encodepc_loop]) { last_statu[0] = encodepc_loop; }
		//	last_statu = CSR(last_statu);
		//	reg = CSR(reg); // Perform cyclic shift right
		//} 
		std::vector<int> reg(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		std::vector<int> reg2(p2, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC2[encodepc_loop] = reg2[0] ^ u[encodepc_loop]; // XOR operation
			reg2[0] = PC2[encodepc_loop];
			reg2 = CSR(reg2); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepc4(llr, PC,PC2, fp, decbuff2, p,p2, ccc, L)) {
			count1++;
		}
		//if (!decoder->decodepc5(llr, PC, fp, decbuff2, p, ccc, L)) {
		//	count1++;
		//}
		if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, L)) {
			count2++;
		}
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorpc0: " << count2 << "  fer0:  "
		<< (float)count2 / a
		<< std::endl;

}

void polartest3()
{
	//不同长度奇偶校验
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 200;
	int K = N - N2;
	const int p = 9;
	const int p2 = 18;
	float ep = 0.02;
	int L = 16;
	std::string path = "conf/N_1024_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));
		std::vector<int> numbers;
		numbers.push_back(56);
		numbers.push_back(64);
		numbers.push_back(68);
		numbers.push_back(70);
		numbers.push_back(72);
		numbers.push_back(96);
		numbers.push_back(128);
		numbers.push_back(129);
		numbers.push_back(130);
		numbers.push_back(131);
		numbers.push_back(132);
		numbers.push_back(133);
		numbers.push_back(134);
		numbers.push_back(136);
		numbers.push_back(160);
		numbers.push_back(162);
		numbers.push_back(192);
		numbers.push_back(256);
		numbers.push_back(257);
		numbers.push_back(258);
		numbers.push_back(259);
		numbers.push_back(260);
		numbers.push_back(261);
		numbers.push_back(262);
		numbers.push_back(264);
		numbers.push_back(288);
		numbers.push_back(320);
		numbers.push_back(512);
		numbers.push_back(513);
		numbers.push_back(514);
		numbers.push_back(515);
		numbers.push_back(516);
		numbers.push_back(517);
		numbers.push_back(518);
		numbers.push_back(520);
		numbers.push_back(544);
		numbers.push_back(576);
	int t = 500;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	while (count1 < t || count2 < t) {

		std::vector<bool>	fp2(N, false);
		for (int num : numbers) {
			fp2[num] =true;
		}
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		std::vector<int> reg(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G

		if (!decoder->decodepc6(llr, PC, u, fp, fp2, decbuff2, p, ccc, L)) {
			count1++;
		}
		if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, L)) {
			count2++;
		}
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorpc0: " << count2 << "  fer0:  "
		<< (float)count2 / a
		<< std::endl;

}

void polartest4()
{
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 200;
	int K = N - N2;
	const int p = 9;
	int p2 = 9;
	float ep = 0.02;
	int L = 16;
	std::string path = "conf/N_1024_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits
	std::vector<int>	PC3(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));
	std::string path2 = "conf/P9_1024_002.pc";
	std::vector<double> z;
	std::ifstream inputFile(path2);
	if (!inputFile.is_open()) {
		std::cout << "Failed to open channel parameter file: " << path2 << std::endl;
		return;
	}
	double value;
	while (inputFile >> value) {
		z.push_back(value);
	}
	inputFile.close();

	int t = 500;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	while (count1 < t || count2 < t) {
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());
		std::vector<int> reg(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, L)) {
			count2++;
		}

		std::vector<int> reg2(p2, 0);
		std::vector<double> z_reg2(p2, 0);
		std::vector<double> fb_reg2;
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
				if (fp[encodepc_loop]) {
					if (z_reg2[0] + z[encodepc_loop] <= 1) {
						PC2[encodepc_loop] = reg2[0] ^ u[encodepc_loop];
						reg2[0] = PC2[encodepc_loop];
						z_reg2[0] = z_reg2[0] + z[encodepc_loop];
						fb_reg2.push_back(z_reg2[0]);
						z_reg2[0] = 0;
					}
					else {
						PC2[encodepc_loop] = u[encodepc_loop];
						fb_reg2.push_back(z[encodepc_loop]);
					}
				}
				else {
					PC2[encodepc_loop] = reg2[0] ^ u[encodepc_loop];
					reg2[0] = PC2[encodepc_loop];
					z_reg2[0] = z_reg2[0] + z[encodepc_loop];
				}
				reg2 = CSR(reg2);
				z_reg2 = CSR(z_reg2);
		}  
		if (!decoder->decodepc7(llr, PC2,z, fp, decbuff2, p2, ccc, L)) {
			count1++;
		}
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorpc0: " << count2 << "  fer0:  "
		<< (float)count2 / a
		<< std::endl;

}

void polartest5()
{
	//200-220tehua
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 202;
	int K = N - N2;
	const int p = 9;
	const int p2 = 9;
	float ep = 0.02;
	int L = 16;
	std::string path = "conf/N_1024_002.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<bool>	fp2(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));
		std::vector<int> column1 = { 48 - 1,93 - 1, 103 - 1, 106 - 1, 154 - 1, 155 - 1, 166 - 1, 225 - 1, 279 - 1, 282 - 1, 292 - 1, 337 - 1, 393 - 1, 535 - 1, 538 - 1, 548 - 1, 585 - 1, 593 - 1, 645 - 1, 770 - 1 };
	

	int t = 500;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	fb_generate->generate(fp2, 200);
	while (count1 < t || count2 < t) {

		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		std::vector<int> reg(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepc(llr, PC, fp, decbuff2, p, ccc, L)) {
			count2++;
		}
		std::vector<int>	PCE(2, 0); // Resulting encoded bits

		for (int encodepc_loop = 0; encodepc_loop < 10; ++encodepc_loop) {
			PCE[0] = PCE[0] ^ u[column1[encodepc_loop]]; // XOR operation
		}  
		for (int encodepc_loop = 0; encodepc_loop < 10; ++encodepc_loop) {
			PCE[1] = PCE[1] ^ u[column1[encodepc_loop+10]]; // XOR operation
		}
		if (!decoder->decodepc8(llr, PC, PCE, column1,fp2, decbuff2, p, ccc, L)) {
			count1++;
		}
	}
	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorpc0: " << count2 << "  fer0:  "
		<< (float)count2 / a
		<< std::endl;

}

void polartestnew1()
{
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 113;
	int K = N - N2;
	const int p = 9;
	const int p2 = 9;
	float ep = 0.01;
	int L = 16;
	std::string path = "conf/N_1024_001.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 1000;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	volatile int count3 = 0;
	std::vector<std::vector<int>> error_site1;
	std::vector<std::vector<int>> error_site2;
	std::vector<std::vector<int>> error_site3;
	while ( count1 < t) {
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());

		//if (!decoder->decodescltest(llr, u, fp, decbuff2, ccc, L, error_site2)) {
		//	count2++;
		//}

		//std::vector<int> reg(p, 0);
		//for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
		//	PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
		//	reg[0] = PC[encodepc_loop];
		//	reg = CSR(reg); // Perform cyclic shift right
		//}  //huawei 5G
		//if (!decoder->decodepctest(llr, PC, u, fp, decbuff2,p, ccc, L, error_site1)) {
		//	count1++;
		//}


		std::vector<int> reg2(p2, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			if (fp[encodepc_loop]) {
				PC2[encodepc_loop] = reg2[0] ^ u[encodepc_loop]; // XOR operation
			}
			else {
				reg2[0] = reg2[0] ^ u[encodepc_loop]; // XOR operation
				PC2[encodepc_loop] = reg2[0];
			}
			reg2 = CSR(reg2); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepctest2(llr, PC2, u, fp, decbuff2, p2, ccc, L, error_site3)) {
			count1++;
		}

	}

	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorpc: " << count1 << "  fer1:  "
		<< (float)count1 / a
		<< "  errorscl: " << count2 << "  fer0:  "
		<< (float)count2 / a
		//<< "  errorpc2: " << count3 << "  fer2:  "
		//<< (float)count3 / a
		<< std::endl;

}

void polartestnew2()
{
	int h = 0xdeadbeef;
	int N = 1024;
	int N2 = 100;
	int K = N - N2;
	const int p = 9;
	const int p2 = 9;
	float ep = 0.01;
	int L = 16;
	std::string path = "conf/N_1024_001.pc";
	std::vector<bool>	fp(N);//Freeze bit position
	std::vector<int>	u(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	decbuff2(N);
	std::vector<int>	PC(N, 0); // Resulting encoded bits
	std::vector<int>	PC2(N, 0); // Resulting encoded bits
	std::vector<int>	PC3(N, 0); // Resulting encoded bits

	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 1000;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;
	volatile int count3 = 0;
	volatile int count4 = 0;
	std::vector<std::vector<int>> error_site1;
	std::vector<std::vector<int>> error_site2;
	std::vector<std::vector<int>> error_site3;
	std::vector<std::vector<int>> error_site4;
	while (count1 < t || count2 < t) {
		a++;
		source->generate(u);
		channel->add_noise(u, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());
		if (!decoder->decodescltest(llr, u, fp, decbuff2, ccc, L, error_site1)) {
			count1++;
		}

		std::vector<int> reg(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			PC[encodepc_loop] = reg[0] ^ u[encodepc_loop]; // XOR operation
			reg[0] = PC[encodepc_loop];
			reg = CSR(reg); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepctest(llr, PC, u, fp, decbuff2, p, ccc, L, error_site2)) {
			count2++;
		}

		std::vector<int> reg2(p2, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			if (fp[encodepc_loop]) {
				PC2[encodepc_loop] = reg2[0] ^ u[encodepc_loop]; // XOR operation
			}
			else {
				reg2[0] = reg2[0] ^ u[encodepc_loop]; // XOR operation
				PC2[encodepc_loop] = reg2[0];
			}
			reg2 = CSR(reg2); // Perform cyclic shift right
		}  //huawei 5G
		if (!decoder->decodepctest2(llr, PC2, u, fp, decbuff2, p2, ccc, L, error_site3)) {
			count3++;
		}

		std::vector<int> regpc(p*2, 0);
		std::vector<int> regpclast(p*2, 0);
		std::vector<int> _pc(p, 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			if (fp[encodepc_loop]) {
				if (_pc[encodepc_loop%p] == 0) {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					// 
				}
				else {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)+p] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					//

				}
				regpclast[encodepc_loop % p] = regpc[0];
				regpclast[encodepc_loop % p + p] = regpc[p];
				_pc[encodepc_loop % p] = 1 - _pc[encodepc_loop % p];
			}
			else {
				regpc[0] = regpc[0]^u[encodepc_loop]; // XOR operation
				PC3[encodepc_loop] = regpc[0];
			}
			regpc = CSR(regpc); 
		}
		if (!decoder->decodepctest3(llr, PC3, u, fp, decbuff2, p, ccc, L, error_site4)) {
			count4++;
		}

	}
	/*saveToTxt(error_site3, "errorpc2test.txt");
	saveToTxt(error_site2, "errorpctest2.txt");
	saveToTxt(error_site1, "errorscltest2.txt");*/
	std::cout << std::endl
		<< "frame_number: " << a
		<< "  errorscl: " << count1 << "  fer0:  "
		<< (float)count1 / a
		<< "  errorpc: " << count2 << "  fer1:  "
		<< (float)count2 / a
		<< "  errorpc2: " << count3 << "  fer2:  "
		<< (float)count3 / a
		<< "  errorpc3: " << count4 << "  fer3:  "
		<< (float)count4 / a
		<< std::endl;

}


const size_t N = 65536;  // 每次读取的字符个数
void readBits(std::ifstream& file, std::vector<int>& buffer, size_t& count) {
	buffer.clear();
	char ch;
	count = 0;
	for (size_t i = 0; i < N && file.get(ch); ++i) {
		if (ch == '0' || ch == '1') {
			buffer.push_back(ch - '0');
			++count;
		}
	}
}

void polartestnew3()
{
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 11600;
	int K = N - N2;// 0.005-5  0.01 -8   0.015 -8  0.02 -9 0.025 -9
	const int p = 9;
	volatile float ep = 0.03;
	int L = 16;
	std::string path = "conf/N_1024_003.pc";
	std::vector<bool>	fp(N);//Freeze bit position


	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	int t = 1000;
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;

	std::ifstream infile1("ksa.txt");
	std::ifstream infile2("ksb.txt");
	if (!infile1 || !infile2) {
		std::cerr << "无法打开文件!" << std::endl;
	}
	std::vector<std::vector<int>> error_site1;
	std::vector<std::vector<int>> error_site2;
	std::vector<std::vector<int>> nftxts;
	std::vector<int>	u(N);
	std::vector<int>	u2(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	PC3(N, 0); // Resulting encoded bits
	std::vector<int> regpc(p * 2, 0);
	std::vector<int> regpclast(p * 2, 0);
	std::vector<int> _pc(p, 0);
	size_t idx1 = 0, idx2 = 0;
	int blength = 10;
	while (infile1 && infile2) {
		std::vector<int>	nftxt;
		readBits(infile1, u, idx1);
		readBits(infile2, u2, idx2);
		size_t differences = 0;
		for (size_t i = 0; i < N; ++i) {
			if (u[i] != u2[i]) {
				++differences;
				
			}
		}
		//std::cout << "differences:" << differences << std::endl;
		if (idx1 < N || idx2 < N ) {
			break;
		}
		if ((differences > N * ep)  ||   (differences< N * (ep-0.005))) { continue; }
		nftxt.push_back(differences);

		a++;
		channel->bittollr(u2, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());




		fb_generate->generate(fp);
		int N3 = N2;
		while(!decoder->decodescltest(llr, u, fp, decbuff, ccc, L, error_site1)) {
			N3 += blength;
			fb_generate->generate(fp, N3);
		}
		nftxt.push_back(N3);
		count1 += N3;



		fb_generate->generate(fp);
		N3 = N2;
		std::fill(regpc.begin(), regpc.end(), 0);
		std::fill(regpclast.begin(), regpclast.end(), 0);
		std::fill(_pc.begin(), _pc.end(), 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			if (fp[encodepc_loop]) {
				if (_pc[encodepc_loop % p] == 0) {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					// 
				}
				else {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)+p] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					//

				}
				regpclast[encodepc_loop % p] = regpc[0];
				regpclast[encodepc_loop % p + p] = regpc[p];
				_pc[encodepc_loop % p] = 1 - _pc[encodepc_loop % p];
			}
			else {
				regpc[0] = regpc[0] ^ u[encodepc_loop]; // XOR operation
				PC3[encodepc_loop] = regpc[0];
			}
			regpc = CSR(regpc);
		}
		while (!decoder->decodepctest3(llr, PC3, u, fp, decbuff, p, ccc, L, error_site2)) {
			N3 += blength;
			fb_generate->generate(fp, N3);
			std::fill(regpc.begin(), regpc.end(), 0);
			std::fill(regpclast.begin(), regpclast.end(), 0);
			std::fill(_pc.begin(), _pc.end(), 0);
			for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
				if (fp[encodepc_loop]) {
					if (_pc[encodepc_loop % p] == 0) {
						PC3[encodepc_loop] = regpclast[encodepc_loop % (p)] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
						// 
					}
					else {
						PC3[encodepc_loop] = regpclast[encodepc_loop % (p)+p] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
						//

					}
					regpclast[encodepc_loop % p] = regpc[0];
					regpclast[encodepc_loop % p + p] = regpc[p];
					_pc[encodepc_loop % p] = 1 - _pc[encodepc_loop % p];
				}
				else {
					regpc[0] = regpc[0] ^ u[encodepc_loop]; // XOR operation
					PC3[encodepc_loop] = regpc[0];
				}
				regpc = CSR(regpc);
			}
		}
		nftxt.push_back(N3);
		count2 += N3;
		nftxts.push_back(nftxt);
	}
	infile1.close();
	infile2.close();
	saveToTxt(nftxts, "nf.txt");
	std::cout << std::endl
		<< "times " << a
		<< "  nf_scl:  "<< (float)count1 / a
		<< "  nf_pc:  "<< (float)count2 / a
		<< std::endl;
}

void polartestnew4()
{
	int h = 0xdeadbeef;
	int N = 65536;
	int N2 = 7700;
	int K = N - N2;// 0.005-5  0.01 -8   0.015 -8  0.02 -9 0.025 -9
	const int p = 9;
	volatile float ep = 0.01;
	int L = 16;
	std::string path = "conf/N_65536_001.pc";
	std::vector<bool>	fp(N);//Freeze bit position


	std::unique_ptr<Source_random<>>source = std::unique_ptr<Source_random<>>(new Source_random<>(N));
	std::unique_ptr<Frozenbits_generator>fb_generate = std::unique_ptr<Frozenbits_generator>(new Frozenbits_generator(K, N, path));
	fb_generate->generate(fp);
	std::unique_ptr<Encoder_polar<>>encoder = std::unique_ptr<Encoder_polar<>>(new Encoder_polar<>(K, N, fp, 1));
	std::unique_ptr<Channel_bsc<>>channel = std::unique_ptr<Channel_bsc<>>(new Channel_bsc<>(ep, 0));
	std::unique_ptr<Decoder_polar<>>decoder = std::unique_ptr<Decoder_polar<>>(new Decoder_polar<>(K, N, fp, 1));

	
	volatile int a = 0;
	volatile int count1 = 0;
	volatile int count2 = 0;

	std::ifstream infile1("ksa.txt");
	std::ifstream infile2("ksb.txt");
	if (!infile1 || !infile2) {
		std::cerr << "无法打开文件!" << std::endl;
	}
	std::vector<std::vector<int>> error_site1;
	std::vector<std::vector<int>> error_site2;
	std::vector<std::vector<int>> nftxts;
	std::vector<int>	u(N);
	std::vector<int>	u2(N);
	std::vector<float>  llr(N);
	std::vector<int>	decbuff(N);
	std::vector<int>	PC3(N, 0); // Resulting encoded bits
	std::vector<int> regpc(p * 2, 0);
	std::vector<int> regpclast(p * 2, 0);
	std::vector<int> _pc(p, 0);
	size_t idx1 = 0, idx2 = 0;
	size_t differences = 0;
	int times = 0;
	while (infile1 && infile2) {
		
		std::vector<int>	nftxt;
		readBits(infile1, u, idx1);
		readBits(infile2, u2, idx2);
		size_t difference = 0;
		for (size_t i = 0; i < N; ++i) {
			if (u[i] != u2[i]) {
				++difference;
			}
		}
		std::cout << times++ << "   ber: " << (float)difference/N << std::endl;
		if (idx1 < N || idx2 < N) {
			break;
		}


		if (difference > N * (ep+0.003)) { continue; }
		differences += difference;
		a++;
		channel->bittollr(u2, llr);
		encoder->encode(u);
		uint32_t ccc = calculateCRC32(u, u.size());


		if (!decoder->decodescltest(llr, u, fp, decbuff, ccc, L, error_site1)) {
			count1++;
		}

		std::fill(regpc.begin(), regpc.end(), 0);
		std::fill(regpclast.begin(), regpclast.end(), 0);
		std::fill(_pc.begin(), _pc.end(), 0);
		for (int encodepc_loop = 0; encodepc_loop < N; ++encodepc_loop) {
			if (fp[encodepc_loop]) {
				if (_pc[encodepc_loop % p] == 0) {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					// 
				}
				else {
					PC3[encodepc_loop] = regpclast[encodepc_loop % (p)+p] ^ regpc[0] ^ regpc[p] ^ u[encodepc_loop];
					//

				}
				regpclast[encodepc_loop % p] = regpc[0];
				regpclast[encodepc_loop % p + p] = regpc[p];
				_pc[encodepc_loop % p] = 1 - _pc[encodepc_loop % p];
			}
			else {
				regpc[0] = regpc[0] ^ u[encodepc_loop]; // XOR operation
				PC3[encodepc_loop] = regpc[0];
			}
			regpc = CSR(regpc);
		}
		if (!decoder->decodepctest3(llr, PC3, u, fp, decbuff, p, ccc, L, error_site2)) {
			count2++;
		}
	}
	infile1.close();
	infile2.close();
	std::cout << std::endl
		<< "times " << a
		<< "  fer_scl:  " << (float)count1 / a
		<< "  fer_pc:  " << (float)count2 / a
		<< "  BER:  " << (float)differences / a / N
		<< std::endl;
}

#endif /*TESTPOLAR_HPP_*/