#include <iostream>


class Statistic
{
public:
	int frame_number;
	int error_number1;
	int error_number2;

public:
	Statistic();

	void update(int frame_number, int error_number1, int error_number2);

	void display();
};

Statistic
::Statistic()
	:frame_number(0), error_number1(0), error_number2(0)
{
}

inline void Statistic
::update(int frame_number, int error_number1, int error_number2){
	this->frame_number = frame_number;
	this->error_number1 = error_number1;
	this->error_number2 = error_number2;
}

inline void Statistic
::display()
{
	std::cout
		<< "frame_number: " << frame_number
		<< "  error1: "		<< error_number1 
		<< "  fer1:  "		<< (float)error_number1 / frame_number
		<< "  error2: "		<< error_number2 
		<< "  fer0:  "		<< (float)error_number2 / frame_number
		<< std::endl;
}
