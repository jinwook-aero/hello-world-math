// main.cpp 
// Main file for testing hello_math
//
// Authors      : Jinwook Lee, Sungwook Lee
// First version: May 16, 2020
// Last update  : May 16, 2020
//

#pragma once
#include <iostream>
#include <chrono> 
#include "hello_matrix.h"

int main()
{
	const size_t matSize = 1024;
	Matrix<> R1{ matSize,matSize,1,false,true }; // non-diagonal random number matrix
	Matrix<> R2{ matSize,matSize,1,false,true }; // non-diagonal random number matrix

	Matrix <> SlowM, FastM;
	for (size_t round = 0 ; round <2 ; ++round){
		auto start = std::chrono::high_resolution_clock::now();
		if (round == 0)
			SlowM = R1.SlowMultiply(R2);
		else
			FastM = R1*R2;
		auto stop = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
		std::cout << "Round " << round << "\t: " << duration.count() << " ms" << std::endl;
	}
	if (SlowM == FastM)
		std::cout << "Results are equal" << std::endl;
	else
		std::cout << "Results are not equal" << std::endl;

	// End of program
	return 0;
}