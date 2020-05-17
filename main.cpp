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
	const size_t matSize = 1000;
	Matrix<double> invalid{};
	Matrix<> I{ matSize,matSize,1,true }; // identity matrix
	Matrix<> R1{ matSize,matSize,1,false,true }; // non-diagonal random number matrix
	Matrix<> R2{ matSize,matSize,1,false,true }; // non-diagonal random number matrix

	auto start = std::chrono::high_resolution_clock::now();
	Matrix<> D = I*R1;
	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
	std::cout << duration.count() << std::endl;
	
	/*
	invalid.Display();
	R1.Display(0, 5, 0, 5);
	R2.Display(0, 5, 0, 5);
	D.Display(0,5,0,5);
	*/

	// End of program
	return 0;
}