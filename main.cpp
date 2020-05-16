// main.cpp 
// Main file for testing hello_math
//
// Authors      : Jinwook Lee, Sungwook Lee
// First version: May 16, 2020
// Last update  : May 16, 2020
//

#pragma once
#include <iostream>
#include "hello_matrix.h"

int main()
{
	Matrix<double> invalidMatrix{};
	Matrix<> I{3,3,1,true}; // identity matrix
	Matrix<> R1{ 3,3,1,false,true }; // non-diagonal  random number matrix
	Matrix<> R2{ 3,3,1,false,true }; // non-diagonal  random number matrix
	Matrix<> D = R1 *R2;

	R1.Display(0, 5, 0, 5);
	R2.Display(0, 5, 0, 5);
	D.Display();

	// End of program
	return 0;
}