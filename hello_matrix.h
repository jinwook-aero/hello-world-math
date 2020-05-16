// hello_matrix.h
// Header file that contains classes of hello_matrix
//
// Authors      : Jinwook Lee, Sungwook Lee
// First version: May 16, 2020
// Last update  : May 16, 2020
//

#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

// Matrix class contains data as 2D vector of T
// and provides matrix operations
template <typename T = double>
class Matrix
{
public:
	Matrix(); // To indicate invalid matrix
	Matrix(size_t nRow, size_t nCol, T initialVal = 0, bool isDiag = false, bool isRand = false); // Random numbers in the range of [0-initialVal]
	Matrix(const Matrix&); // Copy constructor
	~Matrix() = default;
	
	// Element access
	T& operator()(size_t iRow, size_t iCol);

	// Diplay the elements of matrix in provided range
	void Display(); // All range
	void Display(size_t nRowMin, size_t nRowMax); // Row range
	void Display(size_t nRowMin, size_t nRowMax, size_t nColMin, size_t nColMax); // Row-Col range

	// Element-wise operation
	bool operator==(const Matrix&);  // If all elements are identical
	Matrix operator+(const Matrix&); // matrix element-wise addition
	Matrix operator-(const Matrix&); // matrix element-wise subtraction
	Matrix ElemMultiply(const Matrix&); // matrix element-wise multiplication
	Matrix ElemDivide(const Matrix&); // matrix element-wise multiplication
	
	Matrix operator/(const Matrix&) = delete; // It can be a confusing operator, thus removed

	// Matrix operation
	Matrix operator*(const Matrix&); // Matrix multiplication
	Matrix BackSlash(const Matrix&); // A\b (== inv(A)*b) is implemented as A.backSlash(b)
	Matrix Transpose(); // Transpose of matrix
	Matrix Inverse();   // Inverse of matrix

private:
	size_t _nRow, _nCol;
	std::vector<std::vector<T>> _elemData;
	bool _validity;
	void Resize(size_t nRow, size_t n_Col);
	bool IsEqualSize(const Matrix&, const Matrix&);
	bool IsValid();
};

// Constructors
template <typename T>
Matrix<T>::Matrix():
_nRow(0), _nCol(0), _validity(false)
{
	_elemData.clear();
}

template <typename T>
Matrix<T>::Matrix(size_t nRow, size_t nCol, T initialVal, bool isDiag, bool isRand):
_nRow(nRow), _nCol(nCol), _validity(true)
{
	_elemData.clear();
	_elemData.resize(_nRow, std::vector<T>(_nCol, initialVal));
	if (isRand){
		// Random number generation sequence by : https://stackoverflow.com/questions/9878965/rand-between-0-and-1
		std::mt19937_64 rng;
		uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
		rng.seed(ss);
		std::uniform_real_distribution<double> unif(0, 1);

		for (size_t iRow = 0; iRow < _nRow; ++iRow) {
		for (size_t iCol = 0; iCol < _nCol; ++iCol) {
			_elemData[iRow][iCol] = unif(rng) * initialVal; // Random numbers in the range of [0-initialVal]
		}
		}
	}

	if (isDiag && (initialVal != 0)){
		for (size_t iRow=0; iRow<_nRow; ++iRow){
		for (size_t iCol=0; iCol<_nCol; ++iCol){
			if (iRow != iCol)
				_elemData[iRow][iCol] = 0;			
		}
		}
	}	
}

template <typename T>
Matrix<T>::Matrix(const Matrix& rhs) // Copy constructor
{
	_validity = rhs._validity;
	this->Resize(rhs._nRow, rhs._nCol);
	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		_elemData[iRow][iCol] = rhs._elemData[iRow][iCol];
	}
	}
}

// Element access
template <typename T>
T & Matrix<T>::operator()(size_t iRow, size_t iCol)
{
	return _elemData[iRow][iCol];
}

// Diplay the elements of matrix in provided range
template <typename T>
void Matrix<T>::Display(size_t nRowMin, size_t nRowMax,
						size_t nColMin, size_t nColMax)
{
	// Proper display range
	size_t nZero{ 0 };
	size_t nRow0 = std::max(nRowMin, nZero); size_t nRow1 = std::min(nRowMax, _nRow);
	size_t nCol0 = std::max(nColMin, nZero); size_t nCol1 = std::min(nColMax, _nCol);

	for (size_t iRow = nRow0; iRow < nRow1; ++iRow) {
	for (size_t iCol = nCol0; iCol < nCol1; ++iCol) {
		std::cout << _elemData[iRow][iCol] << "\t";
	}
	std::cout << "\n";
	}
	std::cout << "\n";
}

template <typename T>
void Matrix<T>::Display(size_t nRowMin, size_t nRowMax)
{
	Display(nRowMin, nRowMax, 0, _nCol);
}

template <typename T>
void Matrix<T>::Display()
{
	Display(0, _nRow, 0, _nCol);
}

// Element-wise operators
template <typename T>
bool Matrix<T>::operator==(const Matrix& rhs) // If all elements are identical
{
	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		if (this->_elemData[iRow][iCol] != rhs._elemData[iRow][iCol])
			return false;
	}
	}
	return true;
}


template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& rhs) // Element-wise addition
{
	Matrix<T> tempMat{ _nRow,_nCol };
	if (!IsEqualSize(*this, rhs)){
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		tempMat._elemData[iRow][iCol] = this->_elemData[iRow][iCol] + rhs._elemData[iRow][iCol];
	}
	}
	return tempMat;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix& rhs) // Element-wise subtraction
{
	Matrix<T> tempMat{ _nRow,_nCol };
	if (!IsEqualSize(*this, rhs)) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		tempMat._elemData[iRow][iCol] = this->_elemData[iRow][iCol] - rhs._elemData[iRow][iCol];
	}
	}
	return tempMat;
}


template <typename T>
Matrix<T> Matrix<T>::ElemMultiply(const Matrix& rhs) // Element-wise multipllication
{
	Matrix<T> tempMat{ _nRow,_nCol };
	if (!IsEqualSize(*this, rhs)) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		tempMat._elemData[iRow][iCol] = this->_elemData[iRow][iCol] * rhs._elemData[iRow][iCol];
	}
	}
	return tempMat;
}


template <typename T>
Matrix<T> Matrix<T>::ElemDivide(const Matrix& rhs) // Element-wise division
{
	Matrix<T> tempMat{ _nRow,_nCol };
	if (!IsEqualSize(*this, rhs)) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < _nRow; ++iRow) {
	for (size_t iCol = 0; iCol < _nCol; ++iCol) {
		tempMat._elemData[iRow][iCol] = this->_elemData[iRow][iCol] / rhs._elemData[iRow][iCol];
	}
	}
	return tempMat;
}

// Matrix arithmetic
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& rhs) // Matrix multiplication
{
	size_t nRow0 = this->_nRow;
	size_t nSum0 = this->_nCol;
	size_t nCol0 = rhs._nCol;

	Matrix<T> tempMat{ nRow0,nCol0 };
	if (this->_nCol != rhs._nRow) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
	for (size_t iCol = 0; iCol < nCol0; ++iCol) {
	for (size_t iSum = 0; iSum < nSum0; ++iSum) {
		tempMat._elemData[iRow][iCol] += this->_elemData[iRow][iSum] * rhs._elemData[iSum][iCol];
	}
	}
	}
	return tempMat;
}

template <typename T>
Matrix<T> Matrix<T>::Transpose() // Matrix transpose
{
	size_t nRow0 = this->_nCol;
	size_t nCol0 = this->_nRow;
	Matrix<T> tempMat{ nRow0,nCol0 };

	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
	for (size_t iCol = 0; iCol < nCol0; ++iCol) {
		tempMat._elemData[iRow][iCol] = this->_elemData[iCol][iRow];
	}
	}
	return tempMat;
}

// Private functions
template <typename T>
bool Matrix<T>::IsEqualSize(const Matrix& lhs, const Matrix& rhs)
{
	if (lhs._nRow == rhs._nRow)
		if (lhs._nCol == rhs._nCol)
			return true;
	return false;
}

template <typename T>
void Matrix<T>::Resize(size_t nRow, size_t nCol)
{
	_nRow = nRow;
	_nCol = nCol;
	_elemData.clear();
	_elemData.resize(_nRow, std::vector<T>(_nCol, 0));
}

template <typename T>
bool Matrix<T>::IsValid()
{
	return _validity;
}