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

// Matrix class contains data as 1D vector of T
// and provides matrix operations
template <typename T = double>
class Matrix
{
public:
	Matrix(); // To indicate invalid matrix
	Matrix(size_t nRow, size_t nCol, T initialVal = 0, bool isDiag = false, bool isRand = false); // Random numbers in the range of [0-initialVal]
	Matrix(const Matrix&); // Copy constructor
	~Matrix();
	
	// Element access or submatrix
	T& operator()(size_t iRow, size_t iCol);
	Matrix operator()(size_t iRowMin, size_t iRowMax, size_t iColMin, size_t iColMax); // Row-Col range

	// Diplay the elements of matrix in provided range
	void Display(); // All range
	void Display(size_t iRowMin, size_t iRowMax); // Row range
	void Display(size_t iRowMin, size_t iRowMax, size_t iColMin, size_t iColMax); // Row-Col range

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
	// Private member variables
	size_t _nRow, _nCol;
	bool _validity;
	std::vector<T> _elemData;  // 1D vector is used to avoid performance overhead from 2D vector
	//T* _elemData;

	// Private member functions
	void Resize(size_t nRow, size_t n_Col);
	bool IsEqualSize(const Matrix&, const Matrix&);
	bool IsValid();
	Matrix RowMerge(const Matrix&, const Matrix&); // Merge as [LHS; RHS]
	Matrix ColMerge(const Matrix&, const Matrix&); // Merge as [LHS, RHS]
};

// Destructor first
template <typename T>
Matrix<T>::~Matrix()
{
	//delete[] _elemData; _elemData = nullptr;
}

// Constructors
template <typename T>
Matrix<T>::Matrix():
_nRow(0), _nCol(0), _validity(false)
{
	_elemData.clear();
	_elemData.resize(0,0);
	//_elemData = new T[0];
}

template <typename T>
Matrix<T>::Matrix(size_t nRow, size_t nCol, T initialVal, bool isDiag, bool isRand):
_nRow(nRow), _nCol(nCol), _validity(true)
{
	_elemData.clear();
	_elemData.resize(_nRow*_nCol, initialVal);
	//_elemData = new T[_nRow * _nCol];
	if (isRand){
		// Random number generation sequence by : https://stackoverflow.com/questions/9878965/rand-between-0-and-1
		std::mt19937_64 rng;
		uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
		rng.seed(ss);
		std::uniform_real_distribution<double> unif(0, 1);

		for (size_t iRun = 0; iRun < _nRow * _nCol; ++iRun)
			_elemData[iRun] = unif(rng) * initialVal; // Random numbers in the range of [0-initialVal]
	}

	if (isDiag && (initialVal != 0)){
		for (size_t iRow=0; iRow<_nRow; ++iRow){
		for (size_t iCol=0; iCol<_nCol; ++iCol){
			if (iRow != iCol)
				_elemData[iRow * _nCol + iCol] = 0;
		}
		}
	}	
}

template <typename T>
Matrix<T>::Matrix(const Matrix& rhs) // Copy constructor
{
	_validity = rhs._validity;
	this->Resize(rhs._nRow, rhs._nCol);
	for (size_t iRun = 0; iRun < _nRow*_nCol; ++iRun)
		_elemData[iRun] = rhs._elemData[iRun];
}

// Element access or submatrix
template <typename T>
T & Matrix<T>::operator()(size_t iRow, size_t iCol)
{
	return _elemData[iRow * _nCol + iCol];
}

template <typename T>
Matrix<T> Matrix<T>::operator()(size_t iRowMin, size_t iRowMax, size_t iColMin, size_t iColMax)
{
	const size_t nRowLhs0 = iRowMax - iRowMin + 1;
	const size_t nColLhs0 = iColMax - iColMin + 1;
	Matrix<T> tempMat{ nRowLhs0,nColLhs0 };

	for (size_t iRow = iRowMin; iRow <= iRowMax; ++iRow) {
	for (size_t iCol = iColMin; iCol <= iColMax; ++iCol) {
		const size_t iRowLhs  = iRow - iRowMin;
		const size_t iColLhs  = iCol - iColMin;
		tempMat._elemData[iRowLhs*nColLhs0 + iColLhs] = this->_elemData[iRow * _nCol + iCol];
	}
	}
	return tempMat;
}

// Diplay the elements of matrix in provided range
template <typename T>
void Matrix<T>::Display(size_t iRowMin, size_t iRowMax,
						size_t iColMin, size_t iColMax)
{
	// Validity check
	if (!IsValid()){
		std::cout << "NaN\n\n";
		return;
	}

	// Proper display range
	size_t nZero{ 0 };
	size_t nRow0 = std::max(iRowMin, nZero); size_t nRow1 = std::min(iRowMax, _nRow);
	size_t nCol0 = std::max(iColMin, nZero); size_t nCol1 = std::min(iColMax, _nCol);

	for (size_t iRow = nRow0; iRow < nRow1; ++iRow) {
	for (size_t iCol = nCol0; iCol < nCol1; ++iCol) {
		std::cout << _elemData[iRow * _nCol + iCol] << "\t";
	}
	std::cout << "\n";
	}
	std::cout << "\n";
}

template <typename T>
void Matrix<T>::Display(size_t iRowMin, size_t iRowMax)
{
	Display(iRowMin, iRowMax, 0, _nCol);
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
	for (size_t iRun = 0; iRun < _nRow* _nCol; ++iRun)
		if (this->_elemData[iRun] != rhs._elemData[iRun])
			return false;
	return true;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix& rhs) // Element-wise addition
{
	Matrix<T> tempMat{ _nRow,_nCol };
	if (!IsEqualSize(*this, rhs)) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRun = 0; iRun < _nRow * _nCol; ++iRun)
		tempMat._elemData[iRun] = this->_elemData[iRun] + rhs._elemData[iRun];
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

	for (size_t iRun = 0; iRun < _nRow * _nCol; ++iRun)
		tempMat._elemData[iRun] = this->_elemData[iRun] - rhs._elemData[iRun];
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

	for (size_t iRun = 0; iRun < _nRow * _nCol; ++iRun)
		tempMat._elemData[iRun] = this->_elemData[iRun] * rhs._elemData[iRun];
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

	for (size_t iRun = 0; iRun < _nRow * _nCol; ++iRun)
		tempMat._elemData[iRun] = this->_elemData[iRun] / rhs._elemData[iRun];
	return tempMat;
}

// Matrix arithmetic
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix& rhs) 
// Matrix multiplication
// Reference1: What every programmer needs to know about memory (https://akkadia.org/drepper/cpumemory.pdf)
// Reference2: Strassen algorithm (https://en.wikipedia.org/wiki/Strassen_algorithm)
{
	// Size ID
	size_t nRow0 = this->_nRow;
	size_t nSum0 = this->_nCol;
	size_t nCol0 = rhs._nCol;
	Matrix<T> tempMat{ nRow0,nCol0 };

	// Validity check
	if (this->_nCol != rhs._nRow) {
		tempMat._validity = false;
		return tempMat;
	}
	
	// Strassen algorithm for multiplication between square matrices
	const int subMatrixSize = 250; // submatrix split
	if ((nRow0> subMatrixSize) && (this->_nRow == this->_nCol) && (rhs._nRow == rhs._nCol)){
		// Adjustment to 2^n and 2^n
		size_t nRow2 = (nRow0 % 2 == 0) ? nRow0 : nRow0 + 1;
		size_t nCol2 = nRow2; // Square matrix
		size_t nCut  = nRow2 / 2;

		Matrix<T> A{ nRow2,nCol2 };
		Matrix<T> B{ nRow2,nCol2 };
		for (size_t iRow = 0; iRow < nRow2; ++iRow) {
		for (size_t iCol = 0; iCol < nCol2; ++iCol) {
			if ((iRow < nRow0) && (iCol < nCol0)) {
				A._elemData[iRow * nCol2 + iCol] = this->_elemData[iRow * this->_nCol + iCol];
				B._elemData[iRow * nCol2 + iCol] = rhs._elemData[  iRow * rhs._nCol   + iCol];
			}
			else {
				A._elemData[iRow * nCol2 + iCol] = 0;
				B._elemData[iRow * nCol2 + iCol] = 0;
			}
		}
		}
		Matrix<T> A11 = A(   0, nCut  - 1,    0, nCut  - 1);
		Matrix<T> A12 = A(   0, nCut  - 1, nCut, nCol2 - 1);
		Matrix<T> A21 = A(nCut, nRow2 - 1,    0, nCut  - 1);
		Matrix<T> A22 = A(nCut, nRow2 - 1, nCut, nCol2 - 1);
		Matrix<T> B11 = B(   0, nCut  - 1,    0, nCut  - 1);
		Matrix<T> B12 = B(   0, nCut  - 1, nCut, nCol2 - 1);
		Matrix<T> B21 = B(nCut, nRow2 - 1,    0, nCut  - 1);
		Matrix<T> B22 = B(nCut, nRow2 - 1, nCut, nCol2 - 1);

		Matrix<T> M1 = (A11 + A22) * (B11 + B22);
		Matrix<T> M2 = (A21 + A22) * B11;
		Matrix<T> M3 = A11 * (B12 - B22);
		Matrix<T> M4 = A22 * (B21 - B11);
		Matrix<T> M5 = (A11 + A12) * B22;
		Matrix<T> M6 = (A21 - A11) * (B11 + B12);
		Matrix<T> M7 = (A12 - A22) * (B21 + B22);

		const Matrix<T> C11 = M1 + M4 - M5 + M7;
		const Matrix<T> C12 = M3 + M5;
		const Matrix<T> C21 = M2 + M4;
		const Matrix<T> C22 = M1 - M2 + M3 + M6;

		Matrix<T> sum = Matrix<T>::ColMerge(Matrix<T>::RowMerge(C11, C21), Matrix<T>::RowMerge(C12, C22));
		return sum(0,nRow0-1,0,nCol0-1);
	}

	// Transpose rhs to optimize cache access
	Matrix<T> rhsT = rhs;
	rhsT = rhsT.Transpose();

	// Use split summation to aid compiler for instruction pipelining
	const size_t sumSplit = 4; // Increasing beyond 4 is not effective (1:1.34, 2:0.75, 4: 0.64, 5: 0.68, 10: 0.70)
	T tempSum0, tempSum1, tempSum2, tempSum3;
	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
	for (size_t iCol = 0; iCol < nCol0; ++iCol) {
		tempSum0 = 0; tempSum1 = 0; tempSum2 = 0; tempSum3 = 0;
		for (size_t iSum = 0; iSum < nSum0/sumSplit; ++iSum) {
			tempSum0 += this->_elemData[iRow * this->_nCol + iSum * sumSplit + 0] * rhsT._elemData[iCol * rhsT._nCol + iSum * sumSplit + 0];
			tempSum1 += this->_elemData[iRow * this->_nCol + iSum * sumSplit + 1] * rhsT._elemData[iCol * rhsT._nCol + iSum * sumSplit + 1];
			tempSum2 += this->_elemData[iRow * this->_nCol + iSum * sumSplit + 2] * rhsT._elemData[iCol * rhsT._nCol + iSum * sumSplit + 2];
			tempSum3 += this->_elemData[iRow * this->_nCol + iSum * sumSplit + 3] * rhsT._elemData[iCol * rhsT._nCol + iSum * sumSplit + 3];
		}
		for (size_t iSum = (nSum0 / sumSplit)*sumSplit; iSum < nSum0; ++iSum)
			tempSum0 += this->_elemData[iRow * this->_nCol + iSum] * rhsT._elemData[iCol * rhsT._nCol + iSum ];
		tempMat._elemData[iRow * nCol0 + iCol] = tempSum0 +tempSum2 + tempSum3;
	}
	}
	
	return tempMat;
}


template <typename T>
Matrix<T> Matrix<T>::BackSlash(const Matrix& rhs) // A\b (== inv(A)*b) is implemented as A.backSlash(b)
{
	size_t nRow0 = this->_nRow;
	size_t nSum0 = this->_nCol;
	size_t nCol0 = rhs._nCol;

	Matrix<T> tempMat{ nRow0,nCol0 };
	if (this->_nRow != rhs._nRow) {
		tempMat._validity = false;
		return tempMat;
	}

	// To be worked on further

	return tempMat;
}

/* To be worked on
template <typename T>
Matrix<T> Matrix<T>::Inverse() // Inverse of matrix
*/

template <typename T>
Matrix<T> Matrix<T>::Transpose() // Matrix transpose
{
	size_t nRow0 = this->_nCol;
	size_t nCol0 = this->_nRow;
	Matrix<T> tempMat{ nRow0,nCol0 };

	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
	for (size_t iCol = 0; iCol < nCol0; ++iCol) {
		tempMat._elemData[iRow * _nCol + iCol] = this->_elemData[iCol * _nCol + iRow];
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
	_elemData.resize(_nRow*_nCol, 0);
	//_elemData = new T[_nRow * _nCol];
}

template <typename T>
bool Matrix<T>::IsValid()
{
	return _validity;
}

template <typename T>
Matrix<T> Matrix<T>::RowMerge(const Matrix& lhs, const Matrix& rhs)
{
	size_t nRow0 = lhs._nRow + rhs._nRow;
	size_t nCol0 = lhs._nCol;
	Matrix<T> tempMat{ nRow0,nCol0 };
	if (lhs._nCol != rhs._nCol) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
		for (size_t iCol = 0; iCol < nCol0; ++iCol) {
			if (iRow < lhs._nRow){
				const size_t iRowLhs = iRow;
				tempMat._elemData[iRow * nCol0 + iCol] = lhs._elemData[iRowLhs * lhs._nCol + iCol];
			}
			else
			{
				const size_t iRowRhs = iRow-lhs._nRow;
				tempMat._elemData[iRow * nCol0 + iCol] = rhs._elemData[iRowRhs * lhs._nCol + iCol];
			}
		}
	}
	return tempMat;
}

template <typename T>
Matrix<T> Matrix<T>::ColMerge(const Matrix& lhs, const Matrix& rhs)
{
	size_t nRow0 = lhs._nRow;
	size_t nCol0 = lhs._nCol + rhs._nCol;
	Matrix<T> tempMat{ nRow0,nCol0 };
	if (lhs._nRow != rhs._nRow) {
		tempMat._validity = false;
		return tempMat;
	}

	for (size_t iRow = 0; iRow < nRow0; ++iRow) {
		for (size_t iCol = 0; iCol < nCol0; ++iCol) {
			if (iCol < lhs._nCol) {
				const size_t iColLhs = iCol;
				tempMat._elemData[iRow * nCol0 + iCol] = lhs._elemData[iRow * lhs._nCol + iColLhs];
			}
			else
			{
				const size_t iColRhs = iCol - lhs._nCol;
				tempMat._elemData[iRow * nCol0 + iCol] = rhs._elemData[iRow * lhs._nCol + iColRhs];
			}
		}
	}
	return tempMat;
}