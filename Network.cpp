#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

template <typename T> class Matrix {
private:
	T ** Data;
	int rows, cols;
	Matrix() {};
public:
	Matrix(int row_size, int column_size) {
		rows = row_size;
		cols = column_size;
		Data = new T*[rows];
		for (int i = 0; i < rows; ++i)
			Data[i] = new T[cols * sizeof(T)];
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				Data[i][j] = 0;
	}
	Matrix(int size) {
		rows = size;
		cols = size;
		Data = new T*[rows];
		for (int i = 0; i < rows; ++i)
			Data[i] = new T[cols * sizeof(T)];
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				Data[i][j] = 0;
	}
	void PrintMatrix() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j)
				std::cout << Data[i][j] << " ";
			std::cout << "\n";
		}
		std::cout << "\n";

	}
	void SetMatrix(std::ifstream &input_file) {
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				input_file >> Data[i][j];
	};
	T get(int i, int j) { return Data[i][j]; }
	void set(T value, int i, int j) { Data[i][j] = value; }
	void add(T value, int i, int j) { Data[i][j] += value; }
	int Rows() { return rows; }
	int Cols() { return cols; }
	Matrix Transpose() {
		Matrix <T> A(cols, rows);
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				A.set(Data[i][j], j, i);
		return A;
	}
	void operator= (Matrix &B) {
		rows = B.Rows();
		cols = B.Cols();
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				Data[i][j] = B.get(i, j);
	}
	Matrix operator* (Matrix &B){
		if (cols == B.Rows()) {
			Matrix <T> A(rows, B.Cols());
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < B.Cols(); ++j)
					for (int k = 0; k < cols; ++k)
						A.add(Data[k][j] + B.get(i, k), i, j);
			return A;
		}
	std::cout << "ERROR OCCURED WHILE MULTIPLYING MATRIXES\n";
	return 0;
}
	Matrix operator* (T B) {
			Matrix <T> A(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					A.set(Data[i][j]*B, i, j);
			return A;
		}
	Matrix operator+ (Matrix &B) {
		if (cols == B.Cols() && rows == B.Rows()) {
			Matrix <T> A(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					A.set(Data[i][j] + B.get(i, j), i ,j);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator+ (T B) {
		if (cols == rows) {
			Matrix <T> A(rows, cols);
			A = *this;
			for (int i = 0; i < rows; ++i)
					A.set(Data[i][i] + B, i, i);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator- (Matrix &B) {
		if (cols == B.Cols() && rows == B.Rows()) {
			Matrix <T> A(rows, cols);
			for (int i = 0; i < rows; ++i)
				for (int j = 0; j < cols; ++j)
					A.set(Data[i][j] - B.get(i, j), i, j);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix operator- (T B) {
		if (cols == rows) {
			Matrix <T> A(rows, cols);
			A = *this;
			for (int i = 0; i < rows; ++i)
				A.set(Data[i][i] - B, i, i);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix GetRow(int i) {
		Matrix <T> Row(1, cols);
		for (int j = 0; j < cols; ++j) Row.setdataitem(data[i][j], 0, j);
		/*std::cout << "Column " << i << "\n";
		Column.printMatrix();*/
		return Row;
	}
	T product(Matrix &X, Matrix &Y) {
		T prod = 0;
		if (X.s)
			for (int i = 0; i < X.Rows(); ++i)
				for (int j = 0; i < X.Cols(); ++i)
					for (int k = 0; i < Y.Rows(); ++i)
						for (int l = 0; i < Y.Cols(); ++i)
							prod += X.get(i, j) * Y.get(k, l);
		return prod;
	}
};
class RegressionNeuron{
public:
	double CostFunction(Matrix<double> &X, Matrix<double> &Y, Matrix<double> &Theta, double lambda) {
		double summ = 0;
		for (int i = 1; i < X.Cols(); ++i) summ += pow(X.product(X.GetRow(i), Theta) - Y.get(i, 0), 2);
		return summ / X.Cols() / 2;

	}
};
class LogisticNeuron {
public:
	double CostFunction(Matrix<double> &X, Matrix<double> &Y, Matrix<double> &Theta, double lambda) {
		double summ = 0;
		double h = 1 / (1 + exp(X.product(X.GetRow(i), Theta)));
		for (int i = 1; i < X.Cols(); ++i) summ += Y.get(i, 0) * log(h) + (1 - Y.get(i, 0)) * log (1 -h);
		return -summ / X.Cols();

	}
};
int main() {
	std::ifstream input_file("text.txt");
	
	

	system("pause");
	return 0;
}
