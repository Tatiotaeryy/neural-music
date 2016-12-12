#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

template <typename T> class Matrix {
private:
	T ** Data;
	int size_a, size_b;
	Matrix() {};
public:
	Matrix(int size_a, int column_size) {
		this->size_a = size_a;
		this->size_b = column_size;
		Data = new T*[size_a];
		for (int i = 0; i < size_a; ++i)
			Data[i] = new T[size_b * sizeof(T)];
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				Data[i][j] = 0;
	}
	Matrix(int size) {
		size_a = size;
		size_b = size;
		Data = new T*[size_a];
		for (int i = 0; i < size_a; ++i)
			Data[i] = new T[size_b * sizeof(T)];
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				Data[i][j] = 0;
	}
	void PrintMatrix() {
		for (int i = 0; i < size_a; ++i) {
			for (int j = 0; j < size_b; ++j)
				std::cout << Data[i][j] << " ";
			std::cout << "\n";
		}
		std::cout << "\n";

	}
	void SetMatrix(std::ifstream &input_file) {
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				input_file >> Data[i][j];
	};
	T get(int i, int j) { return Data[i][j]; }
	void set(T value, int i, int j) { Data[i][j] = value; }
	void add(T value, int i, int j) { Data[i][j] += value; }
	int get_size_a() { return size_a; }
	int get_size_b() { return size_b; }
	Matrix Transpose() {
		Matrix <T> A(size_b, size_a);
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				A.set(Data[i][j], j, i);
		return A;
	}
	void operator= (Matrix &B) {
		size_a = B.get_size_a();
		size_b = B.get_size_b();
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				Data[i][j] = B.get(i, j);
	}
	Matrix operator* (Matrix &B){
		if (size_b == B.get_size_b()) {
			Matrix <T> A(size_a, B.get_size_a());
			for (int i = 0; i < size_a; ++i)
				for (int j = 0; j < B.get_size_a(); ++j)
					for (int k = 0; k < size_b; ++k)
						A.add(Data[k][j] + B.get(i, k), i, j);
			return A;
		}
	std::cout << "ERROR OCCURED WHILE MULTIPLYING MATRIXES\n";
	return 0;
}
	Matrix operator* (T B) {
			Matrix <T> A(size_a, size_b);
			for (int i = 0; i < size_a; ++i)
				for (int j = 0; j < size_b; ++j)
					A.set(Data[i][j]*B, i, j);
			return A;
		}
	Matrix operator+ (Matrix &B) {
		if (size_b == B.get_size_a() && size_a == B.get_size_b()) {
			Matrix <T> A(size_a, size_b);
			for (int i = 0; i < size_a; ++i)
				for (int j = 0; j < size_b; ++j)
					A.set(Data[i][j] + B.get(i, j), i ,j);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator+ (T B) {
		if (size_b == size_a) {
			Matrix <T> A(size_a, size_b);
			A = *this;
			for (int i = 0; i < size_a; ++i)
					A.set(Data[i][i] + B, i, i);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator- (Matrix &B) {
		if (size_b == B.get_size_a() && size_a == B.get_size_b()) {
			Matrix <T> A(size_a, size_b);
			for (int i = 0; i < size_a; ++i)
				for (int j = 0; j < size_b; ++j)
					A.set(Data[i][j] - B.get(i, j), i, j);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix operator- (T B) {
		if (size_b == size_a) {
			Matrix <T> A(size_a, size_b);
			A = *this;
			for (int i = 0; i < size_a; ++i)
				A.set(Data[i][i] - B, i, i);
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix GetRow(int i) {
		Matrix <T> Row(1, size_b);
		for (int j = 0; j < size_b; ++j) Row.set(Data[i][j], 0, j);
		//std::cout << "Row " << i << "\n";
		//Row.PrintMatrix();
		return Row;
	}
	Matrix GetColumn(int i) {
		Matrix <T> Column(size_a, 1);
		for (int j = 0; j < size_a; ++j) Column.set(Data[j][i], j, 0);
		//std::cout << "Column " << i << "\n";
		//Column.PrintMatrix();
		return Column;
	}
	T product(Matrix &X, Matrix &Y) {
		T prod = 0;
		for (int i = 0; i < X.get_size_a(); ++i)
			for (int j = 0; j < X.get_size_b(); ++j)
				for (int k = 0; k < Y.get_size_a(); ++k)
					for (int l = 0; l < Y.get_size_b(); ++l)
						prod += X.get(i, j) * Y.get(k, l);
		//std::cout << "Product is " << prod << "\n";
		return prod;
	}
};
class RegressionNeuron{
public:
	double CostFunction(Matrix<double> &X, Matrix<double> &Y, Matrix<double> &Theta, double lambda) {
		double summ = 0;
		for (int i = 0; i < X.get_size_a(); ++i) summ += pow(X.product(X.GetRow(i), Theta) - Y.get(i, 0), 2);
		return summ / X.get_size_a() / 2;

	}
};
class LogisticNeuron {
public:
	double CostFunction(Matrix<double> &X, Matrix<double> &Y, Matrix<double> &Theta, double lambda) {
		double summ = 0;
		
		for (int i = 0; i < X.get_size_a(); ++i) {
			double h = 1 / (1 + exp(X.product(X.GetRow(i), Theta)));
			summ += Y.get(i, 0) * log(h) + (1 - Y.get(i, 0)) * log(1 - h); }

		return -summ / X.get_size_a();

	}
};

int main() {
	std::ifstream input_file("text.txt");
	int number_of_examples;
	int number_of_variables;
	
	input_file >> number_of_examples;
	input_file >> number_of_variables;

	Matrix <double> X(number_of_examples, number_of_variables);
	Matrix <double> Y(number_of_examples, 1);
	Matrix <double> Theta(number_of_variables, 1);
	RegressionNeuron Test;
	LogisticNeuron Test2;
	double lambda = 0;

	X.SetMatrix(input_file);
	X.PrintMatrix();
	Y.SetMatrix(input_file);
	Y.PrintMatrix();
	Theta.SetMatrix(input_file);
	Theta.PrintMatrix();
	input_file >> lambda;

	std::cout<<"CostFunction for regression is "<<Test.CostFunction(X, Y, Theta, lambda)<<"\n";
	std::cout << "CostFunction for classification is " << Test2.CostFunction(X, Y, Theta, lambda)<<"\n";
	system("pause");
	return 0;
}
