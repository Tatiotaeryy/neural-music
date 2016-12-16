#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <malloc.h>

template <typename T>
class Matrix
{
private:
	T ** Data;
	int size_a, size_b;

	void free_data()
	{
		if ((unsigned int)Data != 0xcdcdcdcd) {
			for (int i = 0; i < size_a; ++i)
			{
				delete[] Data[i];
			}
			delete[] Data;
		}
	}
public:
	Matrix(int size_a, int  size_b) : size_a(size_a), size_b(size_b)
	{
		Data = new T*[size_a];
		for (int i = 0; i < size_a; ++i)
		{
			Data[i] = new T[size_b];
			for (int j = 0; j < size_b; ++j)
			{
				Data[i][j] = 0;
			}
		}
	}
	Matrix(int size) : size_a(size), size_b(size) {
		Data = new T*[size];
		for (int i = 0; i < size; ++i)
		{
			Data[i] = new T[size];
			for (int j = 0; j < size; ++j)
			{
				Data[i][j] = 0;
			}
		}
	}
	//	Matrix() { size_a = 0; size_b = 0; };
	Matrix(Matrix &B) : size_a(B.get_size_a()), size_b(B.get_size_b()) {
		Data = new T*[size_a];
		for (int i = 0; i < size_a; ++i)
		{
			Data[i] = new T[size_b];
			for (int j = 0; j < size_b; ++j)
			{
				Data[i][j] = B.get(i, j);
			}
		}
	}
	~Matrix() {
		free_data();
	}
	void PrintMatrix() const {
		for (int i = 0; i < size_a; ++i) {
			for (int j = 0; j < size_b; ++j)
			{
				std::cout << Data[i][j] << " ";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}
	void ReadMatrix(std::ifstream &input_file) {
		for (int i = 0; i < size_a; ++i)
		{
			for (int j = 0; j < size_b; ++j)
			{
				input_file >> Data[i][j];
			}
		}
	};
	void ReadRow(std::ifstream &input_file, int i) 
		{
		if (i < 0 || i >= size_a)
		{
			std::cout << "Wrong index in ReadRow()\n";
			return;
		}
			for (int j = 0; j < size_b; ++j)
			{
				input_file >> Data[i][j];
			}
		};
	T get(int i, int j) const
	{
		if (i < 0 || j < 0 || i >= size_a || j >= size_b)
		{
			std::cout << "Wrong index in get()\n";
			return NAN;
		}
		return Data[i][j];
	}
	void set(T value, int i, int j)
	{
		if (i < 0 || j < 0 || i >= size_a || j >= size_b)
		{
			std::cout << "Wrong index in set()\n";
			return;
		}
		Data[i][j] = value;
	}
	void add(T value, int i, int j)
	{
		if (i < 0 || j < 0 || i >= size_a || j >= size_b)
		{
			std::cout << "Wrong index in add()\n";
			return;
		}
		Data[i][j] += value;
	}
	inline int get_size_a() const
	{
		return size_a;
	}
	inline int get_size_b() const
	{
		return size_b;
	}

	Matrix Transpose() const
	{
		Matrix <T> A(size_b, size_a);
		for (int i = 0; i < size_a; ++i)
			for (int j = 0; j < size_b; ++j)
				A.set(Data[i][j], j, i);
		return A;
	}
	Matrix& operator= (Matrix &B)
	{
		free_data();
		size_a = B.get_size_a();
		size_b = B.get_size_b();
		Data = new T*[size_a];
		for (int i = 0; i < size_a; ++i)
		{
			Data[i] = new T[size_b];
			for (int j = 0; j < size_b; ++j)
			{
				Data[i][j] = B.get(i, j);
			}
		}
		return *this;
	}
	Matrix operator* (Matrix &B) const
	{
		if (size_b == B.get_size_a())
		{
			Matrix <T> A(size_a, B.get_size_b());
			for (int i = 0; i < size_a; ++i)
			{
				for (int j = 0; j < B.get_size_b(); ++j)
				{
					for (int k = 0; k < size_b; ++k)
					{
						A.add(Data[i][k] * B.get(k, j), i, j);
					}
				}
			}
			return A;
		}
		std::cout << "ERROR OCCURED WHILE MULTIPLYING MATRIXES\n";
		return 0;
	}
	Matrix operator* (T B) const
	{
		Matrix <T> A(size_a, size_b);
		for (int i = 0; i < size_a; ++i)
		{
			for (int j = 0; j < size_b; ++j)
			{
				A.set(Data[i][j] * B, i, j);
			}
		}
		return A;
	}
	Matrix operator+ (Matrix &B) const
	{
		if (size_a == B.get_size_a() && size_b == B.get_size_b()) {
			Matrix <T> A(size_a, size_b);
			for (int i = 0; i < size_a; ++i)
			{
				for (int j = 0; j < size_b; ++j)
				{
					A.set(Data[i][j] + B.get(i, j), i, j);
				}
			}
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator+ (T B) const {
		if (size_b == size_a) {
			Matrix <T> A(size_a, size_b);
			A = *this;
			for (int i = 0; i < size_a; ++i)
			{
				A.set(Data[i][i] + B, i, i);
			}
			return A;
		}
		std::cout << "ERROR OCCURED WHILE ADDING MATRIXES\n";
		return 0;
	}
	Matrix operator- (Matrix &B) const
	{
		if (size_b == B.get_size_b() && size_a == B.get_size_a())
		{
			Matrix <T> A(size_a, size_b);
			for (int i = 0; i < size_a; ++i)
			{
				for (int j = 0; j < size_b; ++j)
				{
					A.set(Data[i][j] - B.get(i, j), i, j);
				}
			}
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix operator- (T B) const
	{
		if (size_b == size_a)
		{
			Matrix <T> A(size_a, size_b);
			A = *this;
			for (int i = 0; i < size_a; ++i)
			{
				A.set(Data[i][i] - B, i, i);
			}
			return A;
		}
		std::cout << "ERROR OCCURED WHILE SUBSTRACTING MATRIXES\n";
		return 0;
	}
	Matrix GetRow(int i) const {
		Matrix <T> Row(1, size_b);
		for (int j = 0; j < size_b; ++j)
		{
			Row.set(Data[i][j], 0, j);
		}
		//std::cout << "Row " << i << "\n";
		//Row.PrintMatrix();
		return Row;
	}
	Matrix GetColumn(int i) const
	{
		Matrix <T> Column(size_a, 1);
		for (int j = 0; j < size_a; ++j)
		{
			Column.set(Data[j][i], j, 0);
		}
		//std::cout << "Column " << i << "\n";
		//Column.PrintMatrix();
		return Column;
	}
	Matrix GetSubMatrix(int i, int j) const
	{
		Matrix <T> SubMatrix(size_a, j - i);
		for (int k = 0; k < size_a; ++k)
		{
			for (int l = i; l < j; ++l)
			{
				SubMatrix.set(Data[k][l], k, l - i);
			}
		}
		return SubMatrix;
	}
	//ПЕРЕПИШИ
	T product(Matrix &X, Matrix &Y) {
		T prod = 0;
		if (X.get_size_a() == Y.get_size_a() && X.get_size_b() == Y.get_size_b())
			for (int i = 0; i < X.get_size_a(); ++i)
				for (int j = 0; j < X.get_size_b(); ++j)
					prod += X.get(i, j)*Y.get(i, j);
		else
			if (X.get_size_b() == Y.get_size_a() && X.get_size_a() == Y.get_size_b())
				for (int i = 0; i < X.get_size_a(); ++i)
					for (int j = 0; j < X.get_size_b(); ++j)
						prod += X.get(i, j)*Y.get(j, i);
			else {
				std::cout << "ERROR IN CALCULATING SCALAR PRODUCT\n";
				return NAN;
			}
			return prod;
	}
	//НЕЗАБУДЬ
	T vector_norm()
	{
		return sqrt(product(*this, *this));
	}
};
class DATA {
public:
	int num_of_examples;
	int num_of_variables;
	Matrix<double> X;
	Matrix<double> Y;
	DATA(int a, int b, Matrix <double> &X, Matrix <double> &Y) :
		num_of_examples(a), num_of_variables(b), X(X), Y(Y) {};
	DATA(DATA &A) :
		num_of_examples(A.num_of_examples),
		num_of_variables(A.num_of_variables),
		X(A.X), Y(A.Y) {};
	DATA(std::ifstream &input_file, int num_of_examples, int num_of_variables) :
		num_of_examples(num_of_examples), num_of_variables(num_of_variables),
		X(num_of_examples, num_of_variables), Y(num_of_examples, 1)
	{
		for (int i = 0; i < num_of_examples; ++i)
		{
			X.ReadRow(input_file, i);
			Y.ReadRow(input_file, i);
		}
	}
	void ReadData(std::ifstream &input_file, int num_of_examples, int num_of_variables) {
		Matrix <double> NewX(num_of_examples, num_of_variables);
		Matrix <double> NewY(num_of_examples, 1); //НОВОГОДНЯЯ ПЕРЕМЕННАЯ!  =)
		for (int i = 0; i < num_of_examples; ++i)
		{
			NewX.ReadRow(input_file, i);
			NewY.ReadRow(input_file, i);
		}
		X = NewX;
		Y = NewY;
	}
	DATA operator=(DATA &A)
	{
		return DATA(A);
	}
	~DATA() {
		X.Matrix::~Matrix();
		Y.Matrix::~Matrix();
	}

};
class Neuron
{
protected:
	Matrix<double> Theta;
	Matrix<double> * Input;
	int size;
public:
	Neuron(int number_of_variables, Matrix<double> &Input) : Theta(number_of_variables, 1), size(number_of_variables), Input(new Matrix<double>(Input))
	{
		for (int i = 0; i < Theta.get_size_a(); ++i)
		{
			Theta.set(1, i, 0);
		}
	}
	virtual double H(int i) = 0;
	virtual double CostFunction(Matrix<double> &Y, double lambda = 0) = 0;
	virtual void ChangeInput(Matrix<double> &Input) = 0;
	void GradientDescent(Matrix<double> &Y, double alpha, double lambda = 0)
	{
		double buffer = 1;
		int k = 0;
		Matrix <double> Gradient_vector(Theta.get_size_a(), 1);
		while (k < 10000 && abs(buffer) > 0.00001)
		{
			//std::cout << "k is " << k<<"\n";
			buffer = Theta.vector_norm();
			for (int i = 0; i < Theta.get_size_a(); ++i)
			{
				for (int j = 0; j < Input->get_size_a(); ++j)
				{
					Gradient_vector.add((H(i) - Y.get(j, 0))*(Input->get(j, i)), i, 0);
				}
				Gradient_vector.add(Theta.get(i, 0) * lambda, i, 0);
				Gradient_vector = Gradient_vector * (alpha / Theta.get_size_a());
			}
			Theta = Theta - Gradient_vector;
			buffer -= Theta.vector_norm();
			++k;
		}
		std::cout << "Theta is" << "\n";
		Theta.PrintMatrix();
	}
};
class RegressionNeuron : public Neuron
{
public:
	RegressionNeuron(int number_of_variables, Matrix<double> &Input) : Neuron(number_of_variables, Input) {}
	double H(int i) {
		return Input->product(Input->GetRow(i), Theta);
	}
	double CostFunction(Matrix<double> &Y, double lambda = 0) {
		double summ = 0;
		for (int i = 0; i < Input->get_size_a(); ++i)
		{
			summ += pow(H(i) - Y.get(i, 0), 2);
		}
		//Regularization part
		for (int i = 0; i < Theta.get_size_a(); ++i)
		{
			summ += pow(Theta.get(i, 0), 2) * lambda;
		}

		summ /= Input->get_size_a() * 2;
		return summ;
	}
};
class LogisticNeuron : public Neuron
{
public:
	LogisticNeuron(int number_of_variables, Matrix<double> &Input) : Neuron(number_of_variables, Input) {}
	double H(int i) {
		return 1 / (1 + exp(-Input->product(Input->GetRow(i), Theta)));
	}
	double CostFunction(Matrix<double> &Y, double lambda = 0) {
		double summ = 0;
		double h = 0;
		for (int i = 0; i < Input->get_size_a(); ++i) {
			h = H(i);
			summ += Y.get(i, 0) * log(h) + (1 - Y.get(i, 0)) * log(1 - h);
		}
		//Regularization part
		for (int i = 0; i < Theta.get_size_a(); ++i)  summ -= pow(Theta.get(i, 0), 2) * lambda / 2;
		return -summ / Input->get_size_a();

	}
	void ChangeInput(Matrix<double> &Input)
	{
		delete this->Input;
		this->Input = new Matrix<double>(Input);
	}
};

class NeuronLayer
{
	int num_of_examples;
	int num_of_neurons;
	DATA INPUT;
	Matrix <double> Output;
	Neuron ** Data;
public:
	NeuronLayer(int num_of_examples, int num_of_neurons, Matrix<double> &X, Matrix<double> &Y) :
		INPUT(num_of_examples, X.get_size_b(), X, Y),
		Output(num_of_examples, num_of_neurons),
		num_of_examples(num_of_examples),
		num_of_neurons(num_of_neurons)
	{
		Data = (Neuron **)malloc(sizeof(Neuron*)*num_of_neurons);
		//INPUT.X = INPUT.X*(1 / 9.0);
	}
	NeuronLayer(int num_of_examples, int num_of_neurons, DATA &A) :
		INPUT(A),
		Output(num_of_examples, num_of_neurons),
		num_of_examples(num_of_examples),
		num_of_neurons(num_of_neurons)
	{
		Data = (Neuron **)malloc(sizeof(Neuron*)*num_of_neurons);
	}
	void SetLayer(std::ifstream &input_file) {
		int input_size, border_left, border_right;
		input_size = border_left = border_right = 0;
		for (int i = 0; i < num_of_neurons; ++i) {
			input_file >> input_size;
			input_file >> border_left;
			input_file >> border_right;
			Data[i] = new LogisticNeuron(input_size, INPUT.X.GetSubMatrix(border_left, border_right));
		}
	}
	void SetupLayer(double alpha, double lambda)
	{
		for (int i = 0; i < num_of_neurons; ++i)
		{
			Data[i]->GradientDescent(INPUT.Y, alpha, lambda);
		}
	}
	Matrix<double>& get_output()
	{
		return Output;
	}
	void ChangeInput(Matrix<double> &Input)
	{
		std::ifstream input_file("Neuron_Setting.txt");
		int input_size, border_left, border_right;
		input_size = border_left = border_right = 0;
		INPUT.X = Matrix<double>(INPUT.X*(1 / 9.0));
		for (int i = 0; i < num_of_neurons; ++i)
		{
			input_file >> input_size;
			input_file >> border_left;
			input_file >> border_right;
			Data[i]->ChangeInput(this->INPUT.X.GetSubMatrix(border_left, border_right));
		}
		num_of_examples = Input.get_size_a();
		Output = Matrix<double>(num_of_examples, num_of_neurons);
	}
	void Work()
	{
		for (int i = 0; i < num_of_examples; ++i)
		{
			for (int j = 0; j < num_of_neurons; ++j)
			{
				if (Data[j]->H(i) > 0.5)
				{
					Output.set(Data[j]->H(i), i, j);
				}
				else
				{
					Output.set(Data[j]->H(i), i, j);
				}
			}
		}
	}
};
int main() {
	//Ввод обучающих данных
	std::ifstream input_file("text.txt");
	std::ifstream input_neuron("Neuron_Setting.txt");
	int number_of_examples;
	int number_of_variables;
	double alpha, lambda;
	input_file >> number_of_examples;
	input_file >> number_of_variables;
	input_file >> alpha;
	input_file >> lambda;
	DATA MAIN_DATA(input_file, number_of_examples, number_of_variables);

	//Создание и настройка слоя
	NeuronLayer MAIN_LAYER(number_of_examples, 3, MAIN_DATA);
	MAIN_LAYER.SetLayer(input_neuron);
	MAIN_LAYER.SetupLayer(alpha, lambda);
	MAIN_LAYER.Work();
	LogisticNeuron exit(3, MAIN_LAYER.get_output());
	std::cout << "Output is\n";
	MAIN_LAYER.get_output().PrintMatrix();
	exit.GradientDescent(MAIN_DATA.Y, alpha, lambda);
	Matrix <double> test_set(1, number_of_variables);
	int in1, in2, in3;
	while (true)
	{
		std::cout << "Perinita pelease\n";
		std::cin >> in1 >> in2 >> in3;
		test_set.set(in1, 0, 0);
		test_set.set(in2, 0, 1);
		test_set.set(in3, 0, 2);
		MAIN_LAYER.ChangeInput(test_set);
		MAIN_LAYER.Work();
		MAIN_LAYER.get_output().PrintMatrix();
		exit.ChangeInput(MAIN_LAYER.get_output());
		std::cout << "Hypothesis for " << in1 << in2 << in3 << " is " << exit.H(0) << "\n";
		system("pause");
	}

	//Theta.ReadMatrix(input_file);
	//Theta.PrintMatrix();
	//input_file >> lambda;
	system("pause");
	return 0;
}
