#include "pch.h"
#include <iostream>
#include <list>
#include <string>
#include <fstream>
#include <time.h>
#include <algorithm>

void FU(int **D, int V)
{
	double tn = clock();
	int k;

	for (size_t k = 0; k < V; k++)
		for (size_t i = 0; i < V; i++)
			for (size_t j = 0; j < V; j++)
				if (D[i][k] && D[k][j] && i != j)
					if (D[i][k] + D[k][j] < D[i][j] || D[i][j] == 0)
						D[i][j] = D[i][k] + D[k][j];

	double tk = clock();

	std::string name = "out_matr_" + std::to_string(V) + ".txt";
	std::ofstream matrF(name);

	for (size_t i = 0; i < V; i++)
	{
		for (size_t j = 0; j < V; j++) matrF << D[i][j] << " ";
		matrF << std::endl;
	}
	matrF << std::endl << "Time: " << (tk - tn) / 1000.0;
	matrF.close();
}

void generate(int n, int **mas) {
	srand(time(NULL));
	std::string name = "matr_" + std::to_string(n) + ".txt";
	std::ofstream matrF(name);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i == j) {
				mas[i][j] = 0;
				matrF << 0 << " ";
			}
			else {
				int t = rand() % 50;
				mas[i][j] = t;
				matrF << t << " ";
			}
		}
		matrF << std::endl;
	}
	matrF.close();
}


int main()
{
	std::list<int> n{ 3, 100, 250, 500 , 1000 };
	for (int num : n) {
		int **mas = new int *[num];
		for (size_t i = 0; i < num; i++)
			mas[i] = new int[num];
		generate(num, mas);
		FU(mas, num);
		for (int i = 0; i < num; i++)
			delete[]mas[i];
		delete[]mas;
	}
}