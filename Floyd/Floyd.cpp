#include "pch.h"
#include <iostream>
#include <list>
#include <string>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <mpi.h>

constexpr auto SMINF = 999999;

void FU(int **D, int V, int **h)
{
	double tn = clock();
	int k;

	for (size_t i = 0; i < V; i++)
		for (size_t j = 0; j < V; j++)
			if (D[i][j] < 0)
				h[i][j] = 0;
			else
				h[i][j] = j + 1;


	for (size_t k = 0; k < V; k++)
		for (size_t i = 0; i < V; i++)
			for (size_t j = 0; j < V; j++)
				//if (i == j) {
				//	D[i][j] = 0;
				//}
				//else {
				//	if (D[i][j] == 0) D[i][j] = SMINF;
				//	D[i][j] = std::min(D[i][j], D[i][k] + D[k][j]);
				//	//h[i][j] = h[i][k];
				//}
				if (D[i][k] && D[k][j] && i != j)
					if (D[i][k] + D[k][j] < D[i][j] || D[i][j] == 0) {
						D[i][j] = D[i][k] + D[k][j];
						h[i][j] = h[i][k];
					}

	/*if (V == 10) {
		int u = 3, r = 1, ww = u - 1;
		while (ww != r - 1) {
			ww = h[ww][r - 1];
			std::cout << "Next vert - " + std::to_string(ww) << std::endl;
			ww--;
		}
	}*/

	double tk = clock();

	std::string name = "out_matr_" + std::to_string(V) + ".txt";
	std::string nameh = "out_matrH_" + std::to_string(V) + ".txt";
	std::ofstream matrF(name);
	std::ofstream matrH(nameh);

	for (size_t i = 0; i < V; i++)
	{
		for (size_t j = 0; j < V; j++) {
			matrF << D[i][j] << " ";
			matrH << h[i][j] << " ";
		}
		matrF << std::endl;
		matrH << std::endl;
	}
	matrF << std::endl << "Time: " << (tk - tn) / 1000.0;
	matrF.close();
	matrH.close();
}


void load(int **D, int V) {
	std::string name = "matr_" + std::to_string(V) + ".txt";
	std::ifstream matrF(name);

	for (int i = 0; i < V; i++)
		for (int j = 0; j < V; j++) {
			matrF >> D[i][j];
		}
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
	}
	matrF.close();
}

void FloydsAlgorithm(int rank, int *data, int N, int start, int count) {

	int k, i, j, k_here;
	int ij, ik;

	int *rowk = new int[N];

	for (k = 0; k < N; k++) {

		int owner = (int)ceil(k / count);
		// Если k принадлежит данному процессу - заполняем строку
		if (rank == owner) {
			for (j = 0; j < N; j++)
				rowk[j] = data[k*N + j];
		}

		MPI_Bcast(&k, 1, MPI_INT, owner, MPI_COMM_WORLD);

		MPI_Bcast(rowk, N, MPI_INT, owner, MPI_COMM_WORLD);
		for (i = start; i < start + count; i++) {
			for (j = 0; j < N; j++) {

				ij = i * N + j;
				ik = i * N + k;

				if (data[ik] && rowk[j] && i != j)
					if (data[ik] + rowk[j] < data[ij] || data[ij] == 0) {
						data[ij] = data[ik] + rowk[j];
					}
				/*if (i == j) {
					data[ij] = 0;
				}
				else {
					if (data[ij] == 0) data[ij] = SMINF;
					data[ij] = min(data[ij], data[ik] + rowk[j]);
				}*/
			}
		}
		int index;
		/*if (!rank && N == 10) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					index = i * N + j;
					std::cout << data[index] << ' ';
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}*/
	}
}

void Server(int size, std::string name, int N) {
	MPI_Status status;

	//ifstream M_in("matr_100.txt");
	std::ifstream M_in(name);
	int tmp, index;
	int tt = 0;
	//int N = 100;

	// Загружаем данные
	int *data = new int[N*N];
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			M_in >> tmp;
			data[i*N + j] = tmp;
		}
	M_in.close();

	/*for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			std::cout << data[i * N + j] << " ";
		}
		std::cout << endl;
	}*/

	// Передаем N
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// Передаем массив
	MPI_Bcast(data, N*N, MPI_INT, 0, MPI_COMM_WORLD);

	int count = (int)ceil(N / size);
	int start = N;
	int total = N * N;

	FloydsAlgorithm(0, data, N, 0, count);

	int *t = new int[N*N];
	for (int p = 1; p < size; p++) {
		MPI_Recv(t, N*N, MPI_INT, p, 0, MPI_COMM_WORLD, &status);
		for (int v = 0; v < total; v++) {
			data[v] = std::min(data[v], t[v]);
		}
	}

	std::string nameO = "parallel_out_matr_" + std::to_string(N) + ".txt";
	std::ofstream matrF(nameO);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			index = i * N + j;
			matrF << data[index] << ' ';
		}
		matrF << std::endl;
	}
}
// Slave процессы - запрос, Флойд, данные
void Slave(int rank, int S) {
	int N;
	MPI_Status status;

	// Получение N
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int num = N;
	int size = N * N;
	int *data = new int[size];
	//std::cout << N << endl;

	// Получение матрицы
	MPI_Bcast(data, N*N, MPI_INT, 0, MPI_COMM_WORLD);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) 
				data[i*N + j] = 0;
			else {
				if (data[i*N + j] == 0) 
					data[i*N + j] = 999;
			}
		}
	}

	/*if (rank == 4) {


		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				std::cout << data[i * N + j] << " ";
			}
			std::cout << endl;
		}
		std::cout << endl<<endl;
	}*/

	// CРасчет начала и конца полосы
	int count = (int)ceil(num / S);
	int start = rank * count;
	if ((num * start) + (num * count) > size) count = N - start;

	// Флойд
	FloydsAlgorithm(rank, data, num, start, count);

	/*if (rank == 4) {


		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				cout << data[i * N + j] << " ";
			}
			cout << endl;
		}
		cout << endl << endl;
	}*/

	// Рассылка данных на Master процесс
	MPI_Send(data, N*N, MPI_INT, 0, 0, MPI_COMM_WORLD);
}


int main(int argc, char* argv[])
{

	int size, rank, len;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::list<int> n{ 10,100, 250, 500, 1000 };
	if (!rank) {
		std::cout << "Sequential" << std::endl;
		for (int num : n) {
			double startTime = MPI_Wtime();
			int **mas = new int *[num];
			int **h = new int*[num];
			for (size_t i = 0; i < num; i++) {
				mas[i] = new int[num];
				h[i] = new int[num];
			}
			//generate(num, mas);
			load(mas, num);
			FU(mas, num, h);
			std::cout << "N = " << num << " Time: " << MPI_Wtime() - startTime << std::endl;
			for (int i = 0; i < num; i++) {
				delete[]mas[i];
				delete[]h[i];
			}
			delete[]mas;
			delete[]h;
		}
	}

	if (rank == 0)
		std::cout << std::endl << "Parallel" << std::endl;
	for (int N : n) {
		if (rank == 0)
		{
			std::string name = "matr_" + std::to_string(N) + ".txt";
			double startTime = MPI_Wtime();
			Server(size, name, N);
			std::cout << "N = " << N << " Time: " << MPI_Wtime() - startTime << std::endl;
		}
		else
		{
			Slave(rank, size);
		}
	}
	MPI_Finalize();
}