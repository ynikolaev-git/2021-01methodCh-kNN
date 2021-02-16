/*
 * c4prog3.cpp
 * Сравниваем данные полченные методом характеристик и методом ближайших соседей.
 * Сравниваем максимальные относительные и абсолютные отклонения инвариантов Римана.
 *
 * Created on: 26 октября 2020 г.
 *      Author: ЮН
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <vector>

#include "c4mx.hpp"
using namespace std;

const int maxPrecision = 24;

struct prn2 {
	int s;
	int c;
	double t;
	double x;
	double R;
	double L;
};

vector<prn2> net13, net23;

//========= Читаем данные известной сетки ===================================================
void readNet3(const char* fName, const int flag)
// flag - принимает значение 1 (читаем обучающую сеть) и 2 (сеть значения R и L в которой рассчитываем)
{
	ifstream input(fName);
	prn2 p;
	if (input) {
		cout << fName << " is opened" << endl;
	} else {
		cout << fName << " is not opened" << endl;
	}
	string rd;
	getline(input, rd);  // прочитали строку с названиями колонок
	while (getline(input, rd, ';')) {
		p.s = atof(rd.c_str());
		getline(input, rd, ';');
		p.c = atof(rd.c_str());
		getline(input, rd, ';');
		p.t = atof(rd.c_str());
		getline(input, rd, ';');
		p.x = atof(rd.c_str());
		getline(input, rd, ';');
		p.R = atof(rd.c_str());
		if (flag == 1){
			getline(input, rd);
			p.L = atof(rd.c_str());
			net13.push_back(p);
			}
		else if (flag == 2){
			getline(input, rd, ';');
			p.L = atof(rd.c_str());
			getline(input, rd);
			net23.push_back(p);
			}
	}
	if (flag == 1) cout << "прочитали из файла точки сети (данные из метода характеристик)" << endl;
	else if (flag == 2) cout << "прочитали из файла точки сети (обученные данные)" << endl;
}

//  далее основная процедура ======================================================================
int c4main3()
	{
	cout << endl<< endl<< "==== Сравненение результатов счета методом характеристик и методом ближайших соседей ========" << endl;
	cout << fixed << setprecision(maxPrecision);
	readNet3("netTest.csv", 1);
	readNet3("net4NN.csv", 2);
	// прочитали в переменные net13 и net23, переменные одинакового размера.
	double maxdR = 0, maxdL = 0;
	double dR, dL;
	int j1 = 0, j2 = 0, j5 = 0;
	int sizeNet1 = net13.size(), sizeNet2 = net23.size();
	if (sizeNet1 == sizeNet2)
		cout << "Размеры сеток равны " << sizeNet1 << endl;
	for (int i=0; i<sizeNet1;i++)
		{
		dR = abs((net13[i].R-net23[i].R)/net13[i].R)*100;
		dL = abs((net13[i].L-net23[i].L)/net13[i].L)*100;
		if (dR > 1) cout << "s= " << net13[i].s << "; c= " << net13[i].c << "; i= " << i << "; dR= " << dR << endl;
		if (dL > 1) cout << "s= " << net13[i].s << "; c= " << net13[i].c << "; i= " << i << "; dL= " << dL << endl;
		if (dR > maxdR) maxdR = dR;
		if (dL > maxdL)	maxdL = dL;
		if ((dR >1)||(dL >1)) j1++;
		if ((dR >2)||(dL >2)) j2++;
		if ((dR >5)||(dL >5)) j5++;
		}
	cout << "максимальное относительное отклонение R (%), maxdR= " << maxdR << endl;
	cout << "максимальное относительное отклонение L (%), maxdL= " << maxdL << endl;
	cout << "всего точек с относительным отклонением R или L больше 1% = " << j1 << endl;
	cout << "всего точек с относительным отклонением R или L больше 2% = " << j2 << endl;
	cout << "всего точек с относительным отклонением R или L больше 5% = " << j5 << endl;
	cout << "=== Сравнение завершено ===" << endl;
	return 0;
}

