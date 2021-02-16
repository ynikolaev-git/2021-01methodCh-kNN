/*
 * c4prog2.cpp
 * Модуль программы открывает файл 1 с ранее рассчитанной характеристической сеткой. Затем открывает файл 2 с точками
 * новой сетки и по данным файла 1 рассчитывает значения инвариантов R и L в точках файла 2.
 *
 * Created on: 15 октября 2020 г.
 *      Author: ЮН
 */

/*  !!! убрать комментарий здесь и в самом конце
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

const int maxLen2 = 200000;
const int maxPrecision = 24;
double p[maxLen2][4], Layer1[maxLen2][4], Layer2[maxLen2][4];   // t,x,R,L
double procDR = 0.01;  // в %
double procDL = 0.01;
//int kratnost = 100;

struct prn2 {
	double m;    // метрика = расстояние от точки сети до интерполируемой точки
	int s;
	int c;
	double t;
	double x;
	double R;
	double L;
	double DR_cur;
	double DL_cur;

};

vector<prn2> net1, netTest, net4NN;
vector<prn2> nb1(4);

//========= Читаем данные известной сетки ===================================================
void readNet1(const char* fName, const int flag)
// flag - принимает значение 1 (читаем обучающую сеть) и 2 (сеть значения R и L в которой рассчитываем)
{
	ifstream input(fName);
	prn2 p;
	if (input) {cout << fName << " is opened" << endl;}
	else {cout << fName << " is not opened" << endl;}
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
			net1.push_back(p);
			}
		else if (flag == 3){
			getline(input, rd, ';');
			p.L = atof(rd.c_str());
			getline(input, rd, ';');
			p.DR_cur = atof(rd.c_str());
			getline(input, rd);
			p.DL_cur = atof(rd.c_str());
			net4NN.push_back(p);
			}
	}
	if (flag == 1) cout << "прочитали файл точек сети, рассчитанной методом характеристик (net1.csv)" << endl;
	else if (flag == 3) cout << "прочитали файл точек сети, рассчитываемой через соседей (net4NN.csv)" << endl;
}

//=== Считаем значения газодинамических параметров в точке по значениям в точках-соседях ===
prn2 fourNN(const prn2 p, const vector<prn2> nb)
	{
	prn2 res;

	double m1, k;
	res.t = p.t; res.x = p.x;
	m1 = 0;
	if (nb[0].m != 1) m1 += nb[1].m*nb[2].m*nb[3].m;
	if (nb[1].m != 1) m1 += nb[0].m*nb[2].m*nb[3].m;
	if (nb[2].m != 1) m1 += nb[0].m*nb[1].m*nb[3].m;
	if (nb[3].m != 1) m1 += nb[0].m*nb[1].m*nb[2].m;

	k = (nb[1].m*nb[2].m*nb[3].m)/m1;
	res.R =  nb[0].R*k;
	res.L =  nb[0].L*k;
	k = (nb[0].m*nb[2].m*nb[3].m)/m1;
	res.R += nb[1].R*k;
	res.L += nb[1].L*k;
	k = (nb[0].m*nb[1].m*nb[3].m)/m1;
	res.R += nb[2].R*k;
	res.L += nb[2].L*k;
	k = (nb[0].m*nb[1].m*nb[2].m)/m1;
	res.R += nb[3].R*k;
	res.L += nb[3].L*k;
	res.DR_cur = abs((res.R-p.R)/res.R)*100;
	res.DL_cur = abs((res.L-p.L)/res.L)*100;
	return res;
	}

//=== Поиск ближайших соседей ===
vector<prn2> findNN(const prn2 p, const vector<prn2> net)
// p - точка, для которой ищем ближайших соседей среди точек вектора net
	{
	double m;
	vector<prn2> nb(4);
	nb[0].m = 1; nb[1].m = 1; nb[2].m = 1; nb[3].m = 1; // инициализировали значение расстояния
	nb[0].R = 0; nb[1].R = 0; nb[2].R = 0; nb[3].R = 0; // инициализировали значение расстояния
	nb[0].L = 0; nb[1].L = 0; nb[2].L = 0; nb[3].L = 0; // инициализировали значение расстояния
	for (auto pnt:net)
		{
		m = (pnt.t-p.t)*(pnt.t-p.t) + (pnt.x-p.x)*(pnt.x-p.x);
		if ((pnt.t>=p.t)&&(pnt.x>=p.x))
			{if (nb[0].m>m)
				{nb[0].t = pnt.t; nb[0].x = pnt.x; nb[0].R = pnt.R; nb[0].L = pnt.L; nb[0].m=m;}}
		else if ((pnt.t<=p.t)&&(pnt.x>=p.x))
			{if (nb[1].m>m)
				{nb[1].t = pnt.t; nb[1].x = pnt.x; nb[1].R = pnt.R; nb[1].L = pnt.L; nb[1].m=m;}}
		else if ((pnt.t<=p.t)&&(pnt.x<=p.x))
			{if (nb[2].m>m)
				{nb[2].t = pnt.t; nb[2].x = pnt.x; nb[2].R = pnt.R; nb[2].L = pnt.L; nb[2].m=m;}}
		else if ((pnt.t>=p.t)&&(pnt.x<=p.x))
			{if (nb[3].m>m)
				{nb[3].t = pnt.t; nb[3].x = pnt.x; nb[3].R = pnt.R; nb[3].L = pnt.L; nb[3].m=m;}}
		}
	return nb;
	}

//=== Вывод рассчитанной сетки в файл ===
void outputNet(const vector<prn2> net, const char* fName, int flag = 0)
	{
	ofstream outN2;
	outN2.open(fName);
	outN2 << fixed << setprecision(maxPrecision);
	if (fName == "netTest.csv")
		{
		outN2 << "s;c;t;x;R;L" << endl;
		for (auto pnt:net)
			{
			outN2 << pnt.s << ";" << pnt.c << ";" << pnt.t << ";" << pnt.x << ";" << pnt.R << ";" << pnt.L << endl;
			}
		}
	if (fName == "net4NN.csv")
		{
		outN2 << "s;c;t;x;R;L;DR_cur;DL_cur" << endl;
		for (auto pnt:net)
			{
			outN2 << pnt.s << ";" << pnt.c << ";" << pnt.t << ";" << pnt.x << ";" <<
					pnt.R << ";" << pnt.L << ";" << pnt.DR_cur << ";" << pnt.DL_cur << endl;
			}
		}
	outN2.close();
	}

//  далее основная процедура ======================================================================
int c4main2(){
	cout << fixed << setprecision(maxPrecision);
	net4NN.clear();
	readNet1("net1.csv", 1);
	readNet1("net4NN.csv", 3);
	int i = 0;
	int s = net4NN.size();
	for (auto &pnt:net4NN)
		{
		if ((pnt.DR_cur >procDR)||(pnt.DL_cur >procDL)) {
			nb1 = findNN(pnt, net1);
			pnt = fourNN(pnt, nb1);
			}
		i++;
		if (i%100 == 0) cout << "считаем точку " << i << " из "<< s << endl;
		}
	outputNet(net4NN, "net4NN.csv");
	cout << "==  Расчет окончен  =================================================================" << endl;
	return 0;
}

//  инициализируем файл net4NN.csv ======================================================================
int c4main21(){
	vector<prn2> netTest, net4NN;
	prn2 p;
	int num = 0;

// открываем net1.csv для заполнения netTest.csv
	ifstream input1("net1.csv");
	if (input1) {cout << "net1.csv is opened" << endl;}
	else {cout << "net1.csv is not opened" << endl;}
	string rd;
	getline(input1, rd);  // прочитали строку с названиями колонок
	while (getline(input1, rd, ';')) {
		p.s = atof(rd.c_str());  // прочитали s
		getline(input1, rd, ';'); // прочитали c
		p.c = atof(rd.c_str());
		getline(input1, rd, ';'); // прочитали t
		p.t = atof(rd.c_str());
		getline(input1, rd, ';'); // прочитали x
		p.x = atof(rd.c_str());
		getline(input1, rd, ';'); // прочитали R
		p.R = atof(rd.c_str());
		getline(input1, rd);      // прочитали L
		p.L = atof(rd.c_str());
		num++;
		if (num % 30000 == 0) netTest.push_back(p); // выгружаем каждую 20ю точку
		}
	outputNet(netTest, "netTest.csv");
	cout << "== Подготовка netTest.csv выполнена ===========================================" << endl;

// открываем net1.csv для заполнения net4NN.csv
	ifstream input2("net1.csv");
	num = 0;
	if (input2) {cout << "net1.csv is opened" << endl;}
	else {cout << "net1.csv is not opened" << endl;}
	getline(input2, rd);  // прочитали строку с названиями колонок
	while (getline(input2, rd, ';')) {
		p.s = atof(rd.c_str());  // прочитали s
		getline(input2, rd, ';'); // прочитали c
		p.c = atof(rd.c_str());
		getline(input2, rd, ';'); // прочитали t
		p.t = atof(rd.c_str());
		getline(input2, rd, ';'); // прочитали x
		p.x = atof(rd.c_str());
		getline(input2, rd, ';'); // прочитали R
		p.R = 100;
		getline(input2, rd);      // прочитали L
		p.L = 100;
		p.DR_cur = 100;
		p.DL_cur = 100;
		num++;
		if (num % 30000 == 0) net4NN.push_back(p); // выгружаем каждую 20ю точку
		}
	outputNet(net4NN, "net4NN.csv");
	cout << "== Инициализация net4NN.csv выполнена ===========================================" << endl;
	return 0;
}

// == процедура проверки, что есть еще что считать=============================================
int c4main22(){
	cout << fixed << setprecision(maxPrecision);
	net4NN.clear();
	readNet1("net4NN.csv", 3);
	int i = 0;
	int s = net4NN.size();
	cout << "== Начинаем проверку net4NN.csv на готовность ==========================================" << endl;
	for (auto &pnt:net4NN)
		{
		if ((pnt.DR_cur > procDR)||(pnt.DL_cur > procDL)) {
			cout << "требуется пересчитать точку, s= " << pnt.s << "; c= " << pnt.c <<
					"; t= " << pnt.t << "; r= " << pnt.x 	<< "; dr= " << pnt.DR_cur << "; dl= " << pnt.DL_cur << endl;
			i++;
			}
		}
	cout << "== Закончили проверку net4NN.csv на готовность ==========================================" << endl;
	cout << "требуется пересчитать " << i << " точек" << endl;
	return 0;
}

