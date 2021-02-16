/*
 * mx.cpp
 *
 *  Created on: 3 ма€ 2020 г.
 *      Author: 1
 */
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
using namespace std;

double c0 = 1;
double gamma;
double rho0;
double rho_;
double m;
int n0c3p1;      // число точек на 0-слое, лини€ ј¬, см.стать€ 1
int n1c3p1;      // число точек на ступеньке плотности в момент сжати€, тачка ј, см.стать€ 1
int n2c3p1;      // число точек на плюс-харктеристике области 3, лини€ ≈D, см.стать€ 1
int nrc3p2;      // параметр-разбиение линии ј¬ (см.стать€ 2 - счет в пр€мом направлении
// изменени€ времени при известном законе движени€ сжимающего поршн€ поршн€)
int nu=2;
const char* nameFile;
double ts, tf, rf, rs, rw, maxNum, u_max, c_max, rho_max;
double pi = 3.1415926535;

double pistonEnd;
int numOnPiston;
const int maxLenPiston = 300000;
const int maxLenLayer = 300000;

double calc_rw() {
	double res;
	if (nu == 0) {
//        res = rs + m/rho0;
        res = rf + m/rho_;
	}
	if (nu == 1) {
//        res = pow(rs*rs + m/(pi*rho0), 1.0/2.0);
        res = pow(rf*rf + m/(pi*rho_), 1.0/2.0);
	}
	if (nu == 2) {
//        res = pow(rs*rs*rs + (3*m)/(4*pi*rho0), 1.0/3.0);
        res = pow(rf*rf*rf + (3*m)/(4*pi*rho_), 1.0/3.0);
	}
	return res;
}

// double rw = calc_rw();

//=======================================================================
// наклон положительной характеристики
double Cp(double R,double L) {
    return R*(gamma+1)/4+L*(3-gamma)/4;
}

//=======================================================================
// наклон отрицательной характеристики
double Cm(double R,double L){
    return L*(gamma+1)/4+R*(3-gamma)/4;
}

//=======================================================================
double Fp(double x, double R, double L) {
    return (-nu)*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
double Fm(double x, double R, double L) {
    return nu*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
void newPoint1(double p1[4],
		double p2[4],
		double p3[4]) {
// ѕересечение плюс и минус характеристик в р
// »з р1 выпускаем —-, из р2 —+
// считаем первое приближение
    double a = Cm(p1[2],p1[3]);
    double b = Cp(p2[2],p2[3]);
    p3[0] = (a*p1[0]-b*p2[0]+p2[1]-p1[1])/(a-b);
    p3[1] = p1[1]+(p2[1]-p1[1]+b*(p1[0]-p2[0]))*a/(a-b);
    p3[2] = p2[2]+(p3[0]-p2[0])*Fp(p2[1],p2[2],p2[3]);
    p3[3] = p1[3]+(p3[0]-p1[0])*Fm(p1[1],p1[2],p1[3]);
// уточн€ем полученные значени€
    a = (Cm(p1[2],p1[3])+Cm(p3[2],p3[3]))/2;
    b = (Cp(p2[2],p2[3])+Cp(p3[2],p3[3]))/2;
    p3[0] = (a*p1[0]-b*p2[0]+p2[1]-p1[1])/(a-b);
    p3[1] = p1[1]+(p2[1]-p1[1]+b*(p1[0]-p2[0]))*a/(a-b);
    p3[2] = p2[2]+(p3[0]-p2[0])*(Fp(p2[1],p2[2],p2[3])+Fp(p3[1],p3[2],p3[3]))/2;
    p3[3] = p1[3]+(p3[0]-p1[0])*(Fm(p1[1],p1[2],p1[3])+Fm(p3[1],p3[2],p3[3]))/2;
}

//=======================================================================
void newPoint2(double* p1, double* p2) {
// ѕересечение плюс характеристики с пр€мой r=rw
// »з р1 выпускаем —+, в р2 нова€ точка
// считаем первое приближение
    double a = Cp(p1[2],p1[3]);
    p2[1] = rw;
    p2[0] = p1[0]+(p2[1]-p1[1])/a;
    p2[2] = p1[2]+(p2[0]-p1[0])*Fp(p1[1],p1[2],p1[3]);
    p2[3] = -p2[2];
// уточн€ем полученные значени€
    a = (Cp(p1[2],p1[3])+Cp(p2[2],p2[3]))/2;
    p2[0] = p1[0]+(rw-p1[1])/a;
    p2[2] = p1[2]+(p2[0]-p1[0])*(Fp(p1[1],p1[2],p1[3])+Fp(p2[1],p2[2],p2[3]))/2;
    p2[3] = -p2[2];
}


//=======================================================================
void readPiston(
		const char* fName,
		double p[maxLenPiston][4])
{
	ifstream input(fName);
	if (input) {
		cout << fName << " is opened" << endl;
	} else {
		cout << fName << " is not opened" << endl;
	}
	string rd;
	getline(input, rd);  // прочитали строку с названи€ми колонок
	int i = -1;
	while (getline(input, rd, ';')) {
		i ++;
		getline(input, rd, ';');
		p[i][0] = atof(rd.c_str());
//		getline(input, rd, ';');
		getline(input, rd);
		p[i][1] = atof(rd.c_str());
//		getline(input, rd, ';');
//		p[i][2] = atof(rd.c_str());
//		getline(input, rd);
//		p[i][3] = atof(rd.c_str());
	}
	numOnPiston = i;
	pistonEnd = p[0][0];
	ts = p[i][0];
	rs = p[i][1];
	tf = p[0][0];
	rf = p[0][1];
//	u_max = (p[0][2]+p[0][3])/2;
//	c_max = (p[0][2]-p[0][3])*(gamma-1)/4;
//	rho_max = pow(c_max, 2/(gamma-1));
	maxNum = i;
	rw = calc_rw();
	cout << "прочитали поршень из файла" << endl;

}

//=======================================================================
bool newPoint3_prog2(double p1[4], double p2[4], double p3[4], double p4[4]) {
// в р1, р2 отрезок траектории поршн€
// в р3 точка сетки из которой ищем пересечение с траекторией поршн€. ≈сли пересечение нашли (пересечение может быть на другом отрезке траеткории),
// то в p4 сохран€ем новую точку и функци€ возврщает true, иначе функци€ возвращает false.
    bool flag3 = 0;
    double u=(p3[2]+p3[3])/2 - (gamma-1)*(p3[2]-p3[3])/4;
    double a = p2[1]-p1[1]+u*(p1[0]-p2[0]);
    double l;
    double unew;
    if (a == 0) {
        l = 10;
    }
    else {
    	l = (p3[1]+u*p1[0]-u*p3[0]-p1[1])/a;
    }
    if ((0<=l)and(l<=1)) {
//    if (((0<=l)and(l<=1))or(numOnPiston == 100)) {
        p4[0] = p1[0] + l*(p2[0]-p1[0]);
        p4[1] = p1[1] + l*(p2[1]-p1[1]);
//        p4[2] = p1[2] + l*(p2[2]-p1[2]);
//        p4[3] = p1[3] + l*(p2[3]-p1[3]);
//        unew = (p1[2]+p1[3])/2 + l*((p2[2]+p2[3])/2-(p1[2]+p1[3])/2);
        unew = (p1[1]-p2[1])/(p1[0]-p2[0]);
        p4[3] = p3[3] + (p4[0]-p3[0])*Fm(p3[1],p3[2],p3[3]);
        p4[2] = 2*unew-p4[3];
        flag3 = 1;
    }
    return flag3;
}

// ======================================================================================
// расчет погрешности масс
double dm(double rs, double rho0, double rfnew, double rho_, double rw)
{
	double m0, m_;
    if (nu == 0)
    	{
        m0 = (rw - rs)*rho0;
        m_ = (rw - rfnew)*rho_;
    	}
    else if (nu == 1)
		{
        m0 = 2*pi*(rw*rw - rs*rs)*rho0;
        m_ = 2*pi*(rw*rw - rfnew*rfnew)*rho_;
		}
    else if (nu == 2)
		{
        m0 = (4/3)*pi*(rw*rw*rw - rs*rs*rs)*rho0;
        m_ = (4/3)*pi*(rw*rw*rw - rfnew*rfnew*rfnew)*rho_;
		}
//	return 100*abs(m0-m_)/m;
	return 100*(m0-m_)/m;
	// знак "-" значит масса стала больше (чего не может быть),
	// знак "+" значит масса стала меньше (чего логично тоже не может быть)
}
