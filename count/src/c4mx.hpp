/*
 * mx.h
 *
 *  Created on: 3 ��� 2020 �.
 *      Author: 1
 */

#ifndef C3MX_HPP_
#define C3MX_HPP_

extern double rho0;
extern int n0c3p1;      // ����� ����� �� 0-����, ����� ��, ��.������ 1
extern int n1c3p1;      // ����� ����� �� ��������� ��������� � ������ ������, ����� �, ��.������ 1
extern int n2c3p1;      // ����� ����� �� ����-������������� ������� 3, ����� �D, ��.������ 1
extern int nrc3p2;      // ��������-��������� ����� �� (��.������ 2 - ���� � ������ �����������
// ��������� ������� ��� ��������� ������ �������� ���������� ������ ������)
extern double c0;
extern double gamma;
extern double rho_;
extern double m;
extern double pi;
extern double ts, tf, rf, rs, rw, maxNum, u_max, c_max, rho_max;

extern int nu;
extern const char* nameFile;
extern const int maxLenPiston;
extern const int maxLenLayer;
extern double pistonEnd;
extern int numOnPiston;

extern double Cp(double R,double L);
extern double Cm(double R,double L);
extern double Fp(double x, double R, double L);
extern double Fm(double x, double R, double L);
extern void newPoint1(double p1[4], double p2[4], double p3[4]);
extern void newPoint2(double* p1, double* p2);
extern bool newPoint3_prog2(double p1[4], double p2[4], double p3[4],  double p4[4]);
extern double calc_rw();
extern void readPiston(const char* fName, double p[200000][4]);
double dm(double rs, double rho0, double rfnew, double rho_, double rw);

#endif /* C3MX_HPP_ */
