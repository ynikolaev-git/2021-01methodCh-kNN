//============================================================================
// Name        : c4prog1.cpp
// Author      : Nikolaev Yuri
// Version     : 01
// Copyright   :
// Description : ������� �� ������ ��������� ����� � �������� ����������� �������.
// ��������� ��� �������� ��������� �����, �������� �������� ������� ���������� ��� �������� ��������� ��������� ������
// � ������� ������� ������� � ���������� ������.
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "c4mx.hpp"
using namespace std;

const int maxSizeOfLayer = 1000000;
const int maxPrecision = 24;
double r1;
int cC = 0;                           // ����� ����� ������� ����� ����� � ����������� ���.�����
bool msg = 1;
double r_;
double c_;
double t012;
double dt0;
ofstream output;
//int kratnostP = 100, kratnostL = 100;
int kratnostP = 1, kratnostL = 1;

struct pointOfNet {
	int s,c;
	double t,x,R,L;
};

vector<pointOfNet> c3p1Layer1(maxSizeOfLayer);
vector<pointOfNet> c3p1Layer2(maxSizeOfLayer);
vector<pointOfNet> piston;
vector<pointOfNet> net;
pointOfNet p1, p2, c3p1p, pp;

//=======================================================================
// ������ ���������� r ������������ ������
double c3prog1_calc_rw() {
	double res;
    if (nu == 0) res = r1 + m/rho0;
    else if (nu == 1) res = pow(r1*r1 + m/(pi*rho0) , 1.0/2.0);
    else if (nu == 2) res = pow(r1*r1*r1 + (3*m)/(4*pi*rho0) , 1.0/3.0);
    return res;
}

//=======================================================================
// ������ ���������� r_
double calc_r_() {
	double res;
    if (nu == 0) res = rw - m/rho_;
    else if (nu == 1) res = pow(rw*rw - m/(rho_*pi) , 1.0/2.0);
    else if (nu == 2) res = pow(rw*rw*rw - (3*m)/(4*rho_*pi) , 1.0/3.0);
    return res;
}

//=======================================================================
pointOfNet c3p1newPoint1(pointOfNet p1, pointOfNet p2) {
// ����������� ���� � ����� ������������� � �
// �� �1 ��������� �-, �� �2 �+
// ������� ������ �����������
	pointOfNet p;
	double a = Cm(p1.R,p1.L);
    double b = Cp(p2.R,p2.L);
    p.t = (a*p1.t-b*p2.t+p2.x-p1.x)/(a-b);
    p.x = p1.x+(p2.x-p1.x+b*(p1.t-p2.t))*a/(a-b);
    p.R = p2.R+(p.t-p2.t)*Fp(p2.x,p2.R,p2.L);
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
// �������� ���������� ��������
    a = (Cm(p1.R,p1.L)+Cm(p.R,p.L))/2;
    b = (Cp(p2.R,p2.L)+Cp(p.R,p.L))/2;
    p.t = (a*p1.t-b*p2.t+p2.x-p1.x)/(a-b);
    p.x = p1.x+(p2.x-p1.x+b*(p1.t-p2.t))*a/(a-b);
    p.R = p2.R+(p.t-p2.t)*(Fp(p2.x,p2.R,p2.L)+Fp(p.x,p.R,p.L))/2;
    p.L = p1.L+(p.t-p1.t)*(Fm(p1.x,p1.R,p1.L)+Fm(p.x,p.R,p.L))/2;
    return p;
}
//=======================================================================
pointOfNet c3p1newPoint2(const pointOfNet p1) {
// ����������� ����� �������������� � ������ r=rw
// ������� ������ �����������
	pointOfNet p;
    double a = Cm(p1.R,p1.L);
    p.x = rw;
    p.t = p1.t+(p.x-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
    p.R = -p.L;
// �������� ���������� ��������
    a = (Cm(p1.R,p1.L)+Cm(p.R,p.L))/2;
    p.t = p1.t+(rw-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*(Fm(p1.x,p1.R,p1.L)+Fm(p.x,p.R,p.L))/2;
    p.R = -p.L;
    return p;
}

//==========================================================================
bool c3p1cross(pointOfNet p1, pointOfNet p2, pointOfNet& p) {
// �1, �2 ���� ������������������ ����� ����� ������� �������� ��������������.
// � � ����� ���������� �� ������� ���������� �������� ������. ���� flag=true,
// �� � � ����� ����� ����������, ����� � �� ����������}
    bool flagI = 0;
    double l;
    double u=(p.R+p.L)/2;
    double a = p2.x-p1.x+u*(p1.t-p2.t);
    if (a!=0) l= (p.x+u*p1.t-u*p.t-p1.x)/a;
    else l = 10;
    if ((0<=l)and(l<=1)) {
        p.t = p1.t + l*(p2.t-p1.t);
        p.x = p1.x + l*(p2.x-p1.x);
        p.L = p1.L + l*(p2.L-p1.L);
        p.R = p1.R + l*(p2.R-p1.R);
//          print('(0<=l<=1)')
        flagI = 1;
    }
    return flagI;
}

//=======================================================================
// ������ ����������� ����
double c3p1dm() {
    double rp_ = piston[0].x;
    double rp0 = piston[piston.size()-1].x;
    double m0, m_;
    if (nu == 0) {
        m0 = (rw - rp0)*rho0;
        m_ = (rw - rp_)*rho_;
    	}
    else if (nu == 1) {
        m0 = 2*pi*(rw*rw - rp0*rp0)*rho0;
        m_ = 2*pi*(rw*rw - rp_*rp_)*rho_;
    	}
    else if (nu == 2) {
        m0 = (4.0/3.0)*pi*(rw*rw*rw - rp0*rp0*rp0)*rho0;
        m_ = (4.0/3.0)*pi*(rw*rw*rw - rp_*rp_*rp_)*rho_;
    	}
    return 100 * abs(m0-m_)/m;
}
//=======================================================================
// ������� ������ ������
double c3p1Ap() {
    double res0, res1, res = 0;
    int i, k = piston.size();
    if (nu == 0) {
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= (res1+res0)*(piston[i+1].x-piston[i].x)/(2*gamma);
        	}
    	}
    else if (nu == 1) {
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= pi*(res1+res0)*(piston[i+1].x*piston[i+1].x-piston[i].x*piston[i].x)/(2*gamma);
        	}
    	}
    else if (nu == 2)
    	{
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= pi*(4.0/3.0)*(res1+res0)*(piston[i+1].x*piston[i+1].x*piston[i+1].x
                           -piston[i].x*piston[i].x*piston[i].x)/(2*gamma);
        	}
    	}
    return res;
}

//=======================================================================
// ��������� � ������ �����
void c3p1printStart() {
    cout << "     ==== �������� ���� ====" << endl;
//    print('������� ����� � ' + str(l-1))
    cout << "���������� �������� ���� = " << gamma << endl;
    cout << "���������� ��������� ������� = " << nu << endl;
    cout << "���������� �� ������ ��������� �� ����������� ������ = " << rw << endl;
    cout << "���������� �� ������ ��������� �� ������ � ������ ������ = " << r_ << endl;
    cout << "��������� ������� ���� = " << rho_ << endl;
    cout << "�������� ����� � ������ ���� = " << c_ << endl;
}

//=======================================================================
// ��������� � ����� �����
void c3p1printEnd() {
    cout << "                     ==== ��������� ���� ====" << endl;
    cout << "����� ������ (���������������� ����������)= " << piston[piston.size()-1].x << endl;
    cout << "����� ������ (�����)= " << piston[piston.size()-1].t << endl;
    cout << "����������� ����= " << c3p1dm() << endl;
    cout << "����� ������ (���������������� ����������)= " << piston[0].x << endl;
    cout << "����� ������ (�����)= " << piston[0].t << endl;
    cout << "������ ������ = " << c3p1Ap() << endl;
//    print('=== ������������ ����� ��������, ������= ', int(time.time() - start_time))
}

//=======================================================================
// // ��������� ���� ����� � ���� f
void c4toNet(vector<pointOfNet>& l, int start, int end)
	{
	pointOfNet l1;
//	cout << l[end].t << "; " << l[end].x << ";" << endl;
//	if (not((l[1].s > 7550)&&(l[1].s < 7700))) return;
	if (kratnostL == 1)
		{
		kratnostL = 0;
		kratnostL++;  // ����������������, ���� ����� ��������� ���� ������
		kratnostP = 1;
		for (int k=start;k<end; k+=kratnostP)  // ��������� ����� � ����
    		{
			l1 = c3p1Layer1[k];
//			if (k > 1000) kratnostP = 2;      // 1000 ������� ����� ��������, ���������
			// ���������� ����� � ��������, ����� ��������� ����� �����
			output << l[k].s << ";" << l[k].c << ";" << l[k].t << ";" << l[k].x << ";" <<
					l[k].R << ";" << l[k].L << endl;
//			cout << l1.s << ";" << l1.c << ";" << l1.t << ";" << l1.x << ";" <<
//					l1.R << ";" << l1.L << endl;
    		}
		}
	else kratnostL +=1;
	}

//=======================================================================
// ������ �����
void mainFunction() {
	int i, j;
	double u;
    int numL = 0;    // ����� ���� (������������� ��������������)
    int numOnL = 0;  // ����� �� ���� (������ � ����� � ������������ - ������ �� ������ � ������ ������)

//##################### ������ ������� ���� (�� �^- ��������������)
    int k = 1;
    double t, x, R, L;
    t = tf;
    R = c_*2 / (gamma-1);
    L = -R;
    while (t > t012) {
    	x = r_-c_*(t-tf);
        c3p1Layer1[numOnL] = {0, numOnL, t, x, R, L};
        if (k < 500) t -= dt0*0.01;
        else t -= dt0;
        k += 1;
        numOnL += 1;
    	}
    x = r_-c_*(t012-tf);
    c3p1Layer1[numOnL] = {0, numOnL, t012, x, R, L};
    cout << "�������� 0-����, ������������ ����� ����� �� ����= " << numOnL << endl;
    c4toNet(c3p1Layer1, 1, numOnL); // ��������� ����� � ����

//############ ����� �������� ����� � ������� 1, 2
    double c = c_;
    cout << "c_ = " << c_ << endl;
    double dc1 = c_/n1c3p1;
    for (j=0; j<n1c3p1; j++) {      // j - ���������� ����
        c = c - dc1;
        numL += 1;
        u = 2*(c_-c)/(gamma-1);
        c3p1p = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        c3p1Layer2[0] = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        for (i=0; i<numOnL + j; i++) {      // i - ���������� ����� �� �����
            p1 = c3p1Layer2[i];
            p2 = c3p1Layer1[i+1];
            c3p1p = c3p1newPoint1(p1, p2);
            c3p1p.s = j+1;
            c3p1p.c = i+1;
            c3p1Layer2[i+1] = c3p1p;
            cC += 1;
        	}
        // ������ ����������� � ������  rw
        c3p1p = c3p1newPoint2(c3p1p);
        c3p1p.s = j+1;
        c3p1p.c = i+1;
        c3p1Layer2[i+1] = c3p1p;
        cC += 1;
        c3p1Layer1 = c3p1Layer2;
        cout << "��������� ���� ����� " << j+1 << endl;
        c4toNet(c3p1Layer1, 1, i+1); // ��������� ����� � ����
        rho0 = (c3p1p.R - c3p1p.L)*(gamma-1)/4;   // ��� �0
        if (rho0 <=1) {
            cout << "�������� ��������� c= " << rho0 << endl;
            break;
        	}
    	}
    numL = j + 1;    // ��������� ����������� ����
    numOnL = i + 1;    // ��������� ����� �� ������������� ����
    cout << "��������� ������� 1: ����� ����= " << numL << ", ������������ ����� ����� �� ����= " << numOnL << endl;
    cout << "����� �-����� ��������������, ������������ ������� 1 � 2 = " << numL << endl;


// ����� ������ ����� � ������� 3
    int p_numL = numL;
    int p_numOnL = 0;
    pp = c3p1Layer2[0];
    piston.push_back(pp);
    double t0 = c3p1p.t;
    double dt2 = (rw-r1)/n2c3p1;
    int cnt = 0;
    i = 0;
    while (p_numOnL < numOnL) { // ���������� ����, ���������� � �������� �������� �����
        c3p1Layer2[numOnL] = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
//        c3p1Layer2[numOnL] = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2*rho0  / (gamma-1), 2*rho0  / (1-gamma)};
        j = numOnL;
        c3p1p = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
//        c3p1p = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2*rho0 / (gamma-1), 2*rho0 / (1-gamma)};
        while (j > cnt) {   // ������ ���� � ������� ����� �������
            p1 = c3p1p;
            p2 = c3p1Layer1[j-1];
            c3p1p = c3p1newPoint1(p1, p2);
            c3p1p.s = numL+1+i;
            c3p1p.c = j-1;
            c3p1Layer2[j-1] = c3p1p;
            cC += 1;
            j -= 1;
    		}
//       �������� ����������� ������ ������� ������ ������ �����
        bool flagDown = 1;
        bool flag;
        while ((flagDown==1)and(p_numOnL < numOnL)) {
            p1 = c3p1Layer2[p_numOnL+1];
            p2 = c3p1Layer1[p_numOnL+1];
            flag = 0;
            flag = c3p1cross(p1, p2, pp);
            if (flag==1) {
                p_numOnL += 1;
                pp.s = p1.s;
                pp.c = p_numOnL;
                piston.push_back(pp);
            	}
            else flagDown = 0;
        	}
//       �������� ����������� ����� ������� �������� ������
        if (p_numOnL == numOnL) {
            cout << "��������� ������ ���������� ������" << endl;
            cout << "������ ������ = " << piston.size() << endl;
            break;
        	}
        p1 = c3p1Layer2[p_numOnL];
        p2 = c3p1Layer2[p_numOnL+1];
        flag = 0;
        flag = c3p1cross(p1,  p2, pp);
        if (flag==1) {
            p_numL += 1;
            piston.push_back(pp);
            cnt = p_numOnL;
        	}
        else {
            cout << "��������� ����������, �������� ��������� �����" << endl;
            //print (p_numL, p_numOnL, t, x, R, L)
            break;
        	}
        c3p1Layer1 = c3p1Layer2;
        cout << "��������� ���� ����� " << numL+1+i << " �� ���������� ������ " << numOnL - j << " �����" << endl;
        c4toNet(c3p1Layer1,j+1,numOnL); // ��������� ����� � ����
        i += 1;
    	}
    numL += i;

    cout << "��������� ������� 2: ����� ����= " << numL << ", ������������ ����� ����� �� ����= " << numOnL << endl;
    cout << "����� �����= " << cC << endl;
    // ����� ���������� � ����
}


//############################################
//    ### ������ ���� ��������� ����� ###
//############################################
int c4main1() {
	tf = 1;
	r1 = 1;
	rw = c3prog1_calc_rw();
	r_ = calc_r_();
	c_ = pow(rho_, (gamma-1)/2.0);
	t012 = tf - (rw - r_) / c_;
	dt0 = (tf - t012) / n0c3p1;

	// !!! ����� �����
	if (msg==1) c3p1printStart();
	nameFile = "net1.csv";
	output.open(nameFile);
	output << fixed << setprecision(maxPrecision);
	cout << fixed << setprecision(maxPrecision);
	output << "s;c;t;x;R;L;" << endl;
	mainFunction();
	output.close();
//	output.open("piston.csv");
	output.open("piston.csv");
	output << fixed << setprecision(maxPrecision);
	output << "numOfPoint;t;x;" << endl;
    for (int i=0;i<piston.size();++i)
    	{
    	piston[i].t -= piston[piston.size()-1].t;
    	output << i << ";" << piston[i].t << ";" << piston[i].x << endl;
    	}
	output.close();
    cout << "��������� ����� ���������� � t=0. ��������� ���������� ������ � ����" << endl;
    if (msg==1) c3p1printEnd();
// !!! ���������� �����

    double u_pred, u, u_next, c_pred, c, c_next, t_s, t_f;
    t_s = piston[piston.size()-1].t;
    t_f = piston[0].t;
    for (int i=piston.size()-1;i>0;i--)
    	{
    	u_pred = (piston[i-1].R+piston[i-1].L)/2;
    	u = (piston[i].R+piston[i].L)/2;
    	u_next = (piston[i+1].R+piston[i+1].L)/2;
    	c_pred = (gamma-1)*(piston[i-1].R-piston[i-1].L)/4;
    	c = (gamma-1)*(piston[i].R-piston[i].L)/4;
    	c_next = (gamma-1)*(piston[i+1].R-piston[i+1].L)/4;
    	// ���� �������
    	if ((u_pred>u)&(u_next>u))  // ����� ������� ��������
    		cout << "����� min ��������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred>c)&(c_next>c))  // ����� ������� ���������
    		cout << "����� min ���������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	// ���� ��������
    	if ((u_pred<u)&(u_next<u))  // ����� �������� ��������
    		cout << "����� max ��������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred<c)&(c_next<c))  // ����� �������� ���������
    		cout << "����� max ���������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
    		", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	}
    cout << "������������ ����� ����� �� ���������� ������ = " << piston.size()-1 << endl;

return 0;
}
