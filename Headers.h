#ifndef HEADERS_H_INCLUDED
#define HEADERS_H_INCLUDED

#include <iostream>
#include <graphics.h>
#include <locale.h>
#include <math.h>
#include <conio.h>
#define N 4
#define M 4
#define MAXARR 2000
#define MAXLST 2000

using namespace std;

static int   KOD, NWER;
static double* Py;
static double* Px;
static int   IYREB[MAXLST]; /* Макс Y-коорд активных ребер */
static int   IBGIND;        /* Номер след вершины в списке */
static int   IEDG[MAXARR];  /* Y-коорд вершин по возрастан */
static int   INOM[MAXARR];  /* и их номера в исх масс Py   */
static int   IDLSPI;
static double RXREB[MAXLST];
static double RPRIR[MAXLST];
static double RYSL[MAXLST];
double shadow_pyramid[N][M];
const int DX = 5, DY = 5, color = 11;
const double S1 = 1.1;
const double S2 = 0.95;
const double ALPHA = 0.087;
void multyplication(double fig[N][M], double mass[M][M]);
void offset(double fig[N][M], double dx, double dy);
void scale(double fig[N][M], double S);
void rotateZ(double fig[N][M], double angle);
void rotateY(double fig[N][M], double angle);
void rotateX(double fig[N][M], double angle);
double aver(double range[N][M], int a);
void painter(double fig[N][M]);
void FILSTR(int kod, int iy, double ixn, double ixk);
void SORT(int  n, double *iarr);
void V_FP1(int pixel, int kol, double* Px, double* Py);
void Shadow(double fig1[N][M]);

void menu()
{
    cout << "                 МЕНЮ" << endl;
    cout << "-------------------- I - figure ------------------------" << endl;
    cout << "Перемещение - W, A, S, D" << endl;
    cout << "Масштаб - Q, E" << endl;
    cout << "Поворот - 1, 2, 3" << endl;
    cout << "-------------------- II - figure ------------------------" << endl;
    cout << "Перемещение - Y, G, H, J" << endl;
    cout << "Масштаб - T, U" << endl;
    cout << "Поворот - 4, 5, 6" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "Для выбора первой пирамиды нажмите - 9" << endl;
    cout << "Для выбора второй пирамиды нажмите - 0" << endl;
    cout << "8 - Выход после выхода из режима управления" << endl;
    cout << "ESC - Выход из меню управления" << endl;
}

void multyplication(double fig[N][M], double mass[M][M]) {
	double res[N][M] = { 0,0,0,0 };
	for (int k = 0; k < N; k++) {
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < M; j++)
				res[k][i] += fig[k][j] * mass[j][i];
		}
	}
	for (int k = 0; k < N; k++) {
		for (int i = 0; i < M; i++)
			fig[k][i] = res[k][i];
	}
}

void moving(double fig[N][M], double dx, double dy) {
	double dz = 0;
	double DX_DY[M][M] = { {1, 0, 0, 0} ,
							{0, 1, 0, 0},
							{0, 0, 1, 0},
							{dx, dy, dz,1} };
	multyplication(fig, DX_DY);
}

void scale(double fig[N][M], double S) {
	double x = 0, y = 0, z = 0;
	x = aver(fig,0);
	y = aver(fig,1);
	z = aver(fig,2);
	double	Sx_Sy[M][M] = { {S,0,0,0},
			  {0,S,0,0},
			  {0,0,S,0},
			  {x*(1 - S),y*(1 - S),z*(1 - S),1} };
	multyplication(fig, Sx_Sy);
}

void my_figure(double fig1[N][M]) {
	clearviewport();
	setcolor(color);
	line(fig1[0][0], fig1[0][1], fig1[1][0], fig1[1][1]);
	line(fig1[1][0], fig1[1][1], fig1[2][0], fig1[2][1]);
	line(fig1[2][0], fig1[2][1], fig1[0][0], fig1[0][1]);
	setcolor(color + 1);
	line(fig1[0][0], fig1[0][1], fig1[3][0], fig1[3][1]);
	line(fig1[1][0], fig1[1][1], fig1[3][0], fig1[3][1]);
	line(fig1[2][0], fig1[2][1], fig1[3][0], fig1[3][1]);
	setcolor(color + 2);
	line(fig1[0][0], fig1[0][1], fig1[1][0], fig1[1][1]);
	line(fig1[1][0], fig1[1][1], fig1[2][0], fig1[2][1]);
	line(fig1[2][0], fig1[2][1], fig1[3][0], fig1[3][1]);
	painter(fig1);
	Shadow(fig1);
}
void rotateX(double fig[N][M], double angle) {
	double y = 0, z = 0;
	y = aver(fig,1);
	z = aver(fig,2);
	double Angle[M][M] = { {1,0, 0, 0},
			{0 , cos(angle), sin(angle),0},
			{0, -sin(angle), cos(angle), 0},
			{0, y*(1 - cos(angle)) + z * sin(angle), z*(1 - cos(angle)) - y * sin(angle), 1} };
	multyplication(fig, Angle);
}

void rotateY(double fig[N][M], double angle) {
	double x = 0, z = 0;
	x = aver(fig,0);
	z = aver(fig,2);
	double Angle[M][M] = { {cos(angle),0, -sin(angle), 0},
			{0, 1, 0,0},
			{sin(angle), 0, cos(angle), 0},
			{x*(1 - cos(angle)) - z * sin(angle), 0, z*(1 - cos(angle)) + x * sin(angle), 1} };

	multyplication(fig, Angle);
}

void rotateZ(double fig[N][M], double angle) {
	double x = 0, y = 0;
	x=aver(fig, 0);
	y=aver(fig, 1);
	double Angle[M][M] = { {cos(angle), sin(angle), 0, 0},
		   { -sin(angle), cos(angle), 0, 0},
			 {0, 0, 1, 0},
			 {x*(1 - cos(angle)) + y * sin(angle), y*(1 - cos(angle)) - x * sin(angle), 0, 1} };
	multyplication(fig, Angle);
}

double aver(double fig[N][M], int cnt){
	double average=0;
	for (int i = 0; i < N; i++) {
		average += fig[i][cnt];
	}
	return average/N;
}

void FILSTR (int kod, int iy, double ixn, double ixk)
{
   while (ixn <= ixk) putpixel (ixn++, iy, kod);
}

static int  FORSPI (int IYBEG)
{

   int   i,ikledg,intek,intabs,isd;
   int   iyt,ixt,nrebra,inc,inpred,inposl;
   double xt, xc, yt, yc, dy;

   ikledg= 0;
   for (i=IBGIND; i<=NWER; ++i)
      if (IEDG[i] != IYBEG) break; else ++ikledg;


   for (i=1; i<=ikledg; ++i) {
      intek= INOM[IBGIND+i-1];
      intabs= abs (intek);
      xt= Px[intabs];
      yt= Py[intabs];

      if ((inpred= intabs - 1) < 1) inpred= NWER;
      if ((inposl= intabs + 1) > NWER) inposl= 1;

      for (isd=0;  isd<=1; ++isd) {
         if (!isd) nrebra= inc= inpred; else {
            inc= inposl;  nrebra= intabs;
         }
         yc= Py[inc];
         dy= yc - yt;
         if (dy < 0.0) continue;
         xc= Px[inc];
         if (dy != 0.0) goto DYNE0;
            if ((inc= INOM[nrebra]) < 0) continue;
            INOM[nrebra]= -inc;
            iyt= yt;
            inc= xc;
            ixt= xt;
            FILSTR (KOD, iyt, inc, ixt);
            continue;
DYNE0:   ++IDLSPI;
         IYREB[IDLSPI]= yc;
         RXREB[IDLSPI]= xt;
         RPRIR[IDLSPI]= (xc - xt) / dy;
         inc= (!isd) ? inposl : inpred;
         RYSL[IDLSPI]=  Py[inc] - yt;
      }
   }


   if ((i= (IBGIND += ikledg)) > NWER) i= NWER;
   return (IEDG[i]);
}

void SORT(int  n, double *iarr)
{
	int l;
	double k, min;
	for (int i = 0; i < n; ++i) {
		min = iarr[l = i];
		for (int j = i + 1; j < n; ++j)
			if ((k = iarr[j]) < min) {
				l = j;
				min = k;
			}
		if (l != i) {
			iarr[l] = iarr[i];
			iarr[i] = min;
		}
	}
}

void V_FP1 (int pixel, int kol, double* px, double* py)
{
int  i,j,k,l;
int  iytek;
int  iymin;
int  iybeg;
int  iymak;
int  iysled;

int  newysl;
int  ixmin;
int  ixtek;
int  irabx[MAXLST];

   KOD= pixel;
   NWER= kol;
   Px = px;
   Py = py;

   for (i= 1; i<=NWER; ++i) {IEDG[i]= py[i];  INOM[i]= i; }

   for (i= 1; i<=NWER; ++i) {
      iymin= IEDG[i];
      k= 0;
      for (j=i+1; j<=NWER; ++j)
         if ((l= IEDG[j]) < iymin) {iymin= l; k= j; }
      if (k) {
         IEDG[k]= IEDG[i]; IEDG[i]= iymin;
         iymin= INOM[k];
         INOM[k]= INOM[i]; INOM[i]= iymin;
      }
   }

   IDLSPI= 0;
   IBGIND= 1;
   iybeg= IEDG[1];
   iymak= IEDG[NWER];


   iysled= FORSPI (iybeg);
   if (!IDLSPI) goto KOHGFA;


ZALIWKA:

   for (iytek=iybeg; iytek<=iysled; ++iytek) {
      if (iytek == iysled) {
         newysl= FORSPI (iytek);
         if (!IDLSPI) goto KOHGFA;
      }

      l= 0;
      for (i=1; i<=IDLSPI; ++i)
         if (RYSL[i] > 0.0) irabx[++l]= RXREB[i];
         else RYSL[i]= 1.0;

      for (i=1;  i<=l; ++i) {
         ixmin= irabx[i];
         k= 0;
         for (j=i+1;  j<=l; ++j) {
            ixtek= irabx[j];
            if (ixtek < ixmin) {k= j; ixmin= ixtek; }
         }
         if (k) {irabx[k]= irabx[i];  irabx[i]= ixmin; }
      }

      for (j=1;  j<=l-1;  j+= 2)
         FILSTR (KOD,iytek,irabx[j],irabx[j+1]);

      for (j=1;  j<=IDLSPI; ++j)
         RXREB[j]= RXREB[j] + RPRIR[j];
   }

   if (iysled == iymak) goto KOHGFA;

   i= 0;
M1:++i;
M2:if (i > IDLSPI) goto WYBROSILI;
      if (IYREB[i] != iysled) goto M1;
         --IDLSPI;
         for (j=i;  j<=IDLSPI; ++j) {
            IYREB[j]= IYREB[k= j+1];
            RXREB[j]= RXREB[k];
            RPRIR[j]= RPRIR[k];
         }
         goto M2;
WYBROSILI:
   iybeg= iysled + 1;
   iysled= newysl;
   goto ZALIWKA;

KOHGFA:;
}

void painter(double fig[N][M]){
	double Pol[4] = {}, Pol1[5] = {}, min=0;
	int ntek, inom[6] = {};
//среднее значение глубины
	Pol[0] = (fig[0][2] + fig[1][2] + fig[3][2])/4.0;
	Pol[1] = (fig[0][2] + fig[2][2] + fig[3][2])/4.0;
	Pol[2] = (fig[1][2] + fig[2][2] + fig[3][2])/4.0;
	Pol[3] = (fig[0][2] + fig[1][2] + fig[2][2])/4.0;

	for (int i = 1; i <= 4; ++i) {
		Pol1[i] = Pol[i - 1];
		inom[i] = i;
	}

	for (int i = 1; i <= 4; ++i) {
		min = Pol1[i];
		ntek = i;
		for (int j = i + 1; j <= 4; ++j)
			if (Pol1[j] <= min) {
				min = Pol1[j];
				ntek = j;
			}
		if (ntek != i) {
			Pol1[ntek] = Pol1[i];
			Pol1[i] = min;
			min = inom[ntek];
			inom[ntek] = inom[i];
			inom[i] = min;
		}
	}
	int num;

	for (int i = 0; i < 4; i++) {
		num = inom[i];
		switch (num)
		{
		case 0:
		    Px = new double[3 + 1];
		    Py = new double[3 + 1];
			Px[1] = fig[0][0];
			Py[1] = fig[0][1];
			Px[2] = fig[1][0];
			Py[2] = fig[1][1];
			Px[3] = fig[3][0];
			Py[3] = fig[3][1];
			V_FP1(11, 3, Px, Py);
			delete[] Px;
			delete[] Py;
			break;
		case 1:
		    Px = new double[3 + 1];
		    Py = new double[3 + 1];
			Px[1] = fig[0][0];
			Py[1] = fig[0][1];
			Px[2] = fig[3][0];
			Py[2] = fig[3][1];
			Px[3] = fig[2][0];
			Py[3] = fig[2][1];

			V_FP1(12, 3, Px, Py);
			delete[] Px;
			delete[] Py;
			break;
		case 2:
		    Px = new double[3 + 1];
		    Py = new double[3 + 1];
			Px[1] = fig[1][0];
			Py[1] = fig[1][1];
			Px[2] = fig[3][0];
			Py[2] = fig[3][1];
			Px[3] = fig[2][0];
			Py[3] = fig[2][1];

			V_FP1(13, 3, Px, Py);
			delete[] Px;
			delete[] Py;
			break;
		case 3:
		    Px = new double[3 + 1];
		    Py = new double[3 + 1];
			Px[1] = fig[0][0];
			Py[1] = fig[0][1];
			Px[2] = fig[1][0];
			Py[2] = fig[1][1];
			Px[3] = fig[2][0];
			Py[3] = fig[2][1];

			V_FP1(14, 3, Px, Py);
			delete[] Px;
			delete[] Py;
			break;
		}

	}

}

void Shadow(double fig[M][N]) {

    double height[3] = { 0.0, 0.0, 100.0 };

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            shadow_pyramid[i][j] = fig[i][j];
        }
    }

    for (int i = 0; i < M; i++) {
        shadow_pyramid[i][0] = shadow_pyramid[i][0] + (0.50f * 300.0 - shadow_pyramid[i][1]) *
        (shadow_pyramid[i][0] - height[0]) / (shadow_pyramid[i][1] - height[1]);
        shadow_pyramid[i][2] = fig[i][2] + (0.50f * 300.0 - shadow_pyramid[i][1]) *
        (shadow_pyramid[i][2] - height[2]) / (shadow_pyramid[i][1] - height[1]);
        shadow_pyramid[i][1] = 0.90f * 300.0 - shadow_pyramid[i][2] * 0.8f;
    }

    line(shadow_pyramid[0][0],shadow_pyramid[0][1],shadow_pyramid[1][0],shadow_pyramid[1][1]);
    line(shadow_pyramid[1][0],shadow_pyramid[1][1],shadow_pyramid[2][0],shadow_pyramid[2][1]);
    line(shadow_pyramid[2][0],shadow_pyramid[2][1],shadow_pyramid[0][0],shadow_pyramid[0][1]);
    line(shadow_pyramid[0][0],shadow_pyramid[0][1],shadow_pyramid[3][0],shadow_pyramid[3][1]);
    line(shadow_pyramid[1][0],shadow_pyramid[1][1],shadow_pyramid[3][0],shadow_pyramid[3][1]);
    line(shadow_pyramid[2][0],shadow_pyramid[2][1],shadow_pyramid[3][0],shadow_pyramid[3][1]);

    Px = new double[3 + 1];
    Py = new double[3 + 1];
    Px[1] = shadow_pyramid[0][0];
    Py[1] = shadow_pyramid[0][1];
    Px[2] = shadow_pyramid[1][0];
    Py[2] = shadow_pyramid[1][1];
    Px[3] = shadow_pyramid[3][0];
    Py[3] = shadow_pyramid[3][1];
    V_FP1(7, 3, Px,Py);
    delete[] Px;
    delete[] Py;
    Px = new double[3 + 1];
    Py = new double[3 + 1];
    Px[1] = shadow_pyramid[0][0];
    Py[1] = shadow_pyramid[0][1];
    Px[2] = shadow_pyramid[3][0];
    Py[2] = shadow_pyramid[3][1];
    Px[3] = shadow_pyramid[2][0];
    Py[3] = shadow_pyramid[2][1];
    V_FP1(7, 3, Px,Py);
    delete[] Px;
    delete[] Py;
    Px = new double[3 + 1];
    Py = new double[3 + 1];
    Px[1] = shadow_pyramid[1][0];
    Py[1] = shadow_pyramid[1][1];
    Px[2] = shadow_pyramid[3][0];
    Py[2] = shadow_pyramid[3][1];
    Px[3] = shadow_pyramid[2][0];
    Py[3] = shadow_pyramid[2][1];
    V_FP1(7, 3, Px, Py);
    delete[] Px;
    delete[] Py;
    Px = new double[3 + 1];
    Py = new double[3 + 1];
    Px[1] = shadow_pyramid[0][0];
    Py[1] = shadow_pyramid[0][1];
    Px[2] = shadow_pyramid[1][0];
    Py[2] = shadow_pyramid[1][1];
    Px[3] = shadow_pyramid[2][0];
    Py[3] = shadow_pyramid[2][1];
    V_FP1(7, 3, Px, Py);
    delete[] Px;
    delete[] Py;
}

#endif // HEADERS_H_INCLUDED

