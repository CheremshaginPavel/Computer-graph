#include "Headers.h"

int main() {
	setlocale(LC_ALL, "Russian");
	initwindow(1000, 600, "", 100, 100);
	double figure[N][M] = {{85, 85, 0, 1},
		 				   {105, 85, 0, 1},
		   				   {95, 105, 0, 1},
		  				   {85, 85, 105, 1}};
    double figure1[N][M] = {{110, 110, 0, 1},
		 				   {130, 110, 0, 1},
		   				   {200, 130, 0, 1},
		  				   {110, 110, 30, 1}};
	double pro[M][M] = { {1, 0, 0, 0},
						 {0, 1, 0, 0},
						 {0, 0, 1, 0},
						 {0, 0, 0, 1} };
	multyplication(figure, pro);
	multyplication(figure1, pro);
	setcolor(4);
	menu();
	do {
        if (GetAsyncKeyState((unsigned short)'8') & 0x8000) break;

        if (GetAsyncKeyState((unsigned short)'9') & 0x8000) {
            while (true) {

                if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) break;

                if (GetAsyncKeyState((unsigned short)'W') & 0x8000)
                    moving(figure, 0, -DY), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'S') & 0x8000)
                    moving(figure, 0, DY), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'A') & 0x8000)
                    moving(figure, -DX, 0), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'D') & 0x8000)
                    moving(figure, DX, 0), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'Q') & 0x8000)
                    scale(figure, S1), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'E') & 0x8000)
                    scale(figure, S2), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'1') & 0x8000)
                    rotateZ(figure, -ALPHA), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'2') & 0x8000)
                    rotateY(figure, ALPHA), my_figure(figure);

                if (GetAsyncKeyState((unsigned short)'3') & 0x8000)
                    rotateX(figure, -ALPHA), my_figure(figure);
            }
        }
        if (GetAsyncKeyState((unsigned short)'0') & 0x8000) {
            while (true) {
                if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) break;

                if (GetAsyncKeyState((unsigned short)'Y') & 0x8000)
                    moving(figure1, 0, -DY), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'H') & 0x8000)
                    moving(figure1, 0, DY), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'G') & 0x8000)
                    moving(figure1, -DX, 0), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'J') & 0x8000)
                    moving(figure1, DX, 0), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'T') & 0x8000)
                    scale(figure1, S1), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'U') & 0x8000)
                    scale(figure1, S2), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'4') & 0x8000)
                    rotateZ(figure1, -ALPHA), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'5') & 0x8000)
                    rotateY(figure1, ALPHA), my_figure(figure1);

                if (GetAsyncKeyState((unsigned short)'6') & 0x8000)
                    rotateX(figure1, -ALPHA), my_figure(figure1);
            }
        }
        clearviewport();
		delay(5);
	} while (true);
	closegraph();
	return 0;
}

