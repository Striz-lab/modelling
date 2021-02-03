#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <ctime>
#include <string>

using namespace std;

int s = 1000;              //количество частиц
double dt = 0.01;        //время перерасчета
int S = 80;                //половина длины ребра куба
int N = 500000;            //количество итераций
double epsilon = 1;
double sigma = 1;
double T = 0.2;
double Rmin = 0.72;  //минимально допустимое расстояние для генерации (квадрат, линейно -- 0.83)

double modulation(double x){
    if (x > S)
        x = x - 2 * S;
    if (x < -S)
        x = x + 2 * S;
    return x;
}

double rmin(double *x, double *y, double *z, int NOW){
    double rmin = 100;
    for (int i = 0; i < NOW; i++)
        for (int j = 0; j < i; j++){
            double rx = modulation(x[i] - x[j]);
            double ry = modulation(y[i] - y[j]);
            double rz = modulation(z[i] - z[j]);

            double r = rx * rx + ry * ry + rz * rz;
            
            if (r < rmin)
                rmin = r;
        }
    return rmin;
}


int main(){
    time_t start, end;
    time (&start);
    
    ofstream fout1, fout2, fout3, fout4, fout5;
    fout1.open("Coor.xyz");
    fout2.open("Energy.csv");
    fout3.open("Maxwell.csv");
    fout4.open("Diffuziya.csv");
    fout5.open("Temperature.csv");

    double x[s],     y[s],     z[s],            //координаты
           xpr[s],   ypr[s],   zpr[s],          //координаты в предыдущий момент времени для уменьшения ошибки
           X, Y, Z, //пока ничего лучше не придумал
           vx[s],    vy[s],    vz[s],           //скорости
           ax[s],    ay[s],    az[s],           //ускорения
           xnach[s], ynach[s], znach[s];        //начальные координаты (для диффузии)

    double  Ekin = 0, Ep = 0, dif = 0, a,
            r = 0, rx, ry, rz;
    
    array<double, 1000> V;
    array<int, 100> nV;

    for (int i = 0; i < s; i++){
        x[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        y[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        z[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        double rm = 0;
        while (rm < Rmin){
            x[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            y[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            z[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            
            rm = rmin(x,y,z,i);
        }
        cout << rmin(x,y,z,i) << " " << i << endl;
    }
        
//температура
    for (int i = 0; i<s; i++){
         double V = pow(3*T, 0.5);
        double ranx = ((double)rand() / RAND_MAX)*2*S - S;
         double rany = ((double)rand() / RAND_MAX)*2*S - S;
         double ranz = ((double)rand() / RAND_MAX)*2*S - S;
         double ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
         
        xpr[i] = x[i] - V * dt * ranx / ran;
        ypr[i] = y[i] - V * dt * rany / ran;
        zpr[i] = z[i] - V * dt * ranz / ran;

        ax[i] = 0;
        ay[i] = 0;
        az[i] = 0;
        xnach[i] = x[i];
        ynach[i] = y[i];
        znach[i] = z[i];
    }
    
    for (int ch = 0; ch <= N; ch++){
        for (int i = 0; i < s; i++){
            for (int j = 0; j < i; j++){
                rx = modulation(x[i] - x[j]);
                ry = modulation(y[i] - y[j]);
                rz = modulation(z[i] - z[j]);

                r = rx * rx + ry * ry + rz * rz;

                a = 24*epsilon*(-2*pow((sigma/r), 7) + pow((sigma/r), 4))/(sigma*sigma);
                    
                ax[i] = ax[i] - a * rx;
                ay[i] = ay[i] - a * ry;
                az[i] = az[i] - a * rz;

                ax[j] = ax[j] + a * rx;
                ay[j] = ay[j] + a * ry;
                az[j] = az[j] + a * rz;
                    
                Ep = Ep + 24*epsilon*(pow((sigma/r), 6) - pow((sigma/r), 3));
            }
        }
        
        if (ch % 1000 == 0){
            fout1 << s << "\n" << ch << "\n";
            for (int i = 0; i < s; i++)
                fout1 << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
        }
        
        for (int i = 0; i < s; i++){
            
            X = x[i];
            Y = y[i];
            Z = z[i];
            
            vx[i] = vx[i] + ax[i] * dt;
            vy[i] = vy[i] + ay[i] * dt;
            vz[i] = vz[i] + az[i] * dt;
            
            x[i] = 2 * x[i] - xpr[i] + ax[i] * dt * dt;
            y[i] = 2 * y[i] - ypr[i] + ay[i] * dt * dt;
            z[i] = 2 * z[i] - zpr[i] + az[i] * dt * dt;
            
            xpr[i] = X;
            ypr[i] = Y;
            zpr[i] = Z;
            
            x[i] = modulation(x[i]);
            y[i] = modulation(y[i]);
            z[i] = modulation(z[i]);
            
            ax[i] = 0;
            ay[i] = 0;
            az[i] = 0;
            
            dif = dif + ((x[i] - xnach[i])) * ((x[i] - xnach[i])) + ((y[i] - ynach[i])) * ((y[i] - ynach[i])) + ((z[i] - znach[i])) * ((z[i] - znach[i]));
            Ekin = Ekin + modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i]))/(2 * dt * dt);
        }
        if (ch % 10 == 0){
        fout2 << ch << "," << fixed << Ep << "," << fixed << Ekin << "," << fixed << Ekin + Ep << "\n";
        fout4 << ch << "," << fixed << dif / s << "\n";
        fout5 << ch << "," << fixed << Ekin / (3*s) << "\n";
        }
 
        Ekin = 0;
        Ep = 0;
        dif = 0;
    }
    
    for (int i = 0; i < s; i++){ //Максвелл начало (для последнего шага)
        V[i] = pow((modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i]))/( dt * dt)), 0.5);
        nV[i] = 0;
    }
    sort(V.begin(), V.end());
    double A = - (V[0] - V[s - 1])/s;
    double B = 0;
    for (int j = 0; j < s; j++){
        for (int i = 0; i < s; i++)
            if ((V[i] <= B + 50 * A) and (V[i] >= B - 50 * A))
                nV[j] ++;
        B += A;
    }
    for (int i = 0; i < s; i++)
        fout3 << i * A << "," << nV[i] << "\n"; //Максвелл конец

    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    fout5.close();
    
    time (&end);
    cout << difftime(end, start) << endl;
    return 0;
}
