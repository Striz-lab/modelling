#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <ctime>
#include <string>

using namespace std;

    int s = 1000;              //количество частиц
    double dt = 0.0001;        //время перерасчета
    int S = 50;                //половина длины ребра куба
    int N = 100000;            //количество итераций

struct pairxv{
    double a, b;
};

struct pairxv retraction(double x, double v){ //отражение от стенок
    if (x > S){
        v = - v;
        x = S - 1;
    }
    else{
        if (x < - S){
            v = - v;
            x = - S + 1;
        }
    }
    
    if (abs(v) > pow(10, 2))
        v = (rand() % 100) - 50;
    struct pairxv tmp = {x, v};
    return tmp;
}

struct pairxv teleportation(double x, double v){ //телепортация частичек
    if (x > S)
        x = - S;
    if (x < -S)
        x = S;
    
    if (abs(v) > pow(10, 2))
        v = (rand() % 100) - 50;
    struct pairxv tmp = {x, v};
    return tmp;
}


int main(){
    time_t start, end;
    time (&start);
    
    ofstream fout1, fout2, fout3, fout4;
    fout1.open("Coor.xyz");
    fout2.open("Energy.csv");
    fout3.open("Maxwell.csv");
    fout4.open("Diffuziya.csv");

    double x[s],  y[s],  z[s],            //координаты
           vx[s], vy[s], vz[s],           //скорости
           ax[s], ay[s], az[s],           //ускорения
           xnach[s], ynach[s], znach[s];  //начальные координаты (для диффузии)

    double  Ekin = 0, Ep = 0, dif = 0, a,
            r, rx, ry, rz;
    
    array<double, 1000> V;
    array<int, 100> nV;
    
    for (int i = 0; i < s; i++){
        x[i] = (i % 10) * 10 - 45;
        y[i] = ((i / 10) % 10) * 10 - 45;
        z[i] = ((i / 100) % 10) * 10 - 45;
        vx[i] = (rand() % 10000) - 5000;
        vy[i] = (rand() % 10000) - 5000;
        vz[i] = (rand() % 10000) - 5000;
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
                    rx = x[i] - x[j];
                    ry = y[i] - y[j];
                    rz = z[i] - z[j];

                    r = rx * rx + ry * ry + rz * rz;
                
                if (r != 0){
                    a = 24*(-2*pow(r, -7) + pow(r, -4));
                    
                    ax[i] = ax[i] - a * rx;
                    ay[i] = ay[i] - a * ry;
                    az[i] = az[i] - a * rz;

                    ax[j] = ax[j] + a * rx;
                    ay[j] = ay[j] + a * ry;
                    az[j] = az[j] + a * rz;
                    
                    Ep = Ep + 24*(pow(r, -6) - pow(r, -3));
                }
            }
        }
        
        if (ch % 100 == 0){
            fout1 << s << "\n" << ch << "\n";
            for (int i = 0; i < s; i++)
                fout1 << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
        }
        
        for (int i = 0; i < s; i++){
            vx[i] = vx[i] + ax[i] * dt;
            vy[i] = vy[i] + ay[i] * dt;
            vz[i] = vz[i] + az[i] * dt;
            
            x[i] = x[i] + vx[i] * dt + ax[i] * dt * dt / 2;
            y[i] = y[i] + vy[i] * dt + ay[i] * dt * dt / 2;
            z[i] = z[i] + vz[i] * dt + az[i] * dt * dt / 2;
            
            /*struct pairxv pairx = retraction(x[i], vx[i]);
            x[i] = pairx.a;
            vx[i] = pairx.b;
            
            struct pairxv pairy = retraction(y[i], vy[i]);
            y[i] = pairy.a;
            vy[i] = pairy.b;
            
            struct pairxv pairz = retraction(z[i], vz[i]);
            z[i] = pairz.a;
            vz[i] = pairz.b;*/
            
            struct pairxv pairx = teleportation(x[i], vx[i]);
            x[i] = pairx.a;
            vx[i] = pairx.b;
            
            struct pairxv pairy = teleportation(y[i], vy[i]);
            y[i] = pairy.a;
            vy[i] = pairy.b;
            
            struct pairxv pairz = teleportation(z[i], vz[i]);
            z[i] = pairz.a;
            vz[i] = pairz.b;
            
            ax[i] = 0;
            ay[i] = 0;
            az[i] = 0;
            
            Ekin = Ekin + (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i])/2;
            dif = dif + (x[i] - xnach[i]) * (x[i] - xnach[i]) + (y[i] - ynach[i]) * (y[i] - ynach[i]) + (z[i] - znach[i]) * (z[i] - znach[i]);
            Ekin = Ekin + (x[i] - xnach[i]) * (x[i] - xnach[i]) + (y[i] - ynach[i]) * (y[i] - ynach[i]) + (z[i] - znach[i]) * (z[i] - znach[i])/(dt * dt);
        }
    
        fout2 << ch << "," << fixed << Ep << "," << fixed << Ekin << "," << fixed << Ekin + Ep << "\n";

        fout4 << ch << "," << fixed << dif / s << "\n";
        
        Ekin = 0;
        Ep = 0;
    }
    
    for (int i = 0; i < s; i++){ //Максвелл начало (для последнего шага)
        V[i] = pow((vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]), 0.5);
        nV[i] = 0;
    }
    sort(V.begin(), V.end());
    double A = - (V[0] - V[s - 1])/s;
    double B = 0;
    for (int j = 0; j < s; j++){
        for (int i = 0; i < s; i++)
            if ((V[i] <= B + 10 * A) and (V[i] >= B - 10 * A))
                nV[j] ++;
        B += A;
    }
    for (int i = 0; i < s; i++)
        if ((nV[i] != 0) and (C < 0.7*V[999]) and (nV[i] < 200))
            fout3 << i * A << "," << nV[i] << "\n"; //Максвелл конец

    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    
    time (&end);
    cout << difftime(end, start) << endl;
    return 0;
}
