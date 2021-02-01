#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <ctime>
#include <string>

using namespace std;

int s = 1000;              //количество частиц
double dt = 0.0001;        //время перерасчета
int S = 9;                //половина длины ребра куба
int N = 10000;            //количество итераций
double E0 = 0;
double Ep0 = 0;
double epsilon = 1;
double sigma = 1;
double T = 300;

double modulation(double x){
    if (x > S)
        x = x - 2 * S;
    if (x < -S)
        x = x + 2 * S;
    return x;
}

struct pairxv{
    double a, b;
};

struct pairxv teleportation(double x, double v, double E){ //телепортация частичек
    x = modulation(x);
    //if (v < 10)
        //v = v*1.1;
    //if (abs(E) < 0.7*E0)
        //v = 1.1 * v;
    //if (abs(E) > 1.3*E0)
        //v = v / 1.1;
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

    double x[s],     y[s],     z[s],            //координаты
           xpr[s],   ypr[s],   zpr[s],          //координаты в предыдущий момент времени для уменьшения ошибки
           X, Y, Z, //пока ничего лучше не придумал
           vx[s],    vy[s],    vz[s],           //скорости
           ax[s],    ay[s],    az[s],           //ускорения
           xnach[s], ynach[s], znach[s];        //начальные координаты (для диффузии)

    double  Ekin = 0, Ep = 0, dif = 0, a,
            r, rx, ry, rz;
    
    array<double, 1000> V;
    array<int, 100> nV;
    
    for (int i = 0; i < s; i++){
        x[i] = ((double)rand() / RAND_MAX)*2*S - S;
        y[i] = ((double)rand() / RAND_MAX)*2*S - S;
        z[i] = ((double)rand() / RAND_MAX)*2*S - S;
        for (int j = 0; j < i; j++){
            rx = modulation(x[i] - x[j]);
            ry = modulation(y[i] - y[j]);
            rz = modulation(z[i] - z[j]);

            r = rx * rx + ry * ry + rz * rz;
            
            if (r < 0.3){
                x[i] = ((double)rand() / RAND_MAX)*2*S - S;
                y[i] = ((double)rand() / RAND_MAX)*2*S - S;
                z[i] = ((double)rand() / RAND_MAX)*2*S - S;
            }
        }
        
//температура
        
        
         
         double V = pow(T, 0.5);
         double ranx = (rand() % 1000) - 500;
         double rany = (rand() % 1000) - 500;
         double ranz = (rand() % 1000) - 500;
         double ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
         if (ran == 0){
             ranx = 1;
             rany = 1;
             ranz = 1;
             ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
         }
         
        xpr[i] = x[i] - V * dt * ranx / ran;
        ypr[i] = y[i] - V * dt * rany / ran;
        zpr[i] = z[i] - V * dt * ranz / ran;
        /*vx[i] = (rand() % 10) - 5;
        vy[i] = (rand() % 10) - 5;
        vz[i] = (rand() % 10) - 5;*/
        ax[i] = 0;
        ay[i] = 0;
        az[i] = 0;
        xnach[i] = x[i];
        ynach[i] = y[i];
        znach[i] = z[i];
        /*xpr[i] = x[i] - vx[i] * dt;
        ypr[i] = y[i] - vy[i] * dt;
        zpr[i] = z[i] - vz[i] * dt;*/
        //E0 += (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i])/2;
        
        for (int j = 0; j < i; j++){
            rx = modulation(x[i] - x[j]);
            ry = modulation(y[i] - y[j]);
            rz = modulation(z[i] - z[j]);

            r = rx * rx + ry * ry + rz * rz;
        
        if ((r != 0) and (r <= 50))
            Ep0 = Ep0 + 24*epsilon*(pow((r/sigma), -6) - pow((r/sigma), -3));
    
        }
        E0 = Ep0 + E0;
    }
    
    for (int ch = 0; ch <= N; ch++){
        for (int i = 0; i < s; i++){
            for (int j = 0; j < i; j++){
                    rx = modulation(x[i] - x[j]);
                    ry = modulation(y[i] - y[j]);
                    rz = modulation(z[i] - z[j]);

                    r = rx * rx + ry * ry + rz * rz;
                
                if (r <= 50){
                    a = 24*epsilon*(-2*pow((r/sigma), -7) + pow((r/sigma), -4))/(sigma*sigma);
                    
                    ax[i] = ax[i] - a * rx;
                    ay[i] = ay[i] - a * ry;
                    az[i] = az[i] - a * rz;

                    ax[j] = ax[j] + a * rx;
                    ay[j] = ay[j] + a * ry;
                    az[j] = az[j] + a * rz;
                    
                    Ep = Ep + 24*epsilon*(pow((r/sigma), -6) - pow((r/sigma), -3));
                }
            }
        }
        
        if (ch % 100 == 0){
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
            
            //x[i] = x[i] + vx[i] * dt + ax[i] * dt * dt / 2;
            //y[i] = y[i] + vy[i] * dt + ay[i] * dt * dt / 2;
            //z[i] = z[i] + vz[i] * dt + az[i] * dt * dt / 2;
            
            x[i] = 2 * x[i] - xpr[i] + ax[i] * dt * dt;
            y[i] = 2 * y[i] - ypr[i] + ay[i] * dt * dt;
            z[i] = 2 * z[i] - zpr[i] + az[i] * dt * dt;
            
            xpr[i] = X;
            ypr[i] = Y;
            zpr[i] = Z;
            
            struct pairxv pairx = teleportation(x[i], vx[i], Ekin + Ep);
            x[i] = pairx.a;
            vx[i] = pairx.b;
            
            struct pairxv pairy = teleportation(y[i], vy[i], Ekin + Ep);
            y[i] = pairy.a;
            vy[i] = pairy.b;
            
            struct pairxv pairz = teleportation(z[i], vz[i], Ekin + Ep);
            z[i] = pairz.a;
            vz[i] = pairz.b;
            
            ax[i] = 0;
            ay[i] = 0;
            az[i] = 0;
            
            //Ekin = Ekin + (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i])/2;
            //dif = dif + (modulation(x[i] - xnach[i])) * modulation((x[i] - xnach[i])) + modulation((y[i] - ynach[i])) * modulation((y[i] - ynach[i])) + modulation((z[i] - znach[i])) * modulation((z[i] - znach[i]));
            dif = dif + ((x[i] - xnach[i])) * ((x[i] - xnach[i])) + ((y[i] - ynach[i])) * ((y[i] - ynach[i])) + ((z[i] - znach[i])) * ((z[i] - znach[i]));
            Ekin = Ekin + modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i]))/(2 * dt * dt);
            T = Ekin/s;
        }
    
        fout2 << ch << "," << fixed << Ep << "," << fixed << Ekin << "," << fixed << Ekin + Ep << "\n";

        fout4 << ch << "," << fixed << dif / s << "\n";
        
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
        //if (nV[i] < 450)
            fout3 << i * A << "," << nV[i] << "\n"; //Максвелл конец

    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    
    time (&end);
    cout << difftime(end, start) << endl;
    return 0;
}
