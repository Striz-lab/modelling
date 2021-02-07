#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <ctime>

using namespace std;

int s = 1000;               //количество частиц
double dt = 0.0001;         //время перерасчета
double S = 8;               //половина длины ребра куба
int N = 100000;             //количество итераций
double eps = 1, sigma = 1;  //параметры потенциала (лучше не менять)
double T = 3;               //аналог температуры
double Rmin = 0.7;          //минимально допустимое расстояние для генерации (квадрат, линейно -- 0.83)

double modulation(double x){//функция для периодических граничных условий
    if (x > S)              //и для перерасчета расстояний между частицами
        x = x - 2 * S;
    if (x < -S)
        x = x + 2 * S;
    return x;
}

double correction(double x, double y){ //аналогичная функция для коррекции диффузии
    if (x > S)
        y = y + 2 * S;
    if (x < -S)
        y = y - 2 * S;
    return y;
}

double rmin(double *x, double *y, double *z, int NOW){ //нахождение минимального
    double rmin = 100;                                 //расстояния между
        for (int j = 0; j < NOW; j++){                 //уже созданными частицами
            double rx = modulation(x[NOW] - x[j]);
            double ry = modulation(y[NOW] - y[j]);
            double rz = modulation(z[NOW] - z[j]);
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
    fout3.open("Velocities.csv");
    fout4.open("Diffusion.csv");
    fout5.open("Temperature.csv");

    double x[s],     y[s],     z[s],            //координаты
           xpr[s],   ypr[s],   zpr[s],          //координаты в предыдущий момент
                                                //времени для уменьшения ошибки
           X,        Y,        Z,               //пока ничего лучше не придумал
           vx[s],    vy[s],    vz[s],           //скорости
           ax[s],    ay[s],    az[s],           //ускорения
           xnach[s], ynach[s], znach[s],        //начальные координаты (для
           Ekin = 0, Ep = 0, dif = 0, a,        //диффузии)
           r = 0, rx, ry, rz;

    array<double, 1000> V;
    array<int, 100> nV;

    for (int i = 0; i < s; i++){                //случайная генерация частиц
        x[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        y[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        z[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
        double rm = 0;
        while (rm < Rmin){
            x[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            y[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            z[i] = ((double)rand() / RAND_MAX)*2*(S - Rmin) - (S - Rmin);
            
            rm = rmin(x,y,z,i);                 //проверка, чтобы не улетело далеко
        }
        //cout << rmin(x,y,z,i) << " " << i << endl;
    }
    for (int i = 0; i < s; i++){
        double V = pow(3*T, 0.5);              //инициализация направлений
        double ranx = ((double)rand() / RAND_MAX)*2*S - S;
        double rany = ((double)rand() / RAND_MAX)*2*S - S;
        double ranz = ((double)rand() / RAND_MAX)*2*S - S;
        double ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
         
        xpr[i] = x[i] - V * dt * ranx / ran;    //координаты частиц системы
        ypr[i] = y[i] - V * dt * rany / ran;    //на предыдущем шаге для
        zpr[i] = z[i] - V * dt * ranz / ran;    //более точного рассчета
        
        vx[i] = V * ranx / ran;                 //инициализация скоростей
        vy[i] = V * rany / ran;
        vz[i] = V * ranz / ran;

        ax[i] = 0;                              //инициализация ускорений
        ay[i] = 0;
        az[i] = 0;
        
        xnach[i] = x[i];                        //инициализация начальных положений
        ynach[i] = y[i];                        //частичек
        znach[i] = z[i];
    }
    
    for (int ch = 0; ch <= N; ch++){            //рассчет ускорений
        for (int i = 0; i < s; i++){
            for (int j = 0; j < i; j++){
                rx = modulation(x[i] - x[j]);
                ry = modulation(y[i] - y[j]);
                rz = modulation(z[i] - z[j]);
                r = rx * rx + ry * ry + rz * rz;

                a = 24*eps*(-2*pow((sigma/r), 7) + pow((sigma/r), 4))/(sigma*sigma);
                Ep = Ep + 4*eps*(pow((sigma/r), 6) - pow((sigma/r), 3));
                
                ax[i] = ax[i] - a * rx;
                ay[i] = ay[i] - a * ry;
                az[i] = az[i] - a * rz;
                
                ax[j] = ax[j] + a * rx;
                ay[j] = ay[j] + a * ry;
                az[j] = az[j] + a * rz;
            }
        }
        
        if (ch % 1000 == 0){
            fout1 << s << "\n" << ch << "\n";
            for (int i = 0; i < s; i++)
                fout1 << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
        }
        
        for (int i = 0; i < s; i++){        //обновление координат и скоростей
            
            X = x[i];
            Y = y[i];
            Z = z[i];
            
            x[i] = modulation(2 * x[i] - xpr[i] + ax[i] * dt * dt);//более точный
            y[i] = modulation(2 * y[i] - ypr[i] + ay[i] * dt * dt);//рассчет вместе
            z[i] = modulation(2 * z[i] - zpr[i] + az[i] * dt * dt);//с коррекцией
            
            xpr[i] = X;
            ypr[i] = Y;
            zpr[i] = Z;
            
            vx[i] = vx[i] + ax[i] * dt;
            vy[i] = vy[i] + ay[i] * dt;
            vz[i] = vz[i] + az[i] * dt;
            
            xnach[i] = correction(x[i] - xpr[i], xnach[i]); //коррекция
            ynach[i] = correction(y[i] - ypr[i], ynach[i]);
            znach[i] = correction(z[i] - zpr[i], znach[i]);
            
            ax[i] = 0;                                      //обновление
            ay[i] = 0;
            az[i] = 0;
            
            dif = dif + ((x[i] - xnach[i])) * ((x[i] - xnach[i])) + ((y[i] - ynach[i])) * ((y[i] - ynach[i])) + ((z[i] - znach[i])) * ((z[i] - znach[i]));
            Ekin = Ekin + (modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i])))/(2 * dt * dt);
        }
        if (ch % 10 == 0){
        fout2 << ch << "," << fixed << Ep << "," << fixed << Ekin << "," << fixed << Ekin + Ep << "\n";
        fout4 << ch << "," << fixed << dif / s << "\n";
        fout5 << ch << "," << fixed << 2*Ekin / (3*s) << "\n";
        }
        Ekin = 0;
        Ep = 0;
        dif = 0;
    }
    
    for (int i = 0; i < s; i++){ //Максвелл начало (для последнего шага)
        V[i] = pow(vx[i]*vx[i] + vy[i]*vy[i] +vz[i]*vz[i], 0.5);
        nV[i] = 0;
    }
    sort(V.begin(), V.end());
    double A = (V[s - 1] - V[0])/s;
    double B = 0;
    for (int j = 0; j < s; j++){
        for (int i = 0; i < s; i++)
            if ((V[i] != 0) and (V[i] <= B + 60*A) and (V[i] >= B - 60*A))
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
