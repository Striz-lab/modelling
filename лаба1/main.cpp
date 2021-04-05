#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <ctime>

using namespace std;

int s = 1000;               //количество частиц
double dt = 0.0001;         //время перерасчета
double S = 6;               //половина длины ребра куба
int N = 10000;             //количество итераций
double eps = 1, sigma = 1;  //параметры потенциала (лучше не менять)
double T = 0.6;               //аналог температуры
double Rmin = 0.7;          //минимально допустимое расстояние для генерации (квадрат, линейно -- 0.83)
double r_o = 2.5;

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
    
    ofstream fout1, fout2, fout3, fout4;
    
    fout1.open("Residuals.csv");
    fout2.open("Coor.xyz");
    fout3.open("Velocities.csv");
    fout4.open("Fluctuation.csv");

    
    array<double, 1000> V;
    array<int, 10000> nV;
    array<double, 1000> Vx;
    array<int, 10000> nVx;
    array<double, 1000> Vy;
    array<int, 10000> nVy;
    array<double, 1000> Vz;
    array<int, 10000> nVz;
    
    double x[s],     y[s],     z[s],            //координаты
           xpr[s],   ypr[s],   zpr[s],          //координаты в предыдущий момент
                                                //времени для уменьшения ошибки
           X,        Y,        Z,               //пока ничего лучше не придумал
           vx[s],    vy[s],    vz[s],           //скорости
           ax[s],    ay[s],    az[s],           //ускорения
           xnach[s], ynach[s], znach[s],
           vxnach[s], vynach[s], vznach[s],     //начальные координаты (для
           Ekin = 0, Ep = 0, dif = 0, a,        //диффузии)
           F1 = 0, f = 0,
           r = 0, rx, ry, rz,
           deltaR = 0, deltaV = 0,  P = 0, Pm = 0, Eav2 = 0, Eav = 0;
    
    double x1[s],     y1[s],     z1[s],            //координаты
        xpr1[s],   ypr1[s],   zpr1[s],          //координаты в предыдущий момент
                                         //времени для уменьшения ошибки
        X1,        Y1,        Z1,               //пока ничего лучше не придумал
        vx1[s],    vy1[s],    vz1[s],           //скорости
        ax1[s],    ay1[s],    az1[s],           //ускорения
        xnach1[s], ynach1[s], znach1[s],
        vxnach1[s], vynach1[s], vznach1[s],     //начальные координаты (для
         a1,        //диффузии)
        r1 = 0, rx1, ry1, rz1,
        deltaR1 = 0, deltaV1 = 0;

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
        
        x1[i] = x[i];
        y1[i] = y[i];
        z1[i] = z[i];
        
        xpr1[i] = x1[i] - V * dt * ranx / ran;    //координаты частиц системы
        ypr1[i] = y1[i] - V * dt * rany / ran;    //на предыдущем шаге для
        zpr1[i] = z1[i] - V * dt * ranz / ran;    //более точного рассчета
        
        vx[i] = V * ranx / ran;                 //инициализация скоростей
        vy[i] = V * rany / ran;
        vz[i] = V * ranz / ran;
        
        vx1[i] = vx[i];
        vy1[i] = vy[i];
        vz1[i] = vz[i];

        ax[i] = 0;                              //инициализация ускорений
        ay[i] = 0;
        az[i] = 0;
        
        ax1[i] = 0;
        ay1[i] = 0;
        az1[i] = 0;
        
        xnach[i] = x[i];                        //инициализация начальных положений
        ynach[i] = y[i];                        //частичек
        znach[i] = z[i];
        
        xnach1[i] = x1[i];                        //инициализация начальных положений
        ynach1[i] = y1[i];                        //частичек
        znach1[i] = z1[i];
        
        vxnach[i] = vx[i];
        vynach[i] = vy[i];
        vznach[i] = vz[i];
        
    }
    
    for (int ch = 0; ch <= N; ch++){            //рассчет ускорений
        for (int i = 0; i < s; i++){
            for (int j = 0; j < i; j++){
                rx = modulation(x[i] - x[j]);
                ry = modulation(y[i] - y[j]);
                rz = modulation(z[i] - z[j]);
                r = rx * rx + ry * ry + rz * rz;

                //if (r <= r_o){
                    a = 24*eps*(-2*pow((sigma/r), 7) + pow((sigma/r), 4))/(sigma*sigma);
                    Ep = Ep + 4*eps*(pow((sigma/r), 6) - pow((sigma/r), 3));
                    
                    ax[i] = ax[i] - a * rx;
                    ay[i] = ay[i] - a * ry;
                    az[i] = az[i] - a * rz;
                    
                    ax[j] = ax[j] + a * rx;
                    ay[j] = ay[j] + a * ry;
                    az[j] = az[j] + a * rz;
                //}
                
                if ((ch > 4000) and(ch %10 == 0)){
                    P = P + r*(ax[i]*ax[i] + ay[i]*ay[i] + az[i]*az[i]);
                    //cout << P << endl;
                }
                    
            }
        }
        
        
        for (int i = 0; i < s; i++){        //обновление координат и скоростей
            X = x[i];
            Y = y[i];
            Z = z[i];
            if ((ch > 4000) and (ch%10 == 0)){
            
                Pm += (modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i])))/(2 * dt * dt);
                
            }
             
            Ekin = Ekin + (modulation((x[i] - xpr[i])) * modulation((x[i] - xpr[i])) + modulation((y[i] - ypr[i])) * modulation((y[i] - ypr[i])) + modulation((z[i] - zpr[i])) * modulation((z[i] - zpr[i])))/(2 * dt * dt);
            
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
            
            
            //if (ch> 4000){
            if (ch%20 == 0){
                Eav2 += (Ep+Ekin)*(Ep+Ekin);
                Eav += (Ep+Ekin);
            }
            //}
            
            if (ch%10 == 0)
                fout4 << ch << "," << (Eav2 - Eav*Eav)/((ch)*(ch)) << endl;
            
        }
        
        if (ch % 10 == 0){
            fout2 << s << "\n" << ch << "\n";
            for (int i = 0; i < s; i++)
                fout2 << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << "\n";
        }
        
        if (ch % 10 == 0){
            for (int i = 0; i < s; i++){
                for (int j = 0; j < i; j++){
                    rx1 = modulation(x1[i] - x1[j]);
                    ry1 = modulation(y1[i] - y1[j]);
                    rz1 = modulation(z1[i] - z1[j]);
                    r1 = rx1 * rx1 + ry1 * ry1 + rz1 * rz1;

                    a1 = 24*eps*(-2*pow((sigma/r1), 7) + pow((sigma/r1), 4))/(sigma*sigma);
                    
                    ax1[i] = ax1[i] - a1 * rx1;
                    ay1[i] = ay1[i] - a1 * ry1;
                    az1[i] = az1[i] - a1 * rz1;
                    
                    ax1[j] = ax1[j] + a1 * rx1;
                    ay1[j] = ay1[j] + a1 * ry1;
                    az1[j] = az1[j] + a1 * rz1;
                }
            }
            
            
            for (int i = 0; i < s; i++){        //обновление координат и скоростей
                
                
                
                
                X1 = x1[i];
                Y1 = y1[i];
                Z1 = z1[i];
                
                x1[i] = modulation(2 * x1[i] - xpr1[i] + 100 * ax1[i] * dt * dt);//более точный
                y1[i] = modulation(2 * y1[i] - ypr1[i] + 100 * ay1[i] * dt * dt); //рассчет вместе
                z1[i] = modulation(2 * z1[i] - zpr1[i] + 100 * az1[i] * dt * dt);//с коррекцией
                
                xpr1[i] = X1;
                ypr1[i] = Y1;
                zpr1[i] = Z1;
                
                //xnach1[i] = correction(x1[i] - xpr1[i], xnach1[i]); //коррекция
                //ynach1[i] = correction(y1[i] - ypr1[i], ynach1[i]);
                //znach1[i] = correction(z1[i] - zpr1[i], znach1[i]);
                
                vx1[i] = vx1[i] + 10 * ax1[i] * dt;
                vy1[i] = vy1[i] + 10 * ay1[i] * dt;
                vz1[i] = vz1[i] + 10 * az1[i] * dt;
                
                ax1[i] = 0;                                      //обновление
                ay1[i] = 0;
                az1[i] = 0;
                
                //deltaR1 = deltaR1 + pow(pow(x1[i] * x1[i] + y1[i] * y1[i] + z1[i] * z1[i],0.5),2);// - pow(x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 0.5), 2);
                deltaR1 = deltaR1 + ((x[i] - xnach[i])) * ((x[i] - xnach[i])) + ((y[i] - ynach[i])) * ((y[i] - ynach[i])) + ((z[i] - znach[i])) * ((z[i] - znach[i]));//
                //deltaR = deltaR + ((x1[i] - xnach1[i])) * ((x1[i] - xnach1[i])) + ((y1[i] - ynach1[i])) * ((y1[i] - ynach1[i])) + ((z1[i] - znach1[i])) * ((z1[i] - znach1[i]));
                
                deltaV1 = deltaV1 + pow(pow(vx1[i] * vx1[i] + vy1[i] * vy1[i] + vz1[i] * vz1[i], 0.5) - pow(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i], 0.5), 2);
                
                //deltaR = deltaR + pow(pow(x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 0.5), 2);
                
//deltaV = deltaV + pow(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i], 0.5);
                //deltaR = deltaR + pow(x[i] * x[i] + y[i] * y[i] + z[i] * z[i], 0.5);
                
            }
    
            fout1 << ch << "," << fixed << Ekin << "," << fixed << deltaV1/(s*s) << "\n";
            deltaR = 0;
            deltaV = 0;
            deltaR1 = 0;
            deltaV1 = 0;
            Ep = 0;
            Ekin = 0;
       
            //P = 0;
            //fout4 << ch << "," << (Eav2 - Eav*Eav)/((ch-4000)*(ch-4000) << endl;
        }
        if (ch % 1000 == 0)
            cout << (double)ch/N*100 << "%" << endl;
        
    }
    
    //cout << P/(600*3*pow(2*S,3)*1000) << " " << 2*s*Pm/(600*3*pow(2*S,3)) << " " << P/(600*3*pow(2*S,3)*1000) + 2*s*Pm/(600*3*pow(2*S,3)) << " " << (2*s*Pm/(600*3*pow(2*S,3)))/(P/(600*3*pow(2*S,3)*1000) + 2*s*Pm/(600*3*pow(2*S,3)))<< " " << s/(8*S*S*S) << endl;
    
    //cout << P/(600*3*pow(2*S,3)*1000) + 2*s*Pm/(600*3*pow(2*S,3)) << endl;
    
    cout << (Eav2 - Eav*Eav)/(N*N*N) << endl;
    /*
    for (int i = 0; i < s; i++){ //Максвелл начало (для последнего шага)
        V[i] = pow(vx[i]*vx[i] + vy[i]*vy[i] +vz[i]*vz[i], 0.5);
        Vx[i] = vx[i];
        Vy[i] = vy[i];
        Vz[i] = vz[i];
        nV[i] = 0;
        
        nVx[i] = 0;
        nVy[i] = 0;
        nVz[i] = 0;
    }
    sort(V.begin(), V.end());
    sort(Vx.begin(), Vx.end());
    sort(Vy.begin(), Vy.end());
    sort(Vz.begin(), Vz.end());
    double A = (V[s - 1] - V[0])/s, Ax = (Vx[s - 1] - Vx[0])/s, Ay = (Vy[s - 1] - Vy[0])/s, Az = (Vz[s - 1] - Vz[0])/s;
    double B = 0, Bx = 0, By = 0, Bz = 0;
    for (int j = 0; j < s; j++){
        for (int i = 0; i < s; i++){
            if ((V[i] != 0) and (V[i] <= B + 100*A) and (V[i] >= B - 100*A))
                nV[j] ++;
            if ((Vx[i] != 0) and (Vx[i] <= B + 100*Ax) and (Vx[i] >= B - 100*Ax))
                nVx[j] ++;
            if ((Vy[i] != 0) and (Vy[i] <= B + 100*Ay) and (Vy[i] >= B - 100*Ay))
                nVy[j] ++;
            if ((Vz[i] != 0) and (Vz[i] <= B + 100*Az) and (Vz[i] >= B - 100*Az))
                nVz[j] ++;
            }
        B += A;
        Bx += Ax;
        By += Ay;
        Bz += Az;
    }

    for (int i = 0; i < s; i++)
        fout3 << i * A << "," << nV[i] << "," << nVx[i] << "," << nVy[i] << "," << nVz[i] << "," << i * Ax << "," << i * Ay << "," << i * Az << "\n";*/
    
    fout1.close();
    fout2.close();
    fout3.close();
    fout4.close();
    
    time (&end);
    cout << difftime(end, start) << endl;
    return 0;
}
