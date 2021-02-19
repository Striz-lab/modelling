#include <array>

struct proections(int s){
    double x[s];
    double y[s];
    double z[s];
};

class Particle{
    public:
        double x,     y,     z,            //координаты
               xpr,   ypr,   zpr,          //координаты в предыдущий момент
                                           //времени для уменьшения ошибки
               X,     Y,     Z,            //пока ничего лучше не придумал
               vx,    vy,    vz,           //скорости
               ax,    ay,    az,           //ускорения
               xnach, ynach, znach;        //начальные координаты (для дифф.)
    

    
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
    
        struct proections coordinates (double *x, double *y, double *z, int s){
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
        struct proections velocities (double *x, double *y, double *z, int s){
            for (int i = 0; i < s; i++){
                double V = pow(3*T, 0.5);              //инициализация направлений
                double ranx = ((double)rand() / RAND_MAX)*2*S - S;
                double rany = ((double)rand() / RAND_MAX)*2*S - S;
                double ranz = ((double)rand() / RAND_MAX)*2*S - S;
                double ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
            
                vx[i] = V * ranx / ran;                 //инициализация скоростей
                vy[i] = V * rany / ran;
                vz[i] = V * ranz / ran;
            }
        };
            
        struct proections previous (double *x, double *y, double *z, double *vx, double  *vy, double *vz, int s){
            for (int i = 0; i < s; i++){
                double V = pow(3*T, 0.5);              //инициализация направлений
                double ranx = ((double)rand() / RAND_MAX)*2*S - S;
                double rany = ((double)rand() / RAND_MAX)*2*S - S;
                double ranz = ((double)rand() / RAND_MAX)*2*S - S;
                double ran = pow(ranx * ranx + rany * rany + ranz * ranz, 0.5);
             
                xpr[i] = x[i] - vx[i] * dt;    //координаты частиц системы
                ypr[i] = y[i] - vy[i] * dt;    //на предыдущем шаге для
                zpr[i] = z[i] - vz[i] * dt;    //более точного рассчета
                
            
            
            xnach[i] = x[i];                        //инициализация начальных положений
            ynach[i] = y[i];                        //частичек
            znach[i] = z[i];
                }
        };
            
        struct ptoections accseleration (int s){
            for (int i = 0; i < s; i++){
                ax[i] = 0;                              //инициализация ускорений
                ay[i] = 0;
                az[i] = 0;
            }
        };
            
        struct proections first_coords (double *x, double *y, double *z, int s){
            for (int i = 0; i < s; i++){
                xnach[i] = x[i];
                ynach[i] = y[i];
                znach[i] = z[i];
            }
        };
    
        double* Maxwell(double *vx, double *vy = nullptr, double *vz = nullptr, int s){
            array<double, 1000> V;
            array<int, 10000> nV = {0};
            for (int i = 0; i < s; i++) //Максвелл начало (для последнего шага)
                V[i] = pow(vx[i]*vx[i] + vy[i]*vy[i] +vz[i]*vz[i], 0.5);

            double A = (V[s - 1] - V[0])/s, B = 0;
            for (int j = 0; j < s; j++){
                for (int i = 0; i < s; i++)
                    if ((V[i] != 0) and (V[i] <= B + 100*A) and (V[i] >= B - 100*A))
                        nV[j] ++;
                B += A;
            }
        }

};
