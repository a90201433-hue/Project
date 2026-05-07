#ifndef _LEVEL_SET_H_
#define _LEVEL_SET_H_

#include <vector>
#include <cmath>
#include <algorithm>
#include "Types.h"

class LevelSet {
private:
    Field& W;           // ссылка на поле, где хранится phi (индекс 3)
    int Nx_tot, Ny_tot;
    double dx, dy;
    double center_x, center_y, radius;
    double reinit_interval;  // как часто переинициализировать phi
    
public:
    LevelSet(Field& field, int Nx, int Ny, double dx_, double dy_)
        : W(field), Nx_tot(Nx), Ny_tot(Ny), dx(dx_), dy(dy_),
          reinit_interval(10) {}
    
    // Инициализация круглого пузыря
    void initCircle(double cx, double cy, double r, int fict) {
        center_x = cx;
        center_y = cy;
        radius = r;
        
        for (int i = 0; i < Nx_tot; i++) {
            for (int j = 0; j < Ny_tot; j++) {
                double x = i * dx;
                double y = j * dy;
                // Подписанное расстояние до окружности
                W[i][j][IDX_PHI] = sqrt((x - cx)*(x - cx) + (y - cy)*(y - cy)) - r;
            }
        }
    }
    
    // Получение массовой доли из phi (для уравнения состояния)
    double getMassFraction(double phi, double epsilon = 1.5*dx) {
        // Сглаженный Heaviside для плавного перехода
        if (phi < -epsilon) return 0.0;   // внутри пузыря
        if (phi > epsilon) return 1.0;    // снаружи
        // Плавный переход на границе
        return 0.5 * (1.0 + phi/epsilon + sin(M_PI * phi/epsilon) / M_PI);
    }
    
    // Вычисление нормали к границе
    void getNormal(int i, int j, double& nx, double& ny) {
        double phi = W[i][j][IDX_PHI];
        double dx_phi = (W[i+1][j][IDX_PHI] - W[i-1][j][IDX_PHI]) / (2.0 * dx);
        double dy_phi = (W[i][j+1][IDX_PHI] - W[i][j-1][IDX_PHI]) / (2.0 * dy);
        double len = sqrt(dx_phi*dx_phi + dy_phi*dy_phi);
        if (len > 1e-8) {
            nx = dx_phi / len;
            ny = dy_phi / len;
        } else {
            nx = ny = 0.0;
        }
    }
    
    // Вычисление кривизны (для поверхностного натяжения, опционально)
    double getCurvature(int i, int j) {
        double phi = W[i][j][IDX_PHI];
        double phix = (W[i+1][j][IDX_PHI] - W[i-1][j][IDX_PHI]) / (2.0 * dx);
        double phiy = (W[i][j+1][IDX_PHI] - W[i][j-1][IDX_PHI]) / (2.0 * dy);
        double phixx = (W[i+1][j][IDX_PHI] - 2*phi + W[i-1][j][IDX_PHI]) / (dx*dx);
        double phiyy = (W[i][j+1][IDX_PHI] - 2*phi + W[i][j-1][IDX_PHI]) / (dy*dy);
        double phixy = (W[i+1][j+1][IDX_PHI] - W[i+1][j-1][IDX_PHI] 
                      - W[i-1][j+1][IDX_PHI] + W[i-1][j-1][IDX_PHI]) / (4.0 * dx * dy);
        
        double grad2 = phix*phix + phiy*phiy;
        if (grad2 < 1e-8) return 0.0;
        
        return (phixx*phiy*phiy - 2*phix*phiy*phixy + phiyy*phix*phix) / pow(grad2, 1.5);
    }
    
    // Эволюция level set функции: ∂φ/∂t + u·∇φ = 0
    void advect(const Field& u_tilde, const Field& v_tilde, double dt, int fict) {
        Field phi_new = W;  // копируем
        
        for (int i = fict; i < Nx_tot - fict; i++) {
            for (int j = fict; j < Ny_tot - fict; j++) {
                // Интерполяция скорости на границы ячейки
                double u = u_tilde[i][j][0];
                double v = v_tilde[i][j][0];
                
                // Противопоточная схема первого порядка
                double phi_adv = W[i][j][IDX_PHI];
                
                // Адвекция по X
                if (u > 0) {
                    phi_adv -= u * dt / dx * (W[i][j][IDX_PHI] - W[i-1][j][IDX_PHI]);
                } else {
                    phi_adv -= u * dt / dx * (W[i+1][j][IDX_PHI] - W[i][j][IDX_PHI]);
                }
                
                // Адвекция по Y
                if (v > 0) {
                    phi_adv -= v * dt / dy * (W[i][j][IDX_PHI] - W[i][j-1][IDX_PHI]);
                } else {
                    phi_adv -= v * dt / dy * (W[i][j+1][IDX_PHI] - W[i][j][IDX_PHI]);
                }
                
                phi_new[i][j][IDX_PHI] = phi_adv;
            }
        }
        
        // Обновляем поле
        for (int i = 0; i < Nx_tot; i++) {
            for (int j = 0; j < Ny_tot; j++) {
                W[i][j][IDX_PHI] = phi_new[i][j][IDX_PHI];
            }
        }
    }
    
    // Переинициализация φ в функцию расстояния (алгоритм Sussman)
    void reinitialize(int fict, int num_iter = 10) {
        Field phi_new = W;
        
        for (int iter = 0; iter < num_iter; iter++) {
            for (int i = fict; i < Nx_tot - fict; i++) {
                for (int j = fict; j < Ny_tot - fict; j++) {
                    double phi = W[i][j][IDX_PHI];
                    
                    // Разности для ENO схемы
                    double phix_plus = (W[i+1][j][IDX_PHI] - phi) / dx;
                    double phix_minus = (phi - W[i-1][j][IDX_PHI]) / dx;
                    double phiy_plus = (W[i][j+1][IDX_PHI] - phi) / dy;
                    double phiy_minus = (phi - W[i][j-1][IDX_PHI]) / dy;
                    
                    // Выбор производной в зависимости от знака
                    double phix = (phi > 0) ? 
                        std::max(phix_minus, 0.0) + std::min(phix_plus, 0.0) :
                        std::min(phix_minus, 0.0) + std::max(phix_plus, 0.0);
                    
                    double phiy = (phi > 0) ? 
                        std::max(phiy_minus, 0.0) + std::min(phiy_plus, 0.0) :
                        std::min(phiy_minus, 0.0) + std::max(phiy_plus, 0.0);
                    
                    double grad = sqrt(phix*phix + phiy*phiy);
                    if (grad > 1e-8) {
                        phi_new[i][j][IDX_PHI] = phi - sign(phi) * (grad - 1.0) * dx;
                    }
                }
            }
            
            // Обновляем поле
            for (int i = 0; i < Nx_tot; i++) {
                for (int j = 0; j < Ny_tot; j++) {
                    W[i][j][IDX_PHI] = phi_new[i][j][IDX_PHI];
                }
            }
        }
    }
    
    // Функция знака
    double sign(double x) {
        if (x > 0) return 1.0;
        if (x < 0) return -1.0;
        return 0.0;
    }
    
    // Обновление массовой доли в W на основе phi
    void updateMassFraction(Field& W_phys, int fict) {
        double epsilon = 1.5 * dx;
        for (int i = fict; i < Nx_tot - fict; i++) {
            for (int j = fict; j < Ny_tot - fict; j++) {
                double phi = W[i][j][IDX_PHI];
                W_phys[i][j][3] = getMassFraction(phi, epsilon);
            }
        }
    }
    
    // Получение плотности по массовой доле
    double getDensity(double rho1, double rho2, double phi, double epsilon) {
        double W_frac = getMassFraction(phi, epsilon);
        return W_frac * rho1 + (1.0 - W_frac) * rho2;
    }
    
    // Сохранение границы для визуализации
    void saveBoundary(const std::string& filename, double t) {
        std::ofstream file(filename, std::ios::app);
        file << t;
        
        // Ищем нулевую изолинию (φ = 0)
        for (int i = 0; i < Nx_tot - 1; i++) {
            for (int j = 0; j < Ny_tot - 1; j++) {
                double phi00 = W[i][j][IDX_PHI];
                double phi10 = W[i+1][j][IDX_PHI];
                double phi01 = W[i][j+1][IDX_PHI];
                double phi11 = W[i+1][j+1][IDX_PHI];
                
                // Проверяем, есть ли пересечение нуля
                if (phi00 * phi10 < 0) {
                    double x_int = (i + 0.5) * dx;
                    double y_int = j * dy;
                    file << "," << x_int << "," << y_int;
                }
                if (phi00 * phi01 < 0) {
                    double x_int = i * dx;
                    double y_int = (j + 0.5) * dy;
                    file << "," << x_int << "," << y_int;
                }
            }
        }
        file << std::endl;
    }
};

#endif // _LEVEL_SET_H_