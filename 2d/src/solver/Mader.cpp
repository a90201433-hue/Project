#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include "Mader.h"
#include "FileProcessing.h"
#include "BoundCond.h"
#include "Init.h"
#include "Types.h"

extern double gamm, gamm1, Lx, Ly, T_init, R_gas, M, P_min, E_act, Z_freq, VISC, MINWT, GASW, MINGRHO;
extern int Nx, Ny, fict;

// Вспомогательная функция: вычисление эффективного γ по массовой доле (индекс 3)
inline double getGammaMix(double W_frac) {
    // W_frac = 1 → вещество 1 (gamm)
    // W_frac = 0 → вещество 2 (gamm1)
    return W_frac * gamm + (1.0 - W_frac) * gamm1;
}

// Фаза I: Уравнение состояния для двух веществ
void EOS(Field& W) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho = W[i][j][0];
            double W_frac = W[i][j][3];
            double P = W[i][j][4];
            
            double gamma_mix = getGammaMix(W_frac);
            
            // Внутренняя энергия из давления
            double I = P / ((gamma_mix - 1.0) * rho);
            
            // Пересчёт давления (защита от отрицательных значений)
            double P_new = (gamma_mix - 1.0) * rho * I;
            W[i][j][4] = (P_new < P_min) ? P_min : P_new;
        }
    }
}

// Фаза I: Химическая реакция (для нереагирующих веществ отключена)
void Arrenius(Field& W, double dt) {
    // Для НЕРЕАГИРУЮЩИХ веществ: ничего не делаем
    return;
}

// Фаза II: Искусственная вязкость (фон Неймана-Рихтмайера)
void Viscosity(Field& q1, Field& q2, Field& q3, Field& q4, const Field& W) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho = W[i][j][0];
            double u = W[i][j][1];
            double v = W[i][j][2];
            
            // Грань 1 (нижняя, направление Z): сжатие если V(i,j) >= V(i,j+1)
            q1[i][j][0] = (v >= W[i][j + 1][2]) 
                         ? VISC * rho * (v - W[i][j + 1][2]) 
                         : 0.0;
            // Грань 2 (левая, направление R): сжатие если U(i,j) >= U(i+1,j)
            q2[i][j][0] = (u >= W[i + 1][j][1]) 
                         ? VISC * rho * (u - W[i + 1][j][1]) 
                         : 0.0;
        }
    }
    
    // Грань 3 (верхняя) = нижняя грань ячейки (i, j-1)
    // Грань 4 (правая) = левая грань ячейки (i-1, j)
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            q3[i][j][0] = (j > fict) ? q1[i][j - 1][0] : 0.0;
            q4[i][j][0] = (i > fict) ? q2[i - 1][j][0] : 0.0;
        }
    }
}

// Фаза II: Обновление скоростей (полушаг)
void VelocityTilde(Field& u_tilde, Field& v_tilde, const Field& W,
                   const std::vector<double>& x, const std::vector<double>& y,
                   const Field& q1, const Field& q2, const Field& q3, const Field& q4,
                   double dt) {
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho = W[i][j][0];
            double u = W[i][j][1];
            double v = W[i][j][2];
            
            // Давления от соседей (индекс 4)
            double P1 = W[i][j - 1][4];
            double P2 = W[i - 1][j][4];
            double P3 = W[i][j + 1][4];
            double P4 = W[i + 1][j][4];

            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            // Обновление скорости V (вертикальная)
            v_tilde[i][j][0] = v - dt / (rho * 2.0 * dy) 
                              * ((P3 - P1) + (q3[i][j][0] - q1[i][j][0]));
            
            // Обновление скорости U (горизонтальная)
            u_tilde[i][j][0] = u - dt / (rho * 2.0 * dx) 
                              * ((P4 - P2) + (q4[i][j][0] - q2[i][j][0]));
        }
    }
}

// Фаза III: ZIP Energy Equation
void ZIPEnergy(Field& I_tilde, const Field& W, double dt,
               const Field& u_tilde, const Field& v_tilde,
               const Field& q1, const Field& q2, const Field& q3, const Field& q4,
               const std::vector<double>& x, const std::vector<double>& y) {
    
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho = W[i][j][0];
            double P = W[i][j][4];
            double W_frac = W[i][j][3];
            double gamma_mix = getGammaMix(W_frac);
            
            // Внутренняя энергия из давления
            double I = P / ((gamma_mix - 1.0) * rho);
            
            // Средние скорости на гранях
            double U1 = u_tilde[i - 1][j][0];
            double U2 = u_tilde[i + 1][j][0];
            double V1 = v_tilde[i][j - 1][0];
            double V2 = v_tilde[i][j + 1][0];
            double T3 = u_tilde[i][j][0];
            double T1 = v_tilde[i][j][0];
            
            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            I_tilde[i][j][0] = I 
                - dt / (4.0 * rho) 
                * ((P / dx) * (U2 - U1)
                 + (q4[i][j][0] / dx) * (U2 - T3)
                 + (q2[i][j][0] / dx) * (T3 - U1)
                 + (P / dy) * (V2 - V1)
                 + (q3[i][j][0] / dy) * (V2 - T1)
                 + (q1[i][j][0] / dy) * (T1 - V1));
        }
    }
}

// Поправка Шаргатова для контактной границы
double SharpatovCorrection(double W_donor, const Field& W, int i, int j, int Nx_tot, int Ny_tot) {
    const double EPS = 1e-6;
    
    // Если донор уже чистый — переносим его состав
    if (W_donor <= EPS) return 0.0;
    if (W_donor >= 1.0 - EPS) return 1.0;
    
    // Собираем соседей (i,j) — координаты донора
    std::vector<double> neighbors;
    
    // Левый сосед (i-1, j)
    if (i - 1 >= 0) neighbors.push_back(W[i - 1][j][3]);
    // Правый сосед (i+1, j)
    if (i + 1 < Nx_tot) neighbors.push_back(W[i + 1][j][3]);
    // Нижний сосед (i, j-1)
    if (j - 1 >= 0) neighbors.push_back(W[i][j - 1][3]);
    // Верхний сосед (i, j+1)
    if (j + 1 < Ny_tot) neighbors.push_back(W[i][j + 1][3]);
    // Диагональные соседи
    if (i - 1 >= 0 && j - 1 >= 0) neighbors.push_back(W[i - 1][j - 1][3]);
    if (i - 1 >= 0 && j + 1 < Ny_tot) neighbors.push_back(W[i - 1][j + 1][3]);
    if (i + 1 < Nx_tot && j - 1 >= 0) neighbors.push_back(W[i + 1][j - 1][3]);
    if (i + 1 < Nx_tot && j + 1 < Ny_tot) neighbors.push_back(W[i + 1][j + 1][3]);
    
    // Анализируем соседей
    bool has_pure_1 = false, has_pure_2 = false;
    for (double w : neighbors) {
        if (w >= 1.0 - EPS) has_pure_1 = true;
        if (w <= EPS) has_pure_2 = true;
    }
    
    // Если есть чистый сосед с веществом 1 — переносим вещество 2
    if (has_pure_1 && !has_pure_2) return 0.0;
    // Если есть чистый сосед с веществом 2 — переносим вещество 1
    if (has_pure_2 && !has_pure_1) return 1.0;
    
    // Иначе — по составу донора
    return W_donor;
}

// Фаза IV: Перенос массы (Donor-Acceptor) с поправкой Шаргатова
void ChangingFluxes(Field& DM, Field& DE, Field& DW, Field& DPU, Field& DPV,
                    const Field& u_tilde, const Field& v_tilde,
                    const std::vector<double>& x, const std::vector<double>& y,
                    double dt, const Field& W, const Field& I_tilde) {
    
    int Nx_tot = Nx + 2 * fict - 1;
    int Ny_tot = Ny + 2 * fict - 1;
    
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double dx = x[i] - x[i - 1];
            double dy = y[j] - y[j - 1];

            // Параметр переноса для R-направления
            double alpha = (0.5 * (u_tilde[i - 1][j][0] + u_tilde[i][j][0]) * dt / dx) 
                         / (1.0 + (u_tilde[i - 1][j][0] - u_tilde[i][j][0]) * dt / dx);
            
            // Параметр переноса для Z-направления
            double beta = (0.5 * (v_tilde[i][j - 1][0] + v_tilde[i][j][0]) * dt / dy) 
                        / (1.0 + (v_tilde[i][j - 1][0] - v_tilde[i][j][0]) * dt / dy);
 
            // ========== ПЕРЕНОС В НАПРАВЛЕНИИ R (влево-вправо) ==========
            if (alpha >= 0) {
                // Масса течёт слева направо: донор = (i-1,j), акцептор = (i,j)
                double rho_donor = W[i - 1][j][0];
                double DMASS = rho_donor * std::abs(alpha);
                
                double u_tilde_donor = u_tilde[i - 1][j][0];
                double v_tilde_donor = v_tilde[i - 1][j][0];
                double E_donor = I_tilde[i - 1][j][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
                double W_frac_donor = W[i - 1][j][3];
                
                double W_transfer = SharpatovCorrection(W_frac_donor, W, i - 1, j, Nx_tot, Ny_tot);
                
                DM[i][j][0]     += DMASS;
                DM[i - 1][j][0] -= DMASS;
                DE[i][j][0]     += DMASS * E_donor;
                DE[i - 1][j][0] -= DMASS * E_donor;
                DPU[i][j][0]    += DMASS * u_tilde_donor;
                DPU[i - 1][j][0] -= DMASS * u_tilde_donor;
                DPV[i][j][0]    += DMASS * v_tilde_donor;
                DPV[i - 1][j][0] -= DMASS * v_tilde_donor;
                DW[i][j][0]     += DMASS * W_transfer;
                DW[i - 1][j][0] -= DMASS * W_transfer;
            }
            else {
                // Масса течёт справа налево: донор = (i,j), акцептор = (i-1,j)
                double rho_donor = W[i][j][0];
                double DMASS = rho_donor * std::abs(alpha);
                
                double u_tilde_donor = u_tilde[i][j][0];
                double v_tilde_donor = v_tilde[i][j][0];
                double E_donor = I_tilde[i][j][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
                double W_frac_donor = W[i][j][3];
                
                double W_transfer = SharpatovCorrection(W_frac_donor, W, i, j, Nx_tot, Ny_tot);

                DM[i][j][0]     -= DMASS;
                DM[i - 1][j][0] += DMASS;
                DE[i][j][0]     -= DMASS * E_donor;
                DE[i - 1][j][0] += DMASS * E_donor;
                DPU[i][j][0]    -= DMASS * u_tilde_donor;
                DPU[i - 1][j][0] += DMASS * u_tilde_donor;
                DPV[i][j][0]    -= DMASS * v_tilde_donor;
                DPV[i - 1][j][0] += DMASS * v_tilde_donor;
                DW[i][j][0]     -= DMASS * W_transfer;
                DW[i - 1][j][0] += DMASS * W_transfer;
            }
            
            // ========== ПЕРЕНОС В НАПРАВЛЕНИИ Z (вверх-вниз) ==========
            if (beta >= 0) {
                // Масса течёт снизу вверх: донор = (i,j-1), акцептор = (i,j)
                double rho_donor = W[i][j - 1][0];
                double DMASS = rho_donor * std::abs(beta);
                
                double u_tilde_donor = u_tilde[i][j - 1][0];
                double v_tilde_donor = v_tilde[i][j - 1][0];
                double E_donor = I_tilde[i][j - 1][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
                double W_frac_donor = W[i][j - 1][3];
                
                double W_transfer = SharpatovCorrection(W_frac_donor, W, i, j - 1, Nx_tot, Ny_tot);

                DM[i][j][0]     += DMASS;
                DM[i][j - 1][0] -= DMASS;
                DE[i][j][0]     += DMASS * E_donor;
                DE[i][j - 1][0] -= DMASS * E_donor;
                DPU[i][j][0]    += DMASS * u_tilde_donor;
                DPU[i][j - 1][0] -= DMASS * u_tilde_donor;
                DPV[i][j][0]    += DMASS * v_tilde_donor;
                DPV[i][j - 1][0] -= DMASS * v_tilde_donor;
                DW[i][j][0]     += DMASS * W_transfer;
                DW[i][j - 1][0] -= DMASS * W_transfer;
            }
            else {
                // Масса течёт сверху вниз: донор = (i,j), акцептор = (i,j-1)
                double rho_donor = W[i][j][0];
                double DMASS = rho_donor * std::abs(beta);
                
                double u_tilde_donor = u_tilde[i][j][0];
                double v_tilde_donor = v_tilde[i][j][0];
                double E_donor = I_tilde[i][j][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
                double W_frac_donor = W[i][j][3];
                
                double W_transfer = SharpatovCorrection(W_frac_donor, W, i, j, Nx_tot, Ny_tot);

                DM[i][j][0]     -= DMASS;
                DM[i][j - 1][0] += DMASS;
                DE[i][j][0]     -= DMASS * E_donor;
                DE[i][j - 1][0] += DMASS * E_donor;
                DPU[i][j][0]    -= DMASS * u_tilde_donor;
                DPU[i][j - 1][0] += DMASS * u_tilde_donor;
                DPV[i][j][0]    -= DMASS * v_tilde_donor;
                DPV[i][j - 1][0] += DMASS * v_tilde_donor;
                DW[i][j][0]     -= DMASS * W_transfer;
                DW[i][j - 1][0] += DMASS * W_transfer;
            }
        }
    }
    
    // ========== ГРАНИЧНЫЕ УСЛОВИЯ: ПРАВАЯ ГРАНИЦА (CONTINUUM) ==========
    for (int j = fict; j < Ny + fict - 1; j++) {
        int i = Nx + fict - 2; 
        double dx = x[i] - x[i - 1];
        
        double alpha_right = u_tilde[i][j][0] * dt / dx;
        
        if (alpha_right > 0) { 
            double rho_donor = W[i][j][0];
            double DMASS = rho_donor * alpha_right;
            double u_tilde_donor = u_tilde[i][j][0];
            double v_tilde_donor = v_tilde[i][j][0];
            double E_d = I_tilde[i][j][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
            double W_frac_donor = W[i][j][3];
            
            DM[i][j][0]  -= DMASS;
            DE[i][j][0]  -= E_d * DMASS;
            DPU[i][j][0] -= u_tilde_donor * DMASS;
            DPV[i][j][0] -= v_tilde_donor * DMASS;
            DW[i][j][0]  -= W_frac_donor * DMASS;
        }
    }
    
    // ========== ГРАНИЧНЫЕ УСЛОВИЯ: ВЕРХНЯЯ ГРАНИЦА (CONTINUUM) ==========
    for (int i = fict; i < Nx + fict - 1; i++) {
        int j = Ny + fict - 2;
        double dy = y[j] - y[j - 1];
        
        double beta_top = v_tilde[i][j][0] * dt / dy;
        
        if (beta_top > 0) {
            double rho_donor = W[i][j][0];
            double DMASS = rho_donor * beta_top;
            double u_tilde_donor = u_tilde[i][j][0];
            double v_tilde_donor = v_tilde[i][j][0];
            double E_d = I_tilde[i][j][0] + 0.5 * (u_tilde_donor * u_tilde_donor + v_tilde_donor * v_tilde_donor);
            double W_frac_donor = W[i][j][3];
            
            DM[i][j][0]  -= DMASS;
            DE[i][j][0]  -= E_d * DMASS;
            DPU[i][j][0] -= u_tilde_donor * DMASS;
            DPV[i][j][0] -= v_tilde_donor * DMASS;
            DW[i][j][0]  -= W_frac_donor * DMASS;
        }
    }
}

// Фаза V: Перераспределение (Repartition)
void Repartition(Field& W, const Field& I_tilde, const Field& u_tilde, const Field& v_tilde,
                 const Field& DM, const Field& DE, const Field& DW, 
                 const Field& DPU, const Field& DPV) {
    
    for (int i = fict; i < Nx + fict - 1; i++) {
        for (int j = fict; j < Ny + fict - 1; j++) {
            double rho_old = W[i][j][0];
            double rho_new = rho_old + DM[i][j][0];

            // Защита от слишком малой плотности
            if (rho_new <= MINGRHO * 1e-3) {
                W[i][j][0] = 0.0;
                W[i][j][1] = 0.0;
                W[i][j][2] = 0.0;
                W[i][j][3] = 0.0;
                W[i][j][4] = 0.0;
                continue;
            }

            // Обновление скоростей
            double u_new = (rho_old * u_tilde[i][j][0] + DPU[i][j][0]) / rho_new;
            double v_new = (rho_old * v_tilde[i][j][0] + DPV[i][j][0]) / rho_new;
            
            // Обновление полной энергии и внутренней энергии
            double E_total = I_tilde[i][j][0] + 0.5 * (u_tilde[i][j][0] * u_tilde[i][j][0] 
                                                     + v_tilde[i][j][0] * v_tilde[i][j][0]);
            double I_new = (rho_old * E_total + DE[i][j][0]) / rho_new 
                         - 0.5 * (u_new * u_new + v_new * v_new);
            
            // Обновление массовой доли
            double W_frac_old = W[i][j][3];
            double W_frac_new = (rho_old * W_frac_old + DW[i][j][0]) / rho_new;
            W_frac_new = std::max(0.0, std::min(1.0, W_frac_new));
            
            // Пересчёт давления из внутренней энергии
            double gamma_mix = getGammaMix(W_frac_new);
            double P_new = (gamma_mix - 1.0) * rho_new * I_new;
            P_new = std::max(P_min, P_new);

            // Сохраняем новые значения
            W[i][j][0] = rho_new;
            W[i][j][1] = u_new;
            W[i][j][2] = v_new;
            W[i][j][3] = W_frac_new;
            W[i][j][4] = P_new;
        }
    }
}

// Главная функция метода Мейдера
void Mader(Field& W_new, const Field& W, const std::vector<double>& x, 
           const std::vector<double>& y, double dt) {
    
    int Nx_tot = Nx + 2 * fict - 1;
    int Ny_tot = Ny + 2 * fict - 1;

    // Создаём рабочую копию
    Field W_work = W;
    
    // Вспомогательные поля
    Field q1(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q2(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q3(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field q4(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    
    Field u_tilde(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field v_tilde(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field I_tilde(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    
    Field DM(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DE(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DW(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DPU(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));
    Field DPV(Nx_tot, std::vector<State>(Ny_tot, {0.0, 0.0, 0.0, 0.0}));

    // ========== ФАЗА I: Уравнение состояния ==========
    EOS(W_work);
    
    // ========== ФАЗА II: Искусственная вязкость ==========
    Viscosity(q1, q2, q3, q4, W_work);
    
    // ========== ФАЗА II: Обновление скоростей ==========
    VelocityTilde(u_tilde, v_tilde, W_work, x, y, q1, q2, q3, q4, dt);
    
    // ========== ФАЗА III: ZIP Energy Equation ==========
    ZIPEnergy(I_tilde, W_work, dt, u_tilde, v_tilde, q1, q2, q3, q4, x, y);
    
    // ========== ФАЗА IV: Перенос массы ==========
    ChangingFluxes(DM, DE, DW, DPU, DPV, u_tilde, v_tilde, x, y, dt, W_work, I_tilde);
    
    // ========== ФАЗА V: Перераспределение ==========
    Repartition(W_work, I_tilde, u_tilde, v_tilde, DM, DE, DW, DPU, DPV);
    
    // Копируем результат
    W_new = W_work;
}