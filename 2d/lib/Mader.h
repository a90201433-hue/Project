#ifndef _MADER_H_
#define _MADER_H_

#include "Types.h"

// Вспомогательная функция для вычисления эффективного γ по массовой доле
inline double getGammaMix(double W);

// Фаза I: Уравнение состояния
void EOS(Field& W);

// Фаза I: Химическая реакция (Аррениус)
void Arrenius(Field& W, double dt);

// Фаза II: Искусственная вязкость
void Viscosity(Field& q1, Field& q2, Field& q3, Field& q4, const Field& W);

// Фаза II: Обновление скоростей (полушаг)
void VelocityTilde(Field& u_tilde, Field& v_tilde, const Field& W,
                   const std::vector<double>& x, const std::vector<double>& y,
                   const Field& q1, const Field& q2, const Field& q3, const Field& q4,
                   double dt);

// Фаза III: ZIP Energy Equation
void ZIPEnergy(Field& I_tilde, Field& rho_tilde, const Field& W, double dt,
               const Field& u_tilde, const Field& v_tilde,
               const Field& q1, const Field& q2, const Field& q3, const Field& q4,
               const std::vector<double>& x, const std::vector<double>& y);

// Фаза IV: Перенос массы (Donor-Acceptor) с поправкой Шаргатова
void ChangingFluxes(Field& DM, Field& DE, Field& DW, Field& DPU, Field& DPV,
                    const Field& u_tilde, const Field& v_tilde,
                    const std::vector<double>& x, const std::vector<double>& y,
                    double dt, const Field& W, const Field& I_tilde, const Field& rho_tilde);

// Поправка Шаргатова для контактной границы
double SharpatovCorrection(double W_donor, const Field& W, int i, int j, int Nx_tot, int Ny_tot);

// Фаза V: Перераспределение (Repartition)
void Repartition(Field& W, const Field& I_tilde, const Field& u_tilde, const Field& v_tilde,
                 const Field& DM, const Field& DE, const Field& DW, 
                 const Field& DPU, const Field& DPV);

// Главная функция метода Мейдера
void Mader(Field& W_new, const Field& W, const std::vector<double>& x, 
           const std::vector<double>& y, double dt);

#endif // _MADER_H_