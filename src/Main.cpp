#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <filesystem>
#include "lib/Init.h"
#include "lib/BoundCond.h"
#include "lib/RiemannSolver.h"
#include "lib/Acoustic.h"
#include "lib/GeneralFunctions.h"

int N, fo, step_max, bound_case;
int fict = 3;
double L, t_max, x0, gamm, CFL;
std::string x_left_bound, x_right_bound, Soda;

namespace fs = std::filesystem;

void GetDt(std::vector<std::vector<double>> W, std::vector<double> x, double& dt){
	double c = 0, dx = L, u = 0;
	double c_temp;

	for (int i = fict; i < N + fict; i++){
		c_temp = std::sqrt(gamm * W[i][2] / W[i][0]);
		if (c_temp > c) 
			c = c_temp;
		if (std::abs(W[i][1]) > u)
			u = std::abs(W[i][1]);
		if ((x[i + 1] - x[i]) < dx)
			dx = x[i + 1] - x[i];
	}

	dt = CFL * dx/(c + u);
	//std::cout << "dt= " << dt << std::endl;
}

void WriteToCSV(std::vector<std::vector<double>> W, std::vector<double> xc, double t, std::ofstream& file){
	file << "time,x,rho,u,P,e" << std::endl;
	std::vector<double> e1(N + 2 * fict - 1);
	double e;
	for (int i = 0; i < N + 2 * fict - 1; i++) {
		e1[i] = W[i][2] / (W[i][0] * (gamm - 1));
	}
	double e0 = *std::min_element(e1.begin(), e1.end());
	
	for (int i = 0; i < N + 2 * fict - 1; i++){
		if (W[i][0] < 0.03) {
			e = e0;
		}
		else {
			e = W[i][2] / (W[i][0] * (gamm - 1));
		}
		
		
		file << t << "," << xc[i] 
			<< "," << W[i][0] 
			<< "," << W[i][1] 
			<< "," << W[i][2] 
			<< "," << e << std::endl;
	}
}

void MoveToChache(const std::string& source_root_path, const std::string& cache_root_path) {

	fs::path source_path = source_root_path;
    	fs::path cache_path = cache_root_path;
    	bool all_successful = true;

    	try {
        	// 1. Проверка существования исходной папки
        	if (!fs::exists(source_path) || !fs::is_directory(source_path)) {
            	std::cerr << "Ошибка: Исходная папка не существует или не является папкой: " << source_root_path << std::endl;
            	return;
        	}

        	// 2. Убедиться, что папка Кэша существует. Если нет, создать ее.
        	fs::create_directories(cache_path);

		std::cout << "Очистка старых данных в папке Кэша: " << cache_root_path << std::endl;
        	fs::remove_all(cache_path);

        	// 3. Итерация по содержимому исходной папки
        	for (const auto& entry : fs::directory_iterator(source_path)) {
            		fs::path current_item = entry.path();
            		// Получаем имя перемещаемого элемента (например, "Run_2025_A")
            		fs::path item_filename = current_item.filename();
            		// Формируем новый путь в папке Кэша
            		fs::path destination_path = cache_path / item_filename;

           	 	try {
               	 		// Перемещаем элемент
                		fs::rename(current_item, destination_path);
                		std::cout << "Перемещено: " << current_item.string() << " -> " << destination_path.string() << std::endl;

            		} catch (const fs::filesystem_error& e) {
                	// Если возникает ошибка при перемещении одного элемента
                		std::cerr << "Ошибка перемещения " << current_item.string() << ": " << e.what() << std::endl;
                		all_successful = false; // Отмечаем, что не все прошло успешно
            		}
        	}

    	} catch (const std::exception& e) {
        	std::cerr << "Произошла критическая ошибка: " << e.what() << std::endl;
        	return;
    	}

    	return;
}



void CleanDir(std::vector<std::string> folders){
	for (const auto& folder : folders) {
        	if (!fs::exists(folder)) {
            		std::cout << "Папка не существует: " << folder << std::endl;
            		continue;
        	}

        	for (const auto& entry : fs::directory_iterator(folder)) {
            		if (fs::is_regular_file(entry)) {
                		try {
                    			fs::remove(entry);
                    			
                		} catch (const std::exception& e) {
                    			std::cerr << "Ошибка при удалении " << entry.path() 
                              		<< ": " << e.what() << std::endl;
                		}
            		}	
        	}
    	}

	std::cout << "Очистка завершена." << std::endl;
	return;	
}

void GetCommonDt(std::vector<std::vector<double>> W_G, std::vector<std::vector<double>> W_GK, std::vector<std::vector<double>> W_GKR, std::vector<std::vector<double>> W_ENO, std::vector<std::vector<double>> W_WENO, std::vector<double> x, double& dt_common) {

	double dt_G, dt_GK, dt_GKR, dt_ENO, dt_WENO;

	GetDt(W_G, x, dt_G);
	GetDt(W_GK, x, dt_GK);
	GetDt(W_GKR, x, dt_GKR);
	GetDt(W_ENO, x, dt_ENO);
	GetDt(W_WENO, x, dt_WENO);

	dt_common = std::min(std::min(std::min(std::min(dt_G, dt_GK), dt_GKR), dt_ENO), dt_WENO);
	//std::cout << "dt_common= " << dt_common << std::endl;

}

void printProgressBar(double t, double t_max) {
    const int total_cells = 10; // всего 10 ячеек
    int filled_cells = static_cast<int>((t / t_max) * total_cells);

    std::cout << "\r";

    for (int i = 0; i < total_cells; ++i) {
        if (i < filled_cells)
            std::cout << "\033[92m\u25A0\033[0m"; // заполненная ячейка
        else
            std::cout <<"▱"; // пустая ячейка
    }

    int percent = static_cast<int>((t / t_max) * 100);
    std::cout << " " << percent << "%";
    std::cout.flush();
}

#include <iostream>
#include <random>
#include <vector>
#include <algorithm>

void printProgressTree(double t, double t_max) {
    const int H = 6;
    const int trunk_h = 2;
    const int trunk_w = 3;
    const int total_rows = H + trunk_h + 1;
    static bool first_call = true;

    // Задаём количество лампочек внутри функции
    const int num_lamps = 10; // <-- здесь можно менять число лампочек

    // Сколько всего листьев
    int total_leaves = H * H;

    // --- вычисляем прогресс ---
    double p = t / t_max;
    if (p < 0) p = 0;
    if (p > 1) p = 1;
    int percent = static_cast<int>(p * 100);
    int filled = static_cast<int>(p * total_leaves);

    // --- фаза гирлянды (0..10) ---
    int phase = percent / 10;

    // --- генератор случайных чисел ---
    static std::mt19937 rng(std::random_device{}());
    auto randColor = [&](int x) {
        switch (x % 5) {
            case 0: return "\033[31m"; // красный
            case 1: return "\033[33m"; // жёлтый
            case 2: return "\033[36m"; // голубой
            case 3: return "\033[35m"; // фиолетовый
            default: return "\033[91m"; // ярко-красный
        }
    };

    // --- храним гирлянду ---
    static std::vector<bool> lamps_mask(total_leaves, false);
    static std::vector<int> lamps_color(total_leaves, 0);
    static int current_phase = -1;

    // если перешли в новую фазу — перегенерируем гирлянду!
    if (phase != current_phase) {
        current_phase = phase;

        // сначала очищаем маску
        std::fill(lamps_mask.begin(), lamps_mask.end(), false);

        // выбираем случайные позиции для лампочек
        std::vector<int> indices(total_leaves);
        for (int i = 0; i < total_leaves; ++i) indices[i] = i;
        std::shuffle(indices.begin(), indices.end(), rng);

        for (int i = 0; i < num_lamps && i < total_leaves; ++i) {
            lamps_mask[indices[i]] = true;
        }

        // генерируем случайные цвета для всех листьев
        for (int i = 0; i < total_leaves; ++i) {
            lamps_color[i] = rng();
        }
    }

    // --- на первом вызове создаём пустое пространство ---
    if (first_call) {
        for (int i = 0; i < total_rows; ++i)
            std::cout << "\n";
        first_call = false;
    }

    // --- поднимаем курсор ---
    std::cout << "\033[" << total_rows << "A";

    int max_width = 2 * (H - 1) + 1;
    int counter = 0;

    // --- Листья ---
    for (int i = 0; i < H; ++i) {
        int width = 2 * i + 1;
        int pad = (max_width - width) / 2;
        std::cout << std::string(pad, ' ');
        for (int j = 0; j < width; ++j) {
            bool isLamp = lamps_mask[counter];
            int color = lamps_color[counter];
            if (counter < filled) {
                if (isLamp) std::cout << randColor(color) << "\033[42m" << "●" << "\033[0m";
                else std::cout << "\033[32m█\033[0m";
            } else {
                if (isLamp) std::cout << randColor(color) << "○" << "\033[0m";
                else std::cout << "░";
            }
            counter++;
        }
        std::cout << "\n";
    }

    // --- Ствол ---
    int pad = (max_width - trunk_w) / 2;
    for (int i = 0; i < trunk_h; ++i) {
        std::cout << std::string(pad, ' ');
        std::cout << "\033[38;5;130m";
        for (int j = 0; j < trunk_w; ++j) std::cout << "█";
        std::cout << "\033[0m\n";
    }

    // --- Процент ---
    std::cout << "Прогресс: " << percent << "%\n";
    std::cout.flush();
}

int main() {
	//MoveToChache("CSV_files/ActualCalc", "CSV_files/ChacheCalc");
	std::vector<std::string> folders = {"CSV_files/Godunov",
					    "CSV_files/Kolgan",
 					    "CSV_files/Kolgan2",
					    "CSV_files/Rodionov",
					    "CSV_files/Rodionov2",
					    "CSV_files/ENO",
					    "CSV_files/WENO"};
	CleanDir(folders);
	readConfig();
	
	//std::cout << bound_name << std::endl;
	double dx = L / (N - 1);
	std::vector<double> xc(N + 2 * fict - 1);
	std::vector<double> x(N + 2 * fict);	
	Grid(dx, x, xc);

	double rho_L, u_L, P_L;
	double rho_R, u_R, P_R;
		
	std::vector<std::vector<double>> W_0(N + 2*fict - 1);
	InitValues(Soda, rho_L, u_L, P_L, rho_R, u_R, P_R, t_max, x0, xc, W_0);
	
	std::vector<double> W_L_riemann = {rho_L, u_L, P_L};
	std::vector<double> W_R_riemann = {rho_R, u_R, P_R};
	std::vector<double> W_star_riemann = {0, 0};
	std::vector<double> W_riemann = {0, 0, 0};
	double e_riemann;

	std::vector<double> W_star_ac = {0, 0};
	std::vector<double> W_ac = {0, 0, 0};
	std::vector<std::vector<double>> W_ac1(N + 2 * fict -1);
	InitialZeros(W_ac1, 3);
	double e_ac;

	//t_max = t_riemann;	
	

	// I. Решатель Римана
	NewtonForPressure(W_L_riemann, W_R_riemann, W_star_riemann, 1e-6);
	
	std::ofstream file("CSV_files/Riemann.csv");

	file << "x,rho,u,P,e" << std::endl;

	for (int i = fict; i < N + fict; i++){
		W_riemann = GetParamsFromChoosingWave(W_L_riemann, W_R_riemann, W_star_riemann, xc[i] - x0, t_max);
		e_riemann = W_riemann[2] / (W_riemann[0] * (gamm - 1));
		file << xc[i]  << "," << W_riemann[0] << "," << W_riemann[1] << "," << W_riemann[2] << "," << e_riemann << std::endl;
	}
	file.close();


	// II. Акустика
  	AcoustPStarAndUStar(W_L_riemann, W_R_riemann, W_star_ac);
  	
  	std::ofstream fileAc("CSV_files/acoust.csv");
  	
	fileAc << "x,rho,u,P,e" << std::endl;

	std::vector<double> e2(N + 2 * fict - 1);
	for (int i = fict; i < N + fict; i++) {
		W_ac1[i] = ConfigurationAcousticWaves(W_L_riemann, W_R_riemann, W_star_ac, x[i] - x0, t_max);
		e2[i] = W_ac1[i][2] / (W_ac1[i][0] * (gamm - 1));
	}
	double e2_0 = *std::min_element(e2.begin(), e2.end());
	for (int i = fict; i < N + fict; i++){
		if (W_ac1[i][0] < 0.03) {
			e_ac = e2_0;
		}
		else {
			e_ac = W_ac1[i][2] / (W_ac1[i][0] * (gamm - 1));
		}
		fileAc << x[i] << "," <<  W_ac1[i][0] << "," <<  W_ac1[i][1] << "," <<  W_ac1[i][2] << "," << e_ac << std::endl;
  	}
  	fileAc.close();		


	// For all methods
	Grid(dx, x, xc);


	// For Godunov
	std::vector<std::vector<double>> W_G = W_0;
	std::vector<std::vector<double>> W_new_G(N + 2 * fict - 1);
	InitialZeros(W_new_G, 3);
	std::vector<std::vector<double>> W_b_G(N + 2 * fict);
	InitialZeros(W_b_G, 3);
	std::vector<std::vector<double>> F_G(N + 2 * fict);
	InitialZeros(F_G, 2);
	
	// For Kolgan
	std::vector<std::vector<double>> W_GK = W_0;
	std::vector<std::vector<double>> W_new_GK(N + 2 * fict - 1);
	InitialZeros(W_new_GK, 3);
	std::vector<std::vector<double>> W_b_GK(N + 2 * fict);
	InitialZeros(W_b_GK, 3);
	std::vector<std::vector<double>> F_GK(N + 2 * fict);
	InitialZeros(F_GK, 2);
	
	// For Kolgan 2
	std::vector<std::vector<double>> W_GK2 = W_0;
	std::vector<std::vector<double>> W_new_GK2(N + 2 * fict - 1);
	InitialZeros(W_new_GK2, 3);
	std::vector<std::vector<double>> W_b_GK2(N + 2 * fict);
	InitialZeros(W_b_GK2, 3);
	std::vector<std::vector<double>> F_GK2(N + 2 * fict);
	InitialZeros(F_GK2, 2);


	// For Rodionov
	std::vector<std::vector<double>> W_GKR = W_0;	
	std::vector<std::vector<double>> W_new_GKR(N + 2 * fict - 1);
	InitialZeros(W_new_GKR, 3);
	std::vector<std::vector<double>> W_b_GKR(N + 2 * fict);
	InitialZeros(W_b_GKR, 3);
	std::vector<std::vector<double>> F_GKR(N + 2 * fict);
	InitialZeros(F_GKR, 2);

	// For Rodionov 2
	std::vector<std::vector<double>> W_GKR2 = W_0;	
	std::vector<std::vector<double>> W_new_GKR2(N + 2 * fict - 1);
	InitialZeros(W_new_GKR2, 3);
	std::vector<std::vector<double>> W_b_GKR2(N + 2 * fict);
	InitialZeros(W_b_GKR2, 3);
	std::vector<std::vector<double>> F_GKR2(N + 2 * fict);
	InitialZeros(F_GKR2, 2);

	// For ENO
	std::vector<std::vector<double>> W_ENO = W_0;
	std::vector<std::vector<double>> W_new_ENO(N + 2 * fict - 1);
	InitialZeros(W_new_ENO, 3);
	std::vector<std::vector<double>> W_b_ENO(N + 2 * fict);
	InitialZeros(W_b_ENO, 3);
	std::vector<std::vector<double>> F_ENO(N + 2 * fict);
	InitialZeros(F_ENO, 3);

	// For WENO
	std::vector<std::vector<double>> W_WENO = W_0;
	std::vector<std::vector<double>> W_new_WENO(N + 2 * fict - 1);
	InitialZeros(W_new_WENO, 3);
	std::vector<std::vector<double>> W_b_WENO(N + 2 * fict);
	InitialZeros(W_b_WENO, 3);
	std::vector<std::vector<double>> F_WENO(N + 2 * fict);
	InitialZeros(F_WENO, 3);

	// For all methods
	double t = 0.0;
	double dt_common;
	int step = 0;

	while (t <= t_max) {
		
		GetDt(W_WENO, xc, dt_common);
		//GetCommonDt(W_G, W_GK, W_GKR, W_ENO, W_WENO, x, dt_common);
		
		if (t + dt_common > t_max) {
			dt_common = t_max - t;
		}
		t += dt_common;
		//printProgressBar(t, t_max);
		printProgressTree(t, t_max);
		/*UpdateArrays(W_G, W_new_G, W_b_G, F_G, x, dt_common, "Godunov");
		UpdateArrays(W_GK, W_new_GK, W_b_GK, F_GK, x, dt_common, "Kolgan");
		UpdateArrays(W_GK2, W_new_GK2, W_b_GK2, F_GK2, x, dt_common, "Kolgan2");
		UpdateArrays(W_GKR, W_new_GKR, W_b_GKR, F_GKR, x, dt_common, "Rodionov");
		UpdateArrays(W_GKR2, W_new_GKR2, W_b_GKR2, F_GKR2, x, dt_common, "Rodionov2");*/
		//UpdateArrays(W_ENO, W_new_ENO, W_b_ENO, F_ENO, x, dt_common, "ENO");	
		UpdateArrays(W_WENO, W_new_WENO, W_b_WENO, F_WENO, x, dt_common, "WENO", "R3");
		/*BoundCond(W_G);
		BoundCond(W_GK);
		BoundCond(W_GK2);
		BoundCond(W_GKR);
		BoundCond(W_GKR2);*/
		//BoundCond(W_ENO);
		BoundCond(W_WENO);

		if (step % fo == 0){
			/*std::string filename = "CSV_files/Godunov/step_" + std::to_string(step) + ".csv";
			std::ofstream file(filename);
			WriteToCSV(W_G, xc, t, file);	
			file.close();
			
			filename = "CSV_files/Kolgan/step_" + std::to_string(step) + ".csv";
			file.open(filename);
			WriteToCSV(W_GK, xc, t, file);	
			file.close();

			filename = "CSV_files/Kolgan2/step_" + std::to_string(step) + ".csv";
			file.open(filename);
			WriteToCSV(W_GK2, xc, t, file);	
			file.close();
			
			filename = "CSV_files/Rodionov/step_" + std::to_string(step) + ".csv";
			file.open(filename);
			WriteToCSV(W_GKR, xc, t, file);	
			file.close();

			filename = "CSV_files/Rodionov2/step_" + std::to_string(step) + ".csv";
			file.open(filename);
			WriteToCSV(W_GKR2, xc, t, file);	
			file.close();
			
			std::string filename = "CSV_files/ENO/step_" + std::to_string(step) + ".csv";
			std::ofstream file(filename);
			WriteToCSV(W_ENO, xc, t, file);	
			file.close();
			*/
			std::string filename = "CSV_files/ActualCalc/WENO/step_" + std::to_string(step) + ".csv";
			file.open(filename);
			WriteToCSV(W_WENO, xc, t, file);	
			file.close();
		}		
		if (t == t_max){
			//std::cout  << "Max time";
			/*std::ofstream file("CSV_files/Godunov.csv");
			WriteToCSV(W_G, xc, t, file);
			file.close();
			
			file.open("CSV_files/Kolgan.csv");
			WriteToCSV(W_GK, xc, t, file);
			file.close();

			file.open("CSV_files/Kolgan2.csv");
			WriteToCSV(W_GK2, xc, t, file);
			file.close();

			file.open("CSV_files/Rodionov.csv");
			WriteToCSV(W_GKR, xc, t, file);
			file.close();

			file.open("CSV_files/Rodionov2.csv");
			WriteToCSV(W_GKR2, xc, t, file);
			file.close();
			*/
			//std::ofstream file("CSV_files/ENO.csv");
			//file.open("CSV_files/ENO.csv");
			//WriteToCSV(W_ENO, xc, t, file);
			//file.close();
			
			std::ofstream file("CSV_files/WENO.csv");
			WriteToCSV(W_WENO, xc, t, file);
			file.close();

			break;			
		}
		/*if (step > 100)
			return 0;*/
		step++;
	}
	std::cout << std::endl << "Завершено успешно." << std::endl;	
	return 0;

}
