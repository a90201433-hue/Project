#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "Types.h"

extern int Nx, Ny, fict;
extern double gamm;


/*double ComputeMinInternalEnergy(const Field& W) {

	size_t Nx_cells = Nx - 1;
    size_t Ny_cells = Ny - 1;

    double e_min = std::numeric_limits<double>::max();

    for (size_t i = fict; i < Nx_cells + fict; ++i)
    {
        for (size_t j = fict; j < Ny_cells + fict; ++j)
        {
            double rho = W[i][j][0];
            double P   = W[i][j][NEQ - 1];

            double e = P / (rho * (gamm - 1.0));

            if (e < e_min)
                e_min = e;
        }
    }
    return e_min;

}*/

void SaveFieldToCSV(const Field& W,
                    const std::vector<double>& x,
                    const std::vector<double>& y,
                    const double& time,
                    const std::string& filename) {

    size_t Nx_cells = Nx - 1;
    size_t Ny_cells = Ny - 1;

    std::ofstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Невозможно открыть файл" << filename << std::endl;
        return;
    }

    file << "t,x,y,rho,u,v,P,e\n";

    for (size_t i = fict; i < Nx_cells + fict; i++) {

        double xc = 0.5 * (x[i] + x[i + 1]);

        for (size_t j = fict; j < Ny_cells + fict; j++) {
            double yc = 0.5 * (y[j] + y[j + 1]);

			double rho = W[i][j][0];
            double P   = W[i][j][NEQ - 1];
            double e = P / (rho * (gamm - 1.0));

            file << time << ","
                 << xc << ","
                 << yc << ","
                 << W[i][j][0] << "," 	/* rho */
                 << W[i][j][1] << "," 	/* u */
                 << W[i][j][2] << "," 	/* v */
                 << W[i][j][3] << ","	/* P */
				 << e << "\n"; 			/* e */
        }
    }

    file.close();
}

fs::path CreateDirFromPath(const std::string& file_path) {
    try {
        fs::path input_path(file_path);

        // Убираем расширение (.toml)
        fs::path dir_path = input_path.parent_path() /
                            input_path.stem();   // stem = имя без расширения

        // Создаём директорию
        if (!fs::exists(dir_path)) {
            fs::create_directories(dir_path);
        }
		return dir_path;
    }

    catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: "
                  << e.what() << std::endl;
		return fs::path();
    }

	
}

/*
void ClearDirectoryContents(const std::string& target_dir_path) {
    fs::path target_path = target_dir_path;
    

    try {
        // 1. Проверка существования и типа пути
        if (!fs::exists(target_path)) {
            std::cerr << "Ошибка: Папка не существует: " << target_dir_path << std::endl;
            return;
        }
        if (!fs::is_directory(target_path)) {
            std::cerr << "Ошибка: Путь не является папкой: " << target_dir_path << std::endl;
            return;
        }

        // 2. Итерация по всем элементам внутри папки
        std::cout << "Начало очистки содержимого папки: " << target_dir_path << std::endl;
        
        // directory_iterator проходит по всем элементам первого уровня
        for (const auto& entry : fs::directory_iterator(target_path)) {
            try {
                // fs::remove_all() рекурсивно удаляет текущий файл или подпапку
                fs::remove_all(entry.path());
                //std::cout << "   - Удалено: " << entry.path().string() << std::endl;
            } catch (const fs::filesystem_error& e) {
                // Обработка ошибки удаления конкретного элемента
                std::cerr << "   - Ошибка удаления: " << entry.path().string() << " | " << e.what() << std::endl;
                
            }
        }
        
        std::cout << "Очистка завершена." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Критическая ошибка при работе с директорией: " << e.what() << std::endl;
        return;
    }

    return;
}

void MoveToChache(const std::string& source_root_path, const std::string& cache_root_path) {

	fs::path source_path = source_root_path;
    	fs::path cache_path = cache_root_path;
    	

    	try {
        	// 1. Проверка существования исходной папки
        	if (!fs::exists(source_path) || !fs::is_directory(source_path)) {
            	std::cerr << "Ошибка: Исходная папка не существует или не является папкой: " << source_root_path << std::endl;
            	return;
        	}

        	// 2. Убедиться, что папка Кэша существует. Если нет, создать ее.
        	fs::create_directories(cache_path);
		ClearDirectoryContents(cache_root_path);

  	
		

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
                		//std::cout << "Перемещено: " << current_item.string() << " -> " << destination_path.string() << std::endl;

            		} catch (const fs::filesystem_error& e) {
                	// Если возникает ошибка при перемещении одного элемента
                		std::cerr << "Ошибка перемещения " << current_item.string() << ": " << e.what() << std::endl;
                		
            		}
        	}
		std::cout << "Данные прошлых расчетов перемещены в папку " << cache_root_path << "\n" << std::endl;

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

void CreateDir(const std::string& root_folder_path, const std::string& new_folder_name) {
    
    // 1. Объединение путей в единый целевой путь
    // Оператор деления (/) в fs::path корректно добавляет разделитель (/, \)
    fs::path final_path = fs::path(root_folder_path) / new_folder_name;

    try {
        // 2. Создание всех необходимых папок (включая корневую, если ее нет)
        if (fs::create_directories(final_path)) {
            
        } else {
            // Проверка, существовала ли папка уже
            if (fs::exists(final_path) && fs::is_directory(final_path)) {
                
            } else {       
                return;
            }
        }
        return;

    } catch (const fs::filesystem_error& e) {
        // Ловит ошибки (права доступа, недопустимый путь и т.д.)
        return;
    }
}

void WriteToCSV(std::vector<std::vector<double>> W, std::vector<double> xc, double t, std::ofstream& file){
	file << "time,x,rho,u,P,e" << std::endl;
	std::vector<double> e1(Nx + 2 * fict - 1);
	double e;
	for (int i = 0; i < Nx + 2 * fict - 1; i++) {
		e1[i] = W[i][2] / (W[i][0] * (gamm - 1));
	}
	double e0 = *std::min_element(e1.begin(), e1.end());
	
	for (int i = 0; i < Nx + 2 * fict - 1; i++){
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

// Функция для подсчета ошибок и записи в CSV
void SaveAnalysisData(
    std::ofstream& file,
    double t,
    const std::vector<std::vector<double>>& W_num,
    const std::vector<double>& xc,                 
    // Параметры для точного решения:
    const std::vector<double>& W_L,
    const std::vector<double>& W_R,
    const std::vector<double>& W_star,
    double x0,
    std::vector<double>(*TrueSolveF)(std::vector<double>, std::vector<double>, std::vector<double>, double, double),
    double(*AnalisF)(std::vector<double>)
)
 {
   
    std::vector<double> error_diff_rho(Nx + 2 * fict), error_diff_v(Nx + 2 * fict), error_diff_P(Nx + 2 * fict);
    

    for (int i = fict; i < Nx + fict; i++) {
        // 1.точное решение
        std::vector<double> W_exact = TrueSolveF(W_L, W_R, W_star, xc[i] - x0, t);
        
        // 2. разница
        error_diff_rho[i] = W_num[i][0] - W_exact[0];
	error_diff_v[i] = W_num[i][1] - W_exact[1];
       	error_diff_P[i] = W_num[i][2] - W_exact[1];	
    }

   
    double val_rho = AnalisF(error_diff_rho);
    double val_v = AnalisF(error_diff_v);
    double val_P = AnalisF(error_diff_P);

   
    file << t << "," << val_rho << "," << val_v << "," << val_P << std::endl;
}*/

