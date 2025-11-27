#ifndef _FILE_PROCESSING_H_
#define _FILE_PROCESSING_H_

void MoveToChache(const std::string& source_root_path, const std::string& cache_root_path);

void CleanDir(std::vector<std::string> folders);

void WriteToCSV(std::vector<std::vector<double>> W, std::vector<double> xc, double t, std::ofstream& file);

void ClearDirectoryContents(const std::string& target_dir_path);

void CreateDir(const std::string& root_folder_path, const std::string& new_folder_name);

#endif
