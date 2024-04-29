#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <common/cnpy.h>

#define PATH_MAX 4096
#define NUM_POINTS 4096
#define NUM_FEATURES 6

// todo: fix code so that images don't get copied unnecessarily

std::vector<std::string> get_files_in_symlink(const std::string& relative_path) {
    char buff[PATH_MAX];
    getcwd(buff, PATH_MAX);
    std::string symlink_path(buff);
    symlink_path += relative_path;

    char resolved_path[PATH_MAX];
    if (realpath(symlink_path.c_str(), resolved_path) == nullptr) {
        throw std::runtime_error("Error resolving symbolic link: " + symlink_path);
    }

    std::vector<std::string> files;
    DIR* dir = opendir(resolved_path);
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string filename(entry->d_name);
            if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".npy") {
                files.push_back(symlink_path + "/" + filename);
            }
        }
        closedir(dir);
    }
    return files;
}

std::map<std::string, std::vector<double>> read_files(std::vector<std::string> files) {
    std::map<std::string, std::vector<double>> res;
    for (const std::string& file : files) {
        cnpy::NpyArray npy_arr = cnpy::npy_load(file);
        assert(npy_arr.shape.size() == 2 && arr.shape[0] == NUM_POINTS && arr.shape[1] == NUM_FEATURES);
        double* arr = npy_arr.data<double>();
        std::vector<double> vec(arr, arr + npy_arr.num_vals);
        std::cout << "debug " << *(arr + 0) << " " << vec[0] << " " << vec[1] << " " << npy_arr.num_vals << " " << std::endl;
        res[file] = vec;
    }
    return res;
}

int main() {
    std::vector<std::string> files = get_files_in_symlink("/../assets/point-cloud-300M");
    std::map<std::string, std::vector<double>> point_clouds = read_files(files);
    for (auto it = point_clouds.begin(); it != point_clouds.end(); it++) {
        std::cout << it->first << " " << it->second[0] << " " << it->second[1] << std::endl;
    }
    return 0;
}