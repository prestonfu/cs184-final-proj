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

    // char resolved_path[PATH_MAX];

    // if (realpath(symlink_path.c_str(), resolved_path) == nullptr) {
    //     throw std::runtime_error("Error resolving symbolic link: " + symlink_path);
    // }
    // std::cout<< "hello5" << std::endl;

    const char* resolved_path = symlink_path.c_str();
    std::vector<std::string> files;
    DIR* dir = opendir(resolved_path);
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string filename(entry->d_name);
            if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".npz") {
                files.push_back(symlink_path + "/" + filename);
            }
        }
        closedir(dir);
    }
    return files;
}

std::map<std::string, double*> read_files(std::vector<std::string> files) {
    std::map<std::string, double*> res;
    for (const std::string& file : files) {
        cnpy::npz_t npz = cnpy::npz_load(file);
        int cnt = 0;
        // std::cout<<file<<std::endl;
        for (const auto& pair : npz) {
            std::cout<<"here2?"<<std::endl;
            if (cnt >= 1) {
                throw std::runtime_error("Something fishy here");
            }
            std::cout << pair.first << std::endl;
            cnpy::NpyArray npy_arr = pair.second;
            assert(npy_arr.shape.size() == 2 && arr.shape[0] == NUM_POINTS && arr.shape[1] == NUM_FEATURES);
            double* arr = npy_arr.data<double>();
            res[file.substr(0, file.size() - 4)] = arr;
            std::cout << "here\n";
            cnt++;
        }
    }
    return res;
}

int main() {
    std::vector<std::string> files = get_files_in_symlink("/../assets/points");
    std::map<std::string, double*> point_clouds = read_files(files);
    for (auto it = point_clouds.begin(); it != point_clouds.end(); it++) {
        std::cout << it->first << " " << it->second[0];
    }

    return 0;
}