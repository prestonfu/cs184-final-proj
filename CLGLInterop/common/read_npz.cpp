#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <common/cnpy.h>
#include "read_npz.h"

#define PATH_MAXX 4096
#define NUM_POINTS 4096
#define NUM_FEATURES 6

// todo: fix code so that images don't get copied unnecessarily

std::vector<std::string> get_files_in_symlink(const std::string& relative_path) {
    char buff[PATH_MAXX];
    getcwd(buff, PATH_MAXX);
    std::string symlink_path(buff);
    symlink_path += relative_path;

    char resolved_path[PATH_MAXX];
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

// filename -> [4096 x] [4096 y] [4096 z] [4096 r] [4096 g] [4096 b]
std::map<std::string, std::vector<float>> read_files_helper(std::vector<std::string> files) {
    std::map<std::string, std::vector<float>> res;
    for (const std::string& file : files) {
        cnpy::NpyArray npy_arr = cnpy::npy_load(file);
        assert(npy_arr.shape.size() == 2 && npy_arr.shape[0] == NUM_POINTS && npy_arr.shape[1] == NUM_FEATURES);
        float* arr = npy_arr.data<float>();
        std::vector<float> vec(3 * NUM_POINTS);
        for (int i = 0; i < NUM_POINTS; i++)
        {
            vec[3 * i] = arr[i];
            vec[3 * i + 1] = arr[i + 2 * NUM_POINTS];
            vec[3 * i + 2] = -arr[i + 1 * NUM_POINTS];
        }
        //std::vector<float> vec(arr, arr + npy_arr.num_vals);

        std::string name = file;
        name = name.substr(name.find_last_of('/') + 1);
        name = name.substr(0, name.length() - 4);
        res[name] = vec;

        std::cout << name << " loaded" << std::endl;
    }
    return res;
}

// filename -> [x0 y0 z0 x1 y1 z1 ... x4095 y4095 z4095]
std::map<std::string, std::vector<float>> read_files(std::string path)
{
    std::vector<std::string> files = get_files_in_symlink(path);
    return read_files_helper(files);
}

// int main() {
//     std::vector<std::string> files = get_files_in_symlink("/../assets/point-cloud-300M");
//     std::map<std::string, std::vector<float>> point_clouds = read_files(files);
//     for (auto it = point_clouds.begin(); it != point_clouds.end(); it++) {
//         std::cout << it->first << " " << it->second[0] << " " << it->second[1] << std::endl;
//     }
//     return 0;
// }

// ln -s /Users/ralphcao/Git/Berkeley/cs-184/final-project/cs184-final-proj/point-e/point_e/examples/images/pc/base300M/npy /Users/ralphcao/Git/Berkeley/cs-184/final-project/cs184-final-proj/CLGLInterop/Assets/point-cloud-300M