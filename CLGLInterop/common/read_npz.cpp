#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include "cnpy.h"

std::vector<std::string> get_files_in_dir(const std::string& dirpath) {
    std::vector<std::string> files;
    DIR* dir = opendir(dirpath.c_str());
    if (dir) {
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            std::string filename(entry->d_name);
            if (filename != "." && filename != "..") {
                files.push_back(dirpath + "/" + filename);
            }
        }
        closedir(dir);
    }
    return files;
}

int main() {
    // Get the path to the symbolic link
    char buff[PATH_MAX];
    getcwd(buff, PATH_MAX);
    std::string symlink_path(buff);
    symlink_path += "../assets/point_clouds";

    // Resolve the symbolic link
    char resolved_path[PATH_MAX];
    if (realpath(symlink_path.c_str(), resolved_path) == nullptr) {
        std::cerr << "Error resolving symbolic link: " << symlink_path << std::endl;
        return 1;
    }

    // Get the list of files in the resolved directory
    std::vector<std::string> files = get_files_in_dir(resolved_path);

    // Process each .npz file
    for (const std::string& file : files) {
        if (file.substr(file.size() - 4) == ".npz") {
            std::cout << "Processing file: " << file << std::endl;
            cnpy::npz_t npz = cnpy::npz_load(file);
            for (const auto& pair : npz.data) {
                std::cout << "Array name: " << pair.first << std::endl;
                cnpy::NdArrayPtr arr = pair.second;
                // Do something with the NumPy array
            }
        }
    }

    return 0;
}