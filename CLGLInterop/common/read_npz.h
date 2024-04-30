#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <map>

std::map<std::string, std::vector<float>> read_files(std::string path);