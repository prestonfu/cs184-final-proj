#include <iostream>
#include <vector>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <map>

std::map<std::string, std::pair<std::vector<float>, std::vector<float>>> read_files(std::string path);