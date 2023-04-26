//
// Created by admin on 2023/4/25.
//

#include <fstream>
#include <string>
#include "function.h"

void ParseConfigureLine(const std::string_view& line, std::string& key, std::string& value)
{
    if (line.empty()) return;
    bool is_key{true};
    for (auto word : line) {
        if (word == ' ') continue;
        if (word == '=') { is_key = false; continue; }
        if (word == '#') break;
        is_key ? key += word : value += word;
    }
}

void IdentifyConfigure(const std::string& key, const std::string& value, Config& config)
{
    if (key == "posmode") config.pos_mode = stoi(value);
    else if (key == "online") config.online = value == "on";
    else if (key == "elmask") config.el_mask = stod(value);
    else if (key == "snrmask_r") config.snr_mask_r = std::stod(value);
    else if (key == "snrmask_b") config.snr_mask_b = std::stod(value);
    else if (key == "ionoopt") config.iono_opt = std::stoi(value);
    else if (key == "tropopt") config.trop_opt = std::stoi(value);
    else if (key == "ratiothres") config.ratio_thres = std::stod(value);
    else if (key == "solformat") config.sol_format = std::stoi(value);
    else if (key == "outmode") config.out_mode = stoi(value);
    else if (key == "outhead") config.out_head = value == "on";
    else if (key == "outvel") config.out_vel = value == "on";
    else if (key == "basepostype") config.basepos_type = stoi(value);
    else if (key == "basepos1") config.basepos[0] = stod(value);
    else if (key == "basepos2") config.basepos[1] = stod(value);
    else if (key == "basepos3") config.basepos[2] = stod(value);
    else if (key == "infile_r") config.infile_r = value;
    else if (key == "infile_b") config.infile_b = value;
    else if (key == "IP_r") config.IP_r = value;
    else if (key == "IP_b") config.IP_b = value;
    else if (key == "port_r") config.port_r = stoi(value);
    else if (key == "port_b") config.port_b = stoi(value);
    else return;
}

bool ReadConfigureFile(const char* filename, Config& config)
{
    std::ifstream file{filename, std::ios::in};
    std::string line, key, value;
    if (!file.is_open()) return false;
    while (getline(file, line)) {
        if (line.substr(0, 1) == "#") continue;// header of file
        // get key, value of every line
        ParseConfigureLine(line, key, value);
        if (!key.empty() && !value.empty()) IdentifyConfigure(key, value, config);
        key.clear();
        value.clear();
    }
    file.close();
}
