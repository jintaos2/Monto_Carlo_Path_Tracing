#ifndef _OBJ_H_
#define _OBJ_H_

#include <vector>
#include <map>
#include <regex>
#include <glm.hpp>

#include "bmp.h"
#include "hdr.h"
#include "tools.h"

class Picture
{
public:
    bool empty_ = true;
    std::string filename_;
    int w, h;
    std::vector<glm::dvec3> data_;
    Picture() {}
    void load(std::string filename)
    {
        filename_ = filename;
        int i = filename.size() - 3;
        if (filename.substr(i, 3) == "bmp")
            load_bmp(filename);
        else if (filename.substr(i, 3) == "hdr")
            load_hdr(filename);
    }
    void load_bmp(std::string filename)
    {
        Bitmap bmp(filename.c_str());
        if (bmp.GetBits() != NULL)
        {
            w = bmp.GetW();
            h = bmp.GetH();
            for (int i = 0; i < w * h; ++i)
            {
                glm::dvec3 color = bmp.GetColor(i);
                data_.push_back(color);
            }
            empty_ = false;
        }
    }
    void load_hdr(std::string filename)
    {
        HDRLoader hdr(filename.c_str());
        if (hdr.success)
        {
            std::cout << filename << std::endl;
            w = hdr.color.width;
            h = hdr.color.height;
            for (int i = 0; i < w * h; ++i)
            {
                glm::dvec3 color = hdr.GetColor(i);
                data_.push_back(color);
            }
            empty_ = false;
        }
    }
    inline glm::dvec3 Sample2D(glm::dvec2 uv)
    {
        int x1 = uv.x * w;
        int y1 = uv.y * h;
        x1 = between(0, w, x1);
        y1 = between(0, h, y1);
        return data_[y1 * w + x1];
    }
};

struct Material
{
    glm::dvec3 Kd = glm::dvec3(1, 1, 1); // diffuse, 反射光线系数，0 表示吸收所有光线.
    glm::dvec3 Ks = glm::dvec3(0, 0, 0); // specular, 高光反射系数.
    glm::dvec3 Kr = glm::dvec3(0, 0, 0); // 折射透明度, 0 表示不透明.
    glm::dvec3 Le = glm::dvec3(0, 0, 0); // 自发光.
    double kd = 1;                       // diffuse fraction, Monte Carlo
    double ks = 0;                       // specular fraction, Monte Carlo
    double kr = 0;                       // reflect fraction, Monte Carlo
    double Nr = 1;                       // 物质折射率..
    double Ns = 1;                       // Phong 高光反射参数.
    Picture Map_Kd;                      // 纹理贴图.
};

class Model
{
public:
    std::vector<glm::dvec3> vertex_;
    std::vector<glm::dvec3> norm_;
    std::vector<glm::dvec2> uv_coord_;
    std::vector<Material> material_; // material_[0] is default, no matter .mtl file exists or not
    std::vector<glm::imat3x4> face_; // 3 line of {vertex_idx, norm_idx, uv_coord_idx, material_idx}
                                     // idx start from 0

    std::map<std::string, int> material_name; // find material idx by name
    std::string dir_path;                     // 根目录..

    ~Model() {}

    Model(std::string filename)
    {
        std::smatch result;
        // get the path of obj file
        if (regex_search(filename, result, std::regex("(.*/)[^/]+")) && result.size() == 2)
            dir_path = result[1];
        else
            dir_path = "";
        std::ifstream fs;
        fs.open(filename);
        if (!fs.is_open())
        {
            std::cout << "[error] open " << filename << " failed!\n";
            return;
        }
        std::string str((std::istreambuf_iterator<char>(fs)), std::istreambuf_iterator<char>());
        std::vector<std::string> lines;
        int start = 0;
        for (int i = 0; i < str.size(); ++i)
        {
            if (str[i] == '\n')
            {
                lines.push_back(str.substr(start, i - start + 1));
                start = i + 1;
            }
        }

        // the first material is always the default material;
        material_.push_back(Material());

        std::vector<std::string> v_;
        std::vector<std::string> vn_;
        std::vector<std::string> vt_;
        std::vector<std::string> f_;
        for (auto line : lines)
        {
            if (line[0] == 'm') // mtllib
            {
                std::smatch res;
                if (regex_search(line, res, std::regex("mtllib\\s+(\\S+)")))
                    load_material(res[1]);
            }
            else if (line[0] == 'v')
            {
                if (line[1] == ' ') // v
                    v_.push_back(line);
                else if (line[1] == 'n') // vn
                    vn_.push_back(line);
                else if (line[1] == 't') // vt
                    vt_.push_back(line);
            }
            else if (line[0] == 'f' && line[1] == ' ') // f
            {
                f_.push_back(line);
            }
            else if (line[0] == 'u' && line[1] == 's') // usemtl
            {
                f_.push_back(line);
            }
        }

#pragma omp sections
        {
#pragma omp section
            {
                for (auto line : v_)
                {
                    char trash;
                    glm::dvec3 temp;
                    std::istringstream ss(line);
                    ss >> trash >> temp[0] >> temp[1] >> temp[2];
                    vertex_.push_back(temp);
                }
            }
#pragma omp section
            {
                for (auto line : vn_)
                {

                    char trash;
                    glm::dvec3 temp;
                    std::istringstream ss(line);
                    ss >> trash >> trash >> temp[0] >> temp[1] >> temp[2];
                    norm_.push_back(temp);
                }
            }
#pragma omp section
            {
                for (auto line : vt_)
                {
                    char trash;
                    glm::dvec2 temp;
                    std::istringstream ss(line);
                    ss >> trash >> trash >> temp[0] >> temp[1];
                    uv_coord_.push_back(temp);
                }
            }
#pragma omp section
            {
                int curr_material_id = 0;
                for (auto line : f_)
                {
                    if (line[0] == 'u')
                    {
                        std::smatch res;
                        if (regex_search(line, res, std::regex("usemtl\\s+(\\S+)")))
                            curr_material_id = material_name[std::string(res[1])];
                    }
                    else
                    {
                        char trash;
                        glm::imat3x4 face;
                        std::istringstream ss(line);
                        ss >> trash;
                        for (int i = 0; i < 3; ++i)
                        {
                            ss >> face[i][0] >> trash >> face[i][1] >> trash >> face[i][2];
                            face[i][0]--;
                            face[i][1]--;
                            face[i][2]--;
                            face[i][3] = curr_material_id;
                        }
                        face_.push_back(face);
                    }
                }
            }
        }
        std::cout << "read " << filename << std::endl;
    }

    inline void load_material(std::string filename)
    {
        std::ifstream fs;
        fs.open(dir_path + filename);
        if (!fs.is_open())
        {
            std::cout << "[error] open " << filename << " failed!\n";
            return;
        }
        std::string line;
        std::smatch result;
        // read .mtl file
        while (getline(fs, line))
        {
            if (regex_search(line, std::regex("\\s*#")) || regex_match(line, std::regex("\\s*")))
                continue;
            // 'newmtl' find, add a default material to array
            if (regex_search(line, result, std::regex("newmtl\\s+(\\S+)")) && result.size() == 2)
            {
                material_name[std::string(result[1])] = material_.size();
                material_.push_back(Material());
                continue;
            }
            // update material details
            if (regex_search(line, result, std::regex("\\s*Kd\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Kd = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, std::regex("\\s*Ks\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Ks = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, std::regex("\\s*Kr\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Kr = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, std::regex("\\s*Le\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Le = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, std::regex("\\s*Ns\\s+(\\S+)")))
                material_.back().Ns = tofloat(result[1]);
            else if (regex_search(line, result, std::regex("\\s*Nr\\s+(\\S+)")))
                material_.back().Nr = tofloat(result[1]);
            else if (regex_search(line, result, std::regex("map_Kd\\s+(\\S+)")))
                material_.back().Map_Kd.load(dir_path + std::string(result[1])); // load picture
            else
                continue;
        }
        for (Material &i : material_)
        {
            double kd = max3(i.Kd.x, i.Kd.y, i.Kd.z);
            double ks = max3(i.Ks.x, i.Ks.y, i.Ks.z);
            double kr = max3(i.Kr.x, i.Kr.y, i.Kr.z);
            double _sum = kd + ks + kr; // kd + ks + kr <= 1;
            if (_sum > 1)
            {
                kd /= _sum;
                ks /= _sum;
                kr /= _sum;
                i.Kd /= _sum;
                i.Ks /= _sum;
                i.Kr /= _sum;
            }
            // 调整加权系数，使得 Kd * kd + Ks * ks + Kr * kr = (1,1,1)
            if (kd > 0)
                i.Kd /= kd;
            if (ks > 0)
                i.Ks /= ks;
            if (kr > 0)
                i.Kr /= kr;
            i.kd = kd;
            i.ks = ks;
            i.kr = kr;
        }
    }
};

#endif
