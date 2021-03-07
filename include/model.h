#ifndef _OBJ_H_
#define _OBJ_H_

#include <vector>
#include <map>
#include <regex>
#include <glm.hpp>
#include "bmp.h"

#include "hdr.h"

#include "tools.h"

// TODO: load and store jpg file
class Picture
{
public:
    bool empty_ = true;
    string filename_;
    int w, h;
    vector<vec3> data_;
    void load(string filename)
    {
        filename_ = filename;
        int i = filename.size() - 3;
        if (filename.substr(i, 3) == "bmp")
            load_bmp(filename);
        else if (filename.substr(i, 3) == "hdr")
            load_hdr(filename);
    }
    void load_bmp(string filename)
    {
        Bitmap bmp(filename.c_str());
        if (bmp.GetBits() != NULL)
        {
            w = bmp.GetW();
            h = bmp.GetH();
            for (int i = 0; i < w * h; ++i)
            {
                vec3 color = bmp.GetColor(i);
                data_.push_back(color);
                // cout << print_vec3(color) << endl;
            }
            empty_ = false;
        }
    }
    void load_hdr(string filename)
    {
        HDRLoader hdr(filename.c_str());
        if (hdr.success)
        {
            cout << filename << endl;
            w = hdr.color.width;
            h = hdr.color.height;
            for (int i = 0; i < w * h; ++i)
            {
                vec3 color = hdr.GetColor(i);
                data_.push_back(color);
            }
            empty_ = false;
        }
    }
    inline vec3 Sample2D(vec2 uv)
    {
        int x1 = uv.x * w;
        int y1 = uv.y * h;
        x1 = between(0, w, x1);
        y1 = between(0, h, y1);
        return data_[y1 * w + x1];
    }
};

// TODO: default value
struct Material
{
    vec3 Kd = vec3(1, 1, 1); // diffuse, 反射光线系数，0 表示吸收所有光线.
    vec3 Ks = vec3(0, 0, 0); // specular, 高光反射系数.
    vec3 Kr = vec3(0, 0, 0); // 折射透明度, 0 表示不透明.
    vec3 Le = vec3(0, 0, 0); // 自发光.
    float kd = 1;            // diffuse fraction, Monte Carlo
    float ks = 0;            // specular fraction, Monte Carlo
    float kr = 0;            // reflect fraction, Monte Carlo
    float Nr = 1;            // 物质折射率..
    float Ns = 1;            // Phong 高光反射参数.
    Picture Map_Kd;          // 纹理贴图.
};

class Model
{
public:
    vector<vec3> vertex_;
    vector<vec3> norm_;
    vector<vec2> uv_coord_;
    vector<Material> material_; // material_[0] is default, no matter .mtl file exists or not
    vector<imat3x4> face_;      // 3 line of {vertex_idx, norm_idx, uv_coord_idx, material_idx}
                                // idx start from 0

    map<string, int> material_name; // find material idx by name
    string dir_path;                // 根目录..

    ~Model() {}

    Model(string filename)
    {
        smatch result;
        // get the path of obj file
        if (regex_search(filename, result, regex("(.*/)[^/]+")) && result.size() == 2)
            dir_path = result[1];
        else
            dir_path = "";
        ifstream fs;
        fs.open(filename);
        if (!fs.is_open())
        {
            cout << "[error] open " << filename << " failed!\n";
            return;
        }
        string str((istreambuf_iterator<char>(fs)), istreambuf_iterator<char>());
        vector<string> lines;
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

        vector<string> v_;
        vector<string> vn_;
        vector<string> vt_;
        vector<string> f_;
        for (auto line : lines)
        {
            if (line[0] == 'm') // mtllib
            {
                smatch res;
                if (regex_search(line, res, regex("mtllib\\s+(\\S+)")))
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
                    vec3 temp;
                    istringstream ss(line);
                    ss >> trash >> temp[0] >> temp[1] >> temp[2];
                    // cout << "v  " << print_vec3(temp) << endl;
                    vertex_.push_back(temp);
                }
            }
#pragma omp section
            {
                for (auto line : vn_)
                {

                    char trash;
                    vec3 temp;
                    istringstream ss(line);
                    ss >> trash >> trash >> temp[0] >> temp[1] >> temp[2];
                    // cout << "vn " << print_vec3(temp) << endl;
                    norm_.push_back(temp);
                }
            }
#pragma omp section
            {
                for (auto line : vt_)
                {
                    char trash;
                    vec2 temp;
                    istringstream ss(line);
                    ss >> trash >> trash >> temp[0] >> temp[1];
                    // cout << "vt " << temp[0] << ' ' << temp[1] << endl;
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
                        smatch res;
                        if (regex_search(line, res, regex("usemtl\\s+(\\S+)")))
                            curr_material_id = material_name[string(res[1])];
                    }
                    else
                    {
                        char trash;
                        imat3x4 face;
                        istringstream ss(line);
                        ss >> trash;
                        for (int i = 0; i < 3; ++i)
                        {
                            ss >> face[i][0] >> trash >> face[i][1] >> trash >> face[i][2];
                            face[i][0]--;
                            face[i][1]--;
                            face[i][2]--;
                            face[i][3] = curr_material_id;
                        }
                        // cout << "face:\n" << print_mat3x4(face) << endl;
                        face_.push_back(face);
                    }
                }
            }
        }
        cout << "read done" << endl;
    }

    inline void load_material(string filename)
    {
        ifstream fs;
        fs.open(dir_path + filename);
        if (!fs.is_open())
        {
            cout << "[error] open " << filename << " failed!\n";
            return;
        }
        string line;
        smatch result;
        // read .mtl file
        while (getline(fs, line))
        {
            if (regex_search(line, regex("\\s*#")) || regex_match(line, regex("\\s*")))
                continue;
            // newmtl find, add a default material to array
            if (regex_search(line, result, regex("newmtl\\s+(\\S+)")) && result.size() == 2)
            {
                material_name[string(result[1])] = material_.size();
                material_.push_back(Material());
                continue;
            }
            // update material details
            if (regex_search(line, result, regex("\\s*Kd\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Kd = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, regex("\\s*Ks\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Ks = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, regex("\\s*Kr\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Kr = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, regex("\\s*Le\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")))
                material_.back().Le = {tofloat(result[1]), tofloat(result[2]), tofloat(result[3])};
            else if (regex_search(line, result, regex("\\s*Ns\\s+(\\S+)")))
                material_.back().Ns = tofloat(result[1]);
            else if (regex_search(line, result, regex("\\s*Nr\\s+(\\S+)")))
                material_.back().Nr = tofloat(result[1]);
            else if (regex_search(line, result, regex("map_Kd\\s+(\\S+)")))
                material_.back().Map_Kd.load(dir_path + string(result[1])); // load picture
            else
                continue;
        }
        for (Material &i : material_)
        {
            float kd = max3(i.Kd.x, i.Kd.y, i.Kd.z);
            float ks = max3(i.Ks.x, i.Ks.y, i.Ks.z);
            float kr = max3(i.Kr.x, i.Kr.y, i.Kr.z);
            float _sum = kd + ks + kr; // kd + ks + kr <= 1;
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
