#ifndef _RAY_TRACING_H_
#define _RAY_TRACING_H_

#include <glm.hpp>
#include <ctime>
#include <cmath>
#include <random>
#include "omp.h"
#include "model.h"

#define NS 3
#define RAY_DEEPTH 9

// 3D model + pose and scale
struct Obj
{
    Model &model;
    glm::dmat4x4 pose = glm::dmat4x4(1.0);
    double scale = 1.0;

    Obj(Model &m_) : model(m_) {}

    void set_pose(glm::dmat4x4 pose_, double s_)
    {
        pose = pose_;
        scale = s_;
    }
};

class FrameBuffer
{
public:
    uint32_t *fb_;    // RGBA from memory low to high , y: ^ x: ->
    glm::dvec3 *fb2_; // hdr rgb vector
    int w_, h_;       // screen width, height
    FrameBuffer(int w, int h) : w_(w), h_(h)
    {
        fb_ = new uint32_t[w * h];
        fb2_ = new glm::dvec3[w * h];
    }
    ~FrameBuffer()
    {
        if (fb_)
        {
            delete[] fb_;
            fb_ = NULL;
        }
        if (fb2_)
        {
            delete[] fb2_;
            fb2_ = NULL;
        }
    }
    inline void clear()
    {
        for (int i = 0; i < w_ * h_; ++i)
            fb_[i] = 0;
        for (int i = 0; i < w_ * h_; ++i)
            fb2_[i] = glm::dvec3(0, 0, 0);
    }
    inline void set_pixel(int x, int y, uint32_t color) { fb_[y * w_ + x] = color; }
    // 累加光照强度..
    inline void accumulate(int x, int y, glm::dvec3 color) { fb2_[y * w_ + x] += color; }
    // 平均光照强度..
    inline void show(double n, double adapted_lum)
    {
        for (int i = 0; i < w_ * h_; ++i)
            fb_[i] = toRGB(ToneMapping(fb2_[i] / n, adapted_lum));
    }
};

struct BRDF
{
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dis1;
    std::uniform_real_distribution<double> dis2;

    BRDF() : gen(rd()), dis1(0, 1), dis2(-1, 1) {}

    inline double randf1() { return dis1(gen); }
    inline double randf2() { return dis2(gen); }
    inline glm::dvec3 rand_sphere()
    {
        glm::dvec3 ret;
        do
        {
            ret.x = randf2();
            ret.y = randf2();
            ret.z = randf2();
        } while (glm::dot(ret, ret) > 1.0);
        return glm::normalize(ret);
    }
    // 法向量坐标系转世界坐标系.. input: unit vector
    inline glm::dmat3x3 norm_to_world(glm::dvec3 n)
    {
        if (n.z > 0.999999)
            return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        if (n.z < -0.999999)
            return {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
        double cos_theta0 = n.z;
        double sin_theta0 = glm::sqrt(1.0 - cos_theta0 * cos_theta0);
        double cos_phi0 = n.x / sin_theta0;
        double sin_phi0 = n.y / sin_theta0;

        return {{cos_theta0 * cos_phi0, -sin_phi0, n.x},
                {cos_theta0 * sin_phi0, cos_phi0, n.y},
                {-sin_theta0, 0.0, n.z}};
    }
    // phong 散射 重要性采样..
    inline glm::dvec3 diffuse(glm::dvec3 &n)
    {
        double u1 = randf1();
        double sin_theta = glm::sqrt(1.0 - u1);
        double cos_theta = glm::sqrt(u1);

        double phi = 2.0 * std::_Pi * randf1();
        double sin_phi = glm::sin(phi);
        double cos_phi = glm::cos(phi);
        glm::dvec3 out(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta + 1e-20);
        return glm::normalize(out * norm_to_world(n));
    }
    // phong 高光反射 重要性采样.
    inline glm::dvec3 specular(glm::dvec3 &ref, double Ns)
    {
        double cos_theta = std::pow(randf1(), 1.0 / (Ns + 1.0));
        double sin_theta = glm::sqrt(1.0 - cos_theta * cos_theta);

        double phi = 2.0 * std::_Pi * randf1();
        double sin_phi = glm::sin(phi);
        double cos_phi = glm::cos(phi);
        glm::dvec3 out(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta + 1e-20);
        return glm::normalize(out * norm_to_world(ref));
    }
    inline glm::dvec3 refract(glm::dvec3 &in, glm::dvec3 &n, double eta)
    {
        double N_dot_I = glm::dot(n, in);
        double k = 1.0 - eta * eta * (1.0 - N_dot_I * N_dot_I);
        if (k < 0)
            return glm::reflect(in, n);
        else
            return glm::normalize(eta * in - (eta * N_dot_I + glm::sqrt(k)) * n);
    }
};

class BVH;
// vertex + material
struct Face3D
{
    glm::dvec3 v1, v2, v3;    // 3 vertexes
    glm::dvec3 n1, n2, n3;    // 3 norms
    glm::dvec2 uv1, uv2, uv3; // 3 uv_coordinates
    glm::dvec3 A, B;          // bounding box
    Material &m;

    inline void print_vertex()
    {
        std::cout << print_vec3(v1) << print_vec3(v2) << print_vec3(v3) << std::endl;
    }
};
struct HitInfo
{
    BVH *box;
    double t1, t2;
};
struct Ray
{
    glm::dvec3 p;     // start point
    glm::dvec3 d;     // direction
    glm::dvec3 color; // rgb
    int deepth = 0;   // 反射次数..
    Face3D *hit_face = NULL;
    double hit_time = DBL_MAX;
    glm::dvec3 barycentric; // 重心坐标系 系数..
    std::vector<HitInfo> hit_results;
};

class BVH
{
public:
    BVH *left = NULL;
    BVH *right = NULL;                                    // child node
    glm::dvec3 A = glm::dvec3(DBL_MAX, DBL_MAX, DBL_MAX); // bounding box
    glm::dvec3 B = glm::dvec3(-DBL_MAX, -DBL_MAX, -DBL_MAX);
    std::vector<Face3D *> faces;

    BVH() {}

    inline void hit_test(Ray &ray)
    {
        glm::dvec3 temp = ray.d;
        if (ray.d.x == 0)
            ray.d.x = 1e-20;
        if (ray.d.y == 0)
            ray.d.y = 1e-20;
        if (ray.d.z == 0)
            ray.d.z = 1e-20;
        hit_test_box_(ray);
        ray.d = temp;
        std::vector<HitInfo> &result = ray.hit_results;
        if (result.empty())
            return;
        std::sort(std::begin(result), std::end(result), [](HitInfo &a, HitInfo &b) -> bool { return a.t1 < b.t1; });
        for (auto i : result)
        {
            if (i.t1 >= ray.hit_time)
                break;
            hit_test_faces(ray, i.box->faces);
        }
        result.clear();
        return;
    }
    inline void hit_test_box_(Ray &ray)
    {
        glm::dvec3 t1 = (A - ray.p) / ray.d;
        glm::dvec3 t2 = (B - ray.p) / ray.d;
        if (t1.x > t2.x)
            std::swap(t1.x, t2.x);
        if (t1.y > t2.y)
            std::swap(t1.y, t2.y);
        if (t1.z > t2.z)
            std::swap(t1.z, t2.z);
        double box_t1 = max3(t1.x, t1.y, t1.z);
        double box_t2 = min3(t2.x, t2.y, t2.z);
        if (box_t1 <= box_t2 && box_t2 > 0)
        {
            if (left)
                left->hit_test_box_(ray);
            if (right)
                right->hit_test_box_(ray);
            if (!faces.empty()) // leaf node
                ray.hit_results.push_back({this, box_t1, box_t2});
        }
        return;
    }
    // 光线与一些三角形求交..
    inline void hit_test_faces(Ray &ray, std::vector<Face3D *> &faces)
    {
        for (Face3D *face : faces)
        {
            glm::dvec3 d = ray.d;
            glm::dvec3 p0 = face->v1;
            glm::dvec3 e1 = face->v2 - p0;
            glm::dvec3 e2 = face->v3 - p0;
            glm::dvec3 q = glm::cross(d, e2);
            double a = glm::dot(e1, q);
            if (a == 0)
                continue;
            double f = 1 / a;
            glm::dvec3 s = ray.p - p0;
            double u = f * glm::dot(s, q);
            if (u < 0)
                continue;
            glm::dvec3 r = glm::cross(s, e1);
            double v = f * glm::dot(d, r);
            if (v < 0 || 1 - u - v < 0)
                continue;
            double t = f * glm::dot(e2, r);
            if (t > 1e-10 && t <= ray.hit_time)
            {
                ray.hit_time = t;
                ray.hit_face = face; // 如果相交，指针不为 NULL
                ray.barycentric = glm::dvec3(1.0 - u - v, u, v);
            }
        }
    }
    inline void insert_face(Face3D *f)
    {
        faces.push_back(f);
        A.x = min_(f->A.x, A.x);
        A.y = min_(f->A.y, A.y);
        A.z = min_(f->A.z, A.z);
        B.x = max_(f->B.x, B.x);
        B.y = max_(f->B.y, B.y);
        B.z = max_(f->B.z, B.z);
    }
    inline void build_subtree()
    {
        if (faces.size() < 9)
            return;
        left = new BVH();
        right = new BVH();
        glm::dvec3 size = abs(B - A);
        if (size.x > size.y && size.x > size.z)
            std::sort(std::begin(faces), std::end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.x < b->A.x; });
        else if (size.y > size.z)
            std::sort(std::begin(faces), std::end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.y < b->A.y; });
        else
            std::sort(std::begin(faces), std::end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.z < b->A.z; });
        unsigned i = 0;
        for (; i < faces.size() / 2; ++i)
            left->insert_face(faces[i]);
        for (; i < faces.size(); ++i)
            right->insert_face(faces[i]);

        faces.clear();
        left->build_subtree();
        right->build_subtree();
    }
    ~BVH()
    {
        if (left)
        {
            delete left;
            left = NULL;
        }
        if (right)
        {
            delete right;
            right = NULL;
        }
    }
};

class Render
{
public:
    // coordinate transformation, coord_world = rotate * coord_camera + position
    double camera_distance = 100.0;                   // 像距，越大则物体越大..
    glm::dvec3 camera_position = glm::dvec3(0, 0, 0); //
    glm::dmat3x3 camera_rotate = glm::dmat3x3(1.0);   //
    FrameBuffer fb;                                   // frame buffer
    glm::dvec3 background = glm::dvec3(0, 0, 0);      //
    double adapted_lum = 1;                           // hdr to ldr

    std::vector<Obj *> objs;
    std::vector<Face3D> faces;
    BVH bvh;
    BRDF brdf;
    Picture skybox;
    glm::dmat3x3 skybox_rotate = glm::dmat3x3(1.0);

    clock_t timer;
    int n_threads = 1;
    bool _init = true;
    int n_sample = 0;

    Render(int w, int h) : fb(FrameBuffer(w, h))
    {
        n_threads = omp_get_max_threads() - 2;
        omp_set_num_threads(n_threads);
    }
    ~Render() {}

    inline void set_camera(glm::dvec3 position, glm::dvec3 lookat, glm::dvec3 up, double view)
    {
        camera_position = position;
        camera_distance = fb.h_ / 2.0 / tan(view / 180.0 * std::_Pi / 2.0);
        glm::dvec3 z = lookat - position;
        glm::dvec3 y = up;
        glm::dvec3 x = glm::cross(y, z);
        z = glm::cross(x, y);
        z = glm::normalize(z);
        y = glm::normalize(y);
        x = glm::normalize(x);
        camera_rotate[0][0] = x[0];
        camera_rotate[1][0] = x[1];
        camera_rotate[2][0] = x[2];
        camera_rotate[0][1] = y[0];
        camera_rotate[1][1] = y[1];
        camera_rotate[2][1] = y[2];
        camera_rotate[0][2] = z[0];
        camera_rotate[1][2] = z[1];
        camera_rotate[2][2] = z[2];
    }
    inline void add_obj(Obj &Model) { objs.push_back(&Model); }
    void set_skybox(std::string filename) { skybox.load(filename); }

    inline glm::dvec3 skybox_color(glm::dvec3 &d)
    {
        if (skybox.empty_)
            return background;

        glm::vec3 dir_ = d * skybox_rotate;
        double y = 1.0 - glm::acos(dir_.y) / std::_Pi;
        double x = 0;
        if (dir_.z != 0)
            x = dir_.z > 0 ? glm::atan(dir_.x / dir_.z) : std::_Pi + glm::atan(dir_.x / dir_.z);

        x = 0.75 - x / (std::_Pi * 2);
        // std::cout << glm::atan(-1) << std::endl;

        return skybox.Sample2D(glm::dvec2(x, y));
    }
    int render()
    {
        if (_init) // each frame
        {
            timer = clock();
            n_sample = 0;
            _init = false;
            bvh = BVH();
            fb.clear();
            faces.clear();

            for (Obj *i : objs)
                transform_face(i);
            for (Face3D &i : faces)
                bvh.insert_face(&i);
            bvh.build_subtree();
            std::cout << "BVH: " << get_time_ms()
                      << " ms\tfaces:" << faces.size() << '\n';
        }

#pragma omp parallel for
        for (int x = 0; x < fb.w_; ++x) // current line
        {
            for (int y = 0; y < fb.h_; ++y) // current pixel
            {
                glm::dvec3 color(0, 0, 0);
                for (int i = 0; i < NS; ++i)
                    color += ray_casting(x, y);
                fb.accumulate(x, y, color);
            }
        }
        n_sample += NS;
        fb.show(n_sample, adapted_lum);
        std::cout << "[simple] " << n_sample << "\ttime = " << get_time_ms() << " ms\n";
        return n_sample;
    }
    inline glm::dvec3 ray_casting(int x, int y)
    {
        Ray ray;
        ray.p = camera_position;
        double xx = fb.w_ / 2 - x + brdf.randf2() * 0.2;
        double yy = y - fb.h_ / 2 + brdf.randf2() * 0.2;
        double zz = camera_distance;
        glm::dvec3 direction(xx, yy, zz);
        ray.d = glm::normalize(direction * camera_rotate);
        ray_tracing(ray);
        return ray.color;
    }
    inline void ray_tracing(Ray &ray)
    {
        if (ray.deepth > RAY_DEEPTH)
        {
            ray.color = glm::dvec3(0, 0, 0);
            return;
        }
        ray.hit_face = NULL;
        ray.hit_time = DBL_MAX;
        bvh.hit_test(ray);
        // 40% total time
        // bvh.hit_test(ray);
        // bvh.hit_test(ray);
        if (!ray.hit_face)
        {
            ray.color = skybox_color(ray.d);
            return;
        }
        ray.deepth += 1;
        Face3D *f = ray.hit_face;

        Material &m = f->m;
        glm::dvec3 coord = ray.barycentric;
        glm::dvec3 N_ = glm::normalize(coord.x * f->n1 + coord.y * f->n2 + coord.z * f->n3);
        glm::dvec2 uv_ = coord.x * f->uv1 + coord.y * f->uv2 + coord.z * f->uv3;
        uv_.x -= std::floor(uv_.x);
        uv_.y -= std::floor(uv_.y);
        glm::dvec3 V_ = ray.d;                // visual dir
        glm::dvec3 R_ = glm::reflect(V_, N_); // reflect dir

        // 30% time
        glm::dvec3 K;  // 1 - 吸收率..
        glm::dvec3 L_; // simple direction
        double rand_ = brdf.randf1();
        double pkd = m.kd;
        double pks = pkd + m.ks;
        double pkr = pks + m.kr;
        if (rand_ < pkd)
        {
            L_ = brdf.diffuse(N_);
            if (m.Map_Kd.empty_)
                K = m.Kd;
            else
                K = m.Map_Kd.Sample2D(uv_);
        }
        else if (rand_ < pks)
        {
            L_ = brdf.specular(R_, m.Ns);
            K = m.Ks;
        }
        else if (rand_ < pkr)
        {
            K = m.Kr;
            glm::dvec3 N__ = N_;
            double Nr = 1.0 / m.Nr;
            if (glm::dot(V_, N_) > 0)
            {
                N__ = -N_;
                Nr = 1.0 / Nr;
            }
            L_ = brdf.refract(V_, N__, Nr);
        }
        else
        {
            ray.color = m.Le;
            return;
        }

        ray.p = ray.p + ray.hit_time * V_;
        ray.d = L_;
        ray_tracing(ray);

        ray.color = m.Le + ray.color * K;
    }
    void transform_face(Obj *o)
    {
        glm::dmat3x3 rotate = o->pose * o->scale; // 物体缩放..
        glm::dvec3 move;                          // 物体平移..
        move[0] = o->pose[0][3];
        move[1] = o->pose[1][3];
        move[2] = o->pose[2][3];
        std::vector<glm::dvec3> vertex_new;
        std::vector<glm::dvec3> norm_new;
        // 顶点坐标转换.
        for (int i = 0; i < o->model.vertex_.size(); ++i)
            vertex_new.push_back(o->model.vertex_[i] * rotate + move);
        // 顶点法向量转换.
        for (int i = 0; i < o->model.norm_.size(); ++i)
            norm_new.push_back(o->model.norm_[i] * rotate);
        for (glm::imat3x4 i : o->model.face_)
        {
            glm::dvec3 v1 = vertex_new[i[0][0]];
            glm::dvec3 v2 = vertex_new[i[1][0]];
            glm::dvec3 v3 = vertex_new[i[2][0]];
            glm::dvec3 n1 = norm_new[i[0][1]];
            glm::dvec3 n2 = norm_new[i[1][1]];
            glm::dvec3 n3 = norm_new[i[2][1]];
            glm::dvec2 uv1 = o->model.uv_coord_[i[0][2]];
            glm::dvec2 uv2 = o->model.uv_coord_[i[1][2]];
            glm::dvec2 uv3 = o->model.uv_coord_[i[2][2]];
            // bounding box
            glm::dvec3 A = {min3(v1.x, v2.x, v3.x), min3(v1.y, v2.y, v3.y), min3(v1.z, v2.z, v3.z)};
            glm::dvec3 B = {max3(v1.x, v2.x, v3.x), max3(v1.y, v2.y, v3.y), max3(v1.z, v2.z, v3.z)};

            faces.push_back({v1, v2, v3, n1, n2, n3, uv1, uv2, uv3, A, B,
                             o->model.material_[i[0][3]]});
        }
    }
    uint32_t *get_framebuffer() { return fb.fb_; }
    double get_time_ms()
    {
        double ret = (double)(clock() - timer) * 1000.0 / CLOCKS_PER_SEC;
        timer = clock();
        return ret;
    }
};

#endif
