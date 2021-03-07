#ifndef _RAY_TRACING_H_
#define _RAY_TRACING_H_

#include <glm.hpp>
#include <ctime>
#include <random>
#include "omp.h"
#include "model.h"

#define NS 10

// 3D obj + pose and scale
struct Obj
{
    Model &obj_;
    mat4x4 pose = mat4x4(1.0f);
    float scale = 1.0f;

    Obj(Model &m_) : obj_(m_) {}

    void set_pose(mat4x4 pose_, float s_)
    {
        pose = pose_;
        scale = s_;
    }
};

class FrameBuffer
{
public:
    uint32_t *fb_; // RGBA from memory low to high , y: ^ x: ->
    vec3 *fb2_;
    int w_, h_;
    FrameBuffer(int w, int h) : w_(w), h_(h)
    {
        fb_ = new uint32_t[w * h];
        fb2_ = new vec3[w * h];
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
    inline void fill(uint32_t color)
    {
        for (int i = 0; i < w_ * h_; ++i)
            fb_[i] = color;
        for (int i = 0; i < w_ * h_; ++i)
            fb2_[i] = vec3(0, 0, 0);
    }
    inline void set_pixel(int x, int y, uint32_t color)
    {
        fb_[y * w_ + x] = color;
    }
    inline void set_pixel(int x, int y, vec3 color)
    {
        fb2_[y * w_ + x] += color;
    }
    inline void to_fb(float n)
    {
        for (int i = 0; i < w_ * h_; ++i)
            fb_[i] = toRGB(fb2_[i] / n);
    }
};

struct RandomVec
{
    random_device rd;
    mt19937 gen;
    uniform_real_distribution<float> dis;
    uniform_real_distribution<float> dis1;
    RandomVec() : gen(rd()), dis(-1, 1), dis1(0, 1) {}
    RandomVec(float a, float b) : gen(rd()), dis(a, b), dis1(0, 1) {}
    inline float randf()
    {
        return dis(gen);
    }
    inline float randf1()
    {
        return dis1(gen);
    }
    inline vec3 rand_sphere()
    {
        vec3 ret;
        do
        {
            ret.x = randf();
            ret.y = randf();
            ret.z = randf();
        } while (dot(ret, ret) > 1.0f);
        return normalize(ret);
    }
    inline vec3 rand_diffuse(vec3 in, vec3 n)
    {
        vec3 out = rand_sphere();
        if (dot(n, in) * dot(n, out) > 0)
            return -out;
        else
            return out;
    }
    inline mat3x3 norm_to_world(vec3 n)
    {
        // cout << print_vec3(n) << endl;
        if (n.z > 0.9999f)
            return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        if (n.z < -0.9999f)
            return {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}};
        float cos_theta0 = n.z;
        float theta0 = acos(cos_theta0);
        float sin_theta0 = sin(theta0);

        float cos_phi0 = n.x / sin_theta0;
        float sin_phi0 = n.y / sin_theta0;

        return {{cos_theta0 * cos_phi0, -sin_phi0, n.x},
                {cos_theta0 * sin_phi0, cos_phi0, n.y},
                {-sin_theta0, 0.0f, n.z}};
    }
    inline vec3 diffuse(vec3 n)
    {
        float u1 = randf1();
        float sin_theta = sqrtf(1.0f - u1);
        float cos_theta = sqrtf(u1);

        float phi = 2 * _Pi * randf1();
        float sin_phi = sin(phi);
        float cos_phi = cos(phi);
        vec3 out = vec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
        out = normalize(out * norm_to_world(n));
        if (dot(n, out) < 0)
            out = -out;
        return out;
    }
    inline vec3 specular(vec3 ref, vec3 n, int Ns)
    {
        float u1 = pow(randf1(), 1.0f / (Ns + 1.0f));
        float cos_theta = u1;
        float sin_theta = sqrtf(1.0f - u1 * u1);

        float phi = 2 * _Pi * randf1();
        float sin_phi = sin(phi);
        float cos_phi = cos(phi);

        vec3 out = vec3(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
        out = normalize(out * norm_to_world(ref));
        if (dot(n, out) < 0)
            out = -out;
        return out;
    }
    inline vec3 refract(vec3 &in, vec3 &n, float eta)
    {
        float N_dot_I = dot(n, in);
        float k = 1.f - eta * eta * (1.f - N_dot_I * N_dot_I);
        if (k < 0.f)
            return reflect(in, n);
        else
            return eta * in - (eta * N_dot_I + sqrtf(k)) * n;
    }
};

class BVH;
// vertex + material
struct Face3D
{
    vec3 v1, v2, v3;    // 3 vertexes
    vec3 n1, n2, n3;    // 3 norms
    vec2 uv1, uv2, uv3; // 3 uv_coordinates
    vec3 A, B;          // bounding box
    Material &m;

    inline void print_vertex()
    {
        cout << print_vec3(v1) << print_vec3(v2) << print_vec3(v3) << endl;
    }
};
struct HittedBox
{
    BVH *box;
    float t1, t2;
};
struct Ray
{
    vec3 p;         // start point
    vec3 d;         // direction
    vec3 color;     // rgb
    int deepth = 0; // 反射次数..
    Face3D *hit_face = NULL;
    float hit_time = FLT_MAX;
    vec3 barycentric;
    vector<HittedBox> hit_results;
};

class BVH
{
public:
    BVH *left = NULL;
    BVH *right = NULL;                        // child node
    vec3 A = vec3(FLT_MAX, FLT_MAX, FLT_MAX); // bounding box
    vec3 B = vec3(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    vector<Face3D *> faces;

    BVH() {}

    inline void hit_test(Ray &ray)
    {
        vec3 temp = ray.d;
        if (ray.d.x == 0)
            ray.d.x = 1e-20;
        if (ray.d.y == 0)
            ray.d.y = 1e-20;
        if (ray.d.z == 0)
            ray.d.z = 1e-20;
        hit_test_box_(ray);
        ray.d = temp;
        vector<HittedBox> &result = ray.hit_results;
        if (result.empty())
            return;
        sort(begin(result), end(result), [](HittedBox &a, HittedBox &b) -> bool { return a.t1 < b.t1; });
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
        vec3 t1 = (A - ray.p) / ray.d;
        vec3 t2 = (B - ray.p) / ray.d;
        if (t1.x > t2.x)
            swap(t1.x, t2.x);
        if (t1.y > t2.y)
            swap(t1.y, t2.y);
        if (t1.z > t2.z)
            swap(t1.z, t2.z);
        float box_t1 = max3(t1.x, t1.y, t1.z);
        float box_t2 = min3(t2.x, t2.y, t2.z);
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
    inline void hit_test_faces(Ray &ray, vector<Face3D *> &faces)
    {
        for (Face3D *face : faces) // 是否与三角形相交..
        {
            vec3 d = ray.d;
            vec3 p0 = face->v1;
            vec3 e1 = face->v2 - p0;
            vec3 e2 = face->v3 - p0;
            vec3 q = cross(d, e2);
            float a = dot(e1, q);
            if (a == 0)
                continue;
            float f = 1 / a;
            vec3 s = ray.p - p0;
            float u = f * dot(s, q);
            if (u < 0)
                continue;
            vec3 r = cross(s, e1);
            float v = f * dot(d, r);
            if (v < 0 || 1 - u - v < 0)
                continue;
            float t = f * dot(e2, r);
            if (t > 1e-5 && t <= ray.hit_time)
            {
                ray.hit_time = t;
                ray.hit_face = face;
                ray.barycentric = vec3(1.0f - u - v, u, v);
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
        if (faces.size() < 12)
            return;
        left = new BVH();
        right = new BVH();
        vec3 size = abs(B - A);
        if (size.x > size.y && size.x > size.z)
            sort(begin(faces), end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.x < b->A.x; });
        else if (size.y > size.z)
            sort(begin(faces), end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.y < b->A.y; });
        else
            sort(begin(faces), end(faces), [](Face3D *a, Face3D *b) -> bool { return a->A.z < b->A.z; });
        // cout << "current_faces:" << faces.size() << "\t box: " << print_vec3(A) << print_vec3(B) << endl;
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
    mat4x4 camera = mat4x4(1.0f);   // camera pose
    float camera_distance = 100.0f; // 像距，越大则物体越大..
    // coordinate transformation, coord_world = rotate * coord_camera + position
    vec3 camera_position = vec3(0, 0, 0);
    mat3x3 camera_rotate = mat3x3(1.0f);
    FrameBuffer fb;     // frame buffer
    vec3 background;    // rgb, 0-1
    int SAMPLES = 1000; // rays per pixel

    vector<Obj *> objs;
    vector<Face3D> faces;
    BVH faces_tree;

    clock_t timer;
    int n_threads;
    RandomVec rand0;
    bool init = true;
    int n_sample = 0;

    Render(int w, int h) : fb(FrameBuffer(w, h))
    {
        n_threads = omp_get_max_threads() - 2;
        // for (int i = 0; i < n_threads; ++i)
        // {
        //     rands.push_back(RandomVec());
        // }
        omp_set_num_threads(n_threads);
    }
    ~Render() {}
    void set_camera(mat4x4 c, float scale)
    {
        camera = c;
        camera_distance = scale;
    }
    void set_skybox(string filename)
    {
        return;
    }
    vec3 skybox_color(vec3 dir)
    {
        return background;
    }
    void set_camera(vec3 position, vec3 lookat, vec3 up, float view)
    {
        camera_position = position;
        camera_distance = fb.h_ / 2.0f / tan(view / 180.0f * _Pi / 2.0f);
        vec3 z = lookat - position;
        vec3 y = up;
        vec3 x = cross(y, z);
        z = cross(x, y);
        z = normalize(z);
        y = normalize(y);
        x = normalize(x);
        camera_rotate[0][0] = x[0];
        camera_rotate[1][0] = x[1];
        camera_rotate[2][0] = x[2];
        camera_rotate[0][1] = y[0];
        camera_rotate[1][1] = y[1];
        camera_rotate[2][1] = y[2];
        camera_rotate[0][2] = z[0];
        camera_rotate[1][2] = z[1];
        camera_rotate[2][2] = z[2];
        // cout << print_mat3x3(camera_rotate) << view << camera_distance << endl;
    }
    void set_background(vec3 in)
    {
        background = in;
    }
    void add_obj(Obj &Model)
    {
        objs.push_back(&Model);
    }
    double get_time_ms()
    {
        double ret = (double)(clock() - timer) * 1000.0 / CLOCKS_PER_SEC;
        timer = clock();
        return ret;
    }

    int render()
    {
        // cout << print_vec3(refract( vec3(0,0.8,-0.6), vec3(0,0,1), 1.5f));
        // mat3x3 aa = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        // for (int i = 0; i < 100; ++i)
        // {
        //     cout <<print_mat3x3(rand0.norm_to_world(vec3(0,0,1))) << endl;
        //     cout << print_vec3(vec3(1,2,3) * aa) << endl;
        // }
        // while (1)
        // {
        // }
        // fb.fill(toRGB(background));
        // for (int i = 0; i < 100; ++i)
        // {
        //     vec3 xx = rand1.rand_diffuse(vec3(1,1,-1),vec3(0,0,-1));
        //     // vec3 xx = rand1.rand_sphere();
        //     cout << print_vec3(xx) << endl;
        //     xx = (xx + 2.0f) * 50.0f;
        //     fb.set_pixel(xx.x, xx.y, 0xff00ff00);
        // }
        // return;
        if (init)
        {
            timer = clock();
            n_sample = 0;
            init = false;
            faces_tree = BVH();
            fb.fill(toRGB(background));
            faces.clear();

            for (Obj *i : objs)
                transform_face(i);
            for (Face3D &i : faces)
                faces_tree.insert_face(&i);
            faces_tree.build_subtree();
            cout << "time BVH = " << get_time_ms() << " ms\tn_faces:" << faces.size() << endl;
        }
        if (faces.empty() || n_sample > SAMPLES)
            return -1;

        timer = clock();
#pragma omp parallel for
        for (int x = 0; x < fb.w_; ++x) // current line
        {
            for (int y = 0; y < fb.h_; ++y) // current pixel
            {
                vec3 color = vec3(0, 0, 0);
                for (int i = 0; i < NS; ++i) // number of samples
                    color += ray_casting(x, y);
                fb.set_pixel(x, y, color);
            }
        }
        n_sample += NS;
        fb.to_fb(n_sample);
        cout << "n_sample = " << n_sample << "\ttime = " << get_time_ms() << " ms\n";
        return n_sample;
    }
    inline vec3 ray_casting(int x, int y)
    {
        Ray ray;
        ray.p = camera_position;
        float xx = fb.w_ / 2 - x + rand0.randf() * 0.1f;
        float yy = y - fb.h_ / 2 + rand0.randf() * 0.1f;
        float zz = camera_distance;
        vec3 direction = vec3(xx / zz, yy / zz, 1.0f);
        ray.d = normalize(direction * camera_rotate);

        // cout << print_vec3(ray.p) << print_vec3(ray.d) << endl;
        ray_tracing(ray);
        return ray.color;
    }
    inline void ray_tracing(Ray &ray)
    {
        if (ray.deepth > 8)
        {
            ray.color = vec3(0, 0, 0);
            return;
        }
        float Na = ray.deepth == 0;
        ray.hit_face = NULL;
        ray.hit_time = FLT_MAX;
        faces_tree.hit_test(ray); // 50% total time
        // faces_tree.hit_test(ray);
        // faces_tree.hit_test(ray);
        if (!ray.hit_face)
        {
            ray.color = skybox_color(ray.d);
            return;
        }
        ray.deepth += 1;
        Face3D *f = ray.hit_face;
        // ray.color = f->m.Kd;
        // return;
        Material &m = f->m;
        float t = ray.hit_time;
        vec3 coord = ray.barycentric;
        vec3 N_ = normalize(coord.x * f->n1 + coord.y * f->n2 + coord.z * f->n3); // norm at hit point
        vec2 uv_ = coord.x * f->uv1 + coord.y * f->uv2 + coord.z * f->uv3;        // texture coord
        vec3 V_ = ray.d;                                                          // visual dir

        // ray.p = ray.p + t * V_;
        // ray.d = reflect(V_, N_);
        // ray_tracing(ray);
        // ray.color = vec3(1,1,1);
        // if (!m.Map_Kd.empty_)
        // {
        //     ray.color = m.Map_Kd.Sample2D(uv_);
        // }
        // return;

        /******* 30% time ******/
        vec3 K;  // coeff
        vec3 L_; // simple direction
        float rand_ = rand0.randf1();
        float pkd = m.kd;
        float pks = pkd + m.ks;
        float pkr = pks + m.kr;
        // cout << "rand " << rand_ << "\t Kd " << pkd << "\t Ks " << pks << endl;
        if (rand_ < pkd)
        {
            L_ = rand0.diffuse(N_);
            // cout << "N_ " <<  print_vec3(N_) << " L_ " << print_vec3(L_) << endl;
            if (m.Map_Kd.empty_)
                K = m.Kd;
            else
            {
                // cout << "uv: " << uv_.x << ' ' << uv_.y << endl;
                K = m.Map_Kd.Sample2D(uv_);
            }
        }
        else if (rand_ < pks)
        {
            // L_ = reflect(V_, N_) + rand0.rand_sphere() * pow(rand0.randf(), m.Ns * 0.1f);
            // L_ = normalize(L_);
            // cout << m.ks << "   Ks:" << print_vec3(m.Ks) << endl;
            L_ = rand0.specular(reflect(V_, N_), N_, m.Ns);
            K = m.Ks;
        }
        else if (rand_ < pkr)
        {
            K = m.Kr;
            vec3 N__ = N_;
            float Nr = m.Nr;
            float cos_theta1 = -dot(V_, N_);
            if (cos_theta1 < 0)
            {
                N__ = -N_;
                cos_theta1 = -cos_theta1;
                Nr = 1.0f / Nr;
            }
            L_ = rand0.refract(V_, N__, 1.0f / Nr);
            L_ = normalize(L_);
            // if (cos_theta1 > 0.9999f)
            //     L_ = V_;
            // else
            // {
            //     float sin_theta1 = sqrtf(1 - cos_theta1 * cos_theta1);
            //     float sin_theta2 = sin_theta1 / Nr;
            //     if (sin_theta2 > 0.9999f)
            //         L_ = reflect(V_, N_);
            //     else
            //     {
            //         float cos_theta2 = sqrtf(1 - sin_theta2 * sin_theta2);
            //         // cout << "cos1:" << cos_theta1 << " cos2:" << cos_theta2 << " Nr:" << Nr << endl;
            //         vec3 t1 = (V_ + N_ * cos_theta1) / Nr;
            //         vec3 t2 = -N_ * cos_theta2;
            //         L_ = normalize(t1 + t2);
            //     }
            // }
        }
        else
        {
            ray.color = m.Le;
            return;
        }

        ray.p = ray.p + t * V_;
        ray.d = L_;
        ray_tracing(ray);

        ray.color = m.Le + ray.color * K;
        //   (f_min->m.Kd * max_(dot(N_, L_), .0f) +
        //    f_min->m.Ks * pow(max_(dot(H_, N_), .0f), f_min->m.Ns));
    }
    void transform_face(Obj *o)
    {
        // mat4x4 inverse_camera = transpose(camera);
        // mat3x3 rotate = inverse_camera;
        // vec3 move;
        // move[0] = camera[0][3];
        // move[1] = camera[1][3];
        // move[2] = camera[2][3];
        // // move = -rotate * move;
        // move = -move * rotate;
        // inverse_camera[0][3] = move.x;
        // inverse_camera[1][3] = move.y;
        // inverse_camera[2][3] = move.z;
        // inverse_camera[3][0] = 0;
        // inverse_camera[3][1] = 0;
        // inverse_camera[3][2] = 0;
        // inverse_camera = inverse_camera * o->coordinate; // 转换到相机空间.
        // inverse_camera = o->coordinate * inverse_camera;
        mat3x3 rotate = o->pose * o->scale; // 物体缩放..
        vec3 move;                          // 物体平移..
        move[0] = o->pose[0][3];
        move[1] = o->pose[1][3];
        move[2] = o->pose[2][3];
        vector<vec3> vertex_new;
        vector<vec3> norm_new;
        // 顶点坐标转换.
        for (int i = 0; i < o->obj_.vertex_.size(); ++i)
            vertex_new.push_back(o->obj_.vertex_[i] * rotate + move);
        // 顶点法向量转换.
        for (int i = 0; i < o->obj_.norm_.size(); ++i)
            norm_new.push_back(o->obj_.norm_[i] * rotate);
        for (imat3x4 i : o->obj_.face_)
        {
            vec3 v1 = vertex_new[i[0][0]];
            vec3 v2 = vertex_new[i[1][0]];
            vec3 v3 = vertex_new[i[2][0]];
            vec3 n1 = norm_new[i[0][1]];
            vec3 n2 = norm_new[i[1][1]];
            vec3 n3 = norm_new[i[2][1]];
            vec2 uv1 = o->obj_.uv_coord_[i[0][2]];
            vec2 uv2 = o->obj_.uv_coord_[i[1][2]];
            vec2 uv3 = o->obj_.uv_coord_[i[2][2]];
            // bounding box
            vec3 A = {min3(v1.x, v2.x, v3.x), min3(v1.y, v2.y, v3.y), min3(v1.z, v2.z, v3.z)};
            vec3 B = {max3(v1.x, v2.x, v3.x), max3(v1.y, v2.y, v3.y), max3(v1.z, v2.z, v3.z)};

            faces.push_back({v1, v2, v3, n1, n2, n3, uv1, uv2, uv3, A, B,
                             o->obj_.material_[i[0][3]]});
        }
    }
    uint32_t *get_framebuffer()
    {
        return fb.fb_;
    }
};

#endif
