# 数据结构 


#### mtl 材质文件

```py
# 定义一个名为 'xxx'的材质
newmtl xxx
# 材质的环境光（ambient color）
Ka 0 0 0
# 散射光（diffuse color）用Kd
Kd 0.784314 0.784314 0.784314
# 镜面光（specular color）用Ks
Ks 0 0 0
# 折射值 可在0.001到10之间进行取值。若取值为1.0，光在通过物体的时候不发生弯曲。玻璃的折射率为1.5。
Ni 1
# 反射指数 定义了反射高光度。该值越高则高光越密集，一般取值范围在0~1000。
Ns 400
# 滤光透射率
Tf 1 1 1
# 渐隐指数描述 参数factor表示物体融入背景的数量，取值范围为0.0~1.0，取值为1.0表示完全不透明，取值为0.0时表示完全透明。
d 1
# 为漫反射指定颜色纹理文件
map_Kd test_vt.bmp
```
newmtl：代表材质，以下皆为该材质的属性参数
Ns：高光反射系数，值越高则高光越密集
NI：指定材质表面的光密度，即折射值
d：表示物体融入背景的数量，取值范围为0.0~1.0，取值为1.0表示完全不透明，取值为0.0时表示完全透明
Tr：定义材质的alpha透明度
Tf：材质的透射滤波（transmission filter），对应数据为r，g，b值
illum： 照明度（illumination），后面可接0~10范围内的数字参数
Ka: 环境光（ambient color）
Kd: 散射光（diffuse color）
Ks: 镜面光（specular color）
Ke：放射光（emissive color）
map_Ka：环境光所采样的纹理贴图路径，在.obj模型文件的根目录下
map_Kd：漫反射光所采样的纹理贴图路径

### obj 文件 

```py
# obj对应的材质文件
# mtllib testvt.mtl
# 组名称
g default
# o 对象名称(Object name)
o testvt.obj
# 顶点
v -0.5 -0.5 0.1
v -0.5 -0.5 -0.1
v 0 0.5 0.1
v 0 0.5 -0.1
v 0.5 -0.5 0.1
v 0.5 -0.5 -0.1
# 纹理坐标
vt 0 1
vt 1 1
vt 0.5 0
# 顶点法线
vn 0 0 1
vn 0 0 -1
# 当前图元所用材质
usemtl Default
# s Smooth shading across polygons is enabled by smoothing groups.
# Smooth shading can be disabled as well.
s off
# v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3(索引起始于1)    
f 1/1/1 5/2/1 3/3/1
f 6/2/2 2/1/2 4/3/2
```

### 求交 

$$
p(t) = p_1 + tu \\ n \cdot (p-p_2) = 0 \\ n \cdot (p_1 - p_2 + tu) = 0 \\
t (n \cdot u ) = n \cdot (p_2 - p_1) 
$$


### 重心坐标系求交 

光线: $p+td$
三角: $p_0, p_1, p_2$

$$
p + td = (1-u-v) p_0  + u p_1 + vp_2
$$

$e_1 = p_1 - p_0$
$e_2 = p_2 - p_0$
$q = d \times e_2$
$a = e_1 \cdot q$ ,  if $a = 0$ 平行
$f = 1/a$
$s = p - p_0$
$u = f (s\cdot q)$
$r = s \times e_1$
$v = f (d\cdot r)$
$t = f(e_2 \cdot r)$


### 漫反射

半球面均匀分布随机向量:

- 生成边长为 2 的立方体内随机点
- 得到半径为 1 的球体内的随机点 
- 如果与入射光线不在同一侧，取反



### AABBCC 求交 

光线依次穿过三个区间，

$p.x + t.x * d.x = A.x$


### 法向量坐标系转换

法向量坐标系: 法向量为 $z$ 轴的球坐标系

该球坐标系中某点 $(\phi, \theta, 1)$ 的笛卡尔坐标 $(\sin\theta\cos\phi, \sin\theta\sin\phi,\cos\theta)$

原坐标系变换成法向量坐标系的其中一种转换：绕 $z$ 轴旋转 $\phi_0$, 再绕新 $y$ 轴旋转 $\theta_0$ ，旋转矩阵 $$\begin{bmatrix}
  \cos\theta_0\cos\phi_0 & -\sin\phi_0 & \sin\theta_0\cos\phi_0 \\
  \cos\theta_0\sin\phi_0 & \cos\phi_0 & \sin\theta_0\sin\phi_0  \\
  -\sin\theta_0  & 0 &  \cos\theta_0  
\end{bmatrix}$$


设法向量 $(x, y, z)$, $\theta_0 = \arccos z$, $\phi_0 = \arctan \frac{y}{x}$

由 
$$ \left \{
\begin{array}{l}
  x = \sin\theta_0\cos\phi_0 \\
  y = \sin\theta_0\sin\phi_0 \\
  z = \cos\theta_0 
\end{array} \right .
$$
得 
$$
\left\{
  \begin{array}{l}
    \sin\theta_0 = \pm\sqrt{1-z^2} \\
    \cos\phi_0 = x/ \sin\theta_0 \\
    \sin\phi_0 = y / \sin\theta_0
  \end{array}
\right .
$$



### 折射

### hdr 转 ldr 

tone mapping 

### skybox 

球面坐标

