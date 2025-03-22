#pragma once
#include <vector>
#include <glm/glm.hpp>

#define FromDegToRad 0.0174532925f
#define BG_COLOR 1.0f, 1.0f, 0.0f
#define OUTPUT_FILE_NAME "result.ppm"

struct Material {
	glm::vec3 ka, kd, ks;
	float n;
	float kr, kt, index; // need more?
	Material() {
		ka = kd = ks = glm::vec3(0.0f, 0.0f, 0.0f);
		n = kr = kt = index = 0.0f;
	}
	Material(glm::vec3 ka, glm::vec3 kd, glm::vec3 ks, float n, float kr, float kt, float index) : 
		ka{ ka }, kd{ kd }, ks{ ks }, n{ n }, kr{ kr }, kt{ kt }, index{ index } {}
};

struct Light {
	glm::vec4 position; // position.w = 0: directional light, position.w != 0: point light
	glm::vec3 la, ld, ls;
	// and more
};

struct Rectanglee {
	glm::vec3 a, u, v;
	Material material;
	Rectanglee(glm::vec3 a, glm::vec3 u, glm::vec3 v, Material material) : a{ a }, u{ u }, v{ v }, material{ material } {}
};

struct Sphere {
	glm::vec3 center;
	float radius;
	Material material;
	Sphere(glm::vec3 center, float radius, Material material) : center{ center }, radius{ radius }, material{ material } {}
};

struct Box {
	glm::vec3 min, max;
	Material material;
	Box(glm::vec3 min, glm::vec3 max, Material material) : min{ min }, max{ max }, material{ material } {}
};

struct TriangularMesh {
	int n_triangles;
	float* triangle_list; // n_triangles * 3 * 6 * sizeof(float)
	// x y z nx ny nz x y z nx ny nz x y z nx ny nz ...
	Material material;
	TriangularMesh(int n_triangles, float* triangle_list, Material material) : n_triangles{ n_triangles }, 
		triangle_list{ triangle_list }, material{ material } {}
};

struct Scene {
	glm::vec3 global_ambient; 
	std::vector<Light> lights;
	std::vector<Rectanglee> rectangles;
	std::vector<Sphere> spheres;
	std::vector<Box> boxes;   
	std::vector<TriangularMesh> triangular_meshes;
	Scene() : global_ambient{ glm::vec3(0.2f, 0.2f, 0.2f) }, lights{ },
		spheres{}, boxes{}, triangular_meshes{ } {}
	void add_sphere(glm::vec3 center, float radius, Material material);
	void add_rectangles(glm::vec3 a, glm::vec3 u, glm::vec3 v, Material material);
	void add_boxes(glm::vec3 min, glm::vec3 max, Material material);
	void add_triangular_meshes(int n_triangles, float* triangle_list, Material material);
	void initialize();
};

struct Camera {
	glm::vec3 e, u, v, n;
	float fovy, aspect;
	Camera() : e{ glm::vec3{0.0f, 0.0f, 0.0f} }, u{ glm::vec3{1.0f, 0.0f, 0.0f} }, v{ glm::vec3{0.0f, 1.0f, 0.0f} },
		n{ glm::vec3{0.0f, 0.0f, 1.0f} }, fovy{ 45.0f }, aspect{ 1.0f } {}
	Camera(glm::vec3 e, glm::vec3 u, glm::vec3 v, glm::vec3 n, float fovy, float aspect) : e{ e }, u{ u }, v{ v }, n{ n },
		fovy{ fovy }, aspect{ aspect } {}
};

struct Window {
	int width, height;
	float d, wh, ww;
	glm::vec3 wc, a, b, c;
	unsigned char* pixmap;
	Window() : width{ 0 }, height{ 0 } {
		Camera camera;
		prepare(camera);
	}
	Window(int width, int height, Camera camera) : width{ width }, height{ height } {
		prepare(camera);
	}
	void prepare(Camera camera);
};

struct Ray {
	glm::vec3 orig, dir;
	float tmin, tmax;
};

class Render {
	Scene* scene;
	Camera camera;
	Window window;
	glm::vec3 compute_ray_dir(int i, int j);
	void trace_ray(Ray& ray, glm::vec3* shaded_color);
	Material* find_first_hit(Ray& ray, float* t, glm::vec3* position, glm::vec3* normal);
	float intersect_sphere(Ray& ray, Sphere& sphere);
	std::pair<float, glm::vec3> intersect_rectangle(Ray& ray, Rectanglee& box);
	float intersect_box(Ray& ray, Box& box);
	float intersect_triangular_mesh(Ray& ray, TriangularMesh& box);
	void shade_first_hit(glm::vec3& pos, glm::vec3& norm, Material* material, glm::vec3* shaded_color);
	void dump_pixmap();
public:
	Render() : scene{ nullptr }, camera{ }, window{ } {}
	Render(Scene* scene, Camera camera, Window window) : scene{ scene }, camera{ camera }, window{ window} {}
	void render();
};