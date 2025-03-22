#define _CRT_SECURE_NO_WARNINGS

#include "RayTracing.h"

void Scene::add_sphere(glm::vec3 center, float radius, Material material) {
	spheres.push_back(Sphere(center, radius, material));
}

void Scene::add_rectangles(glm::vec3 a, glm::vec3 u, glm::vec3 v, Material material) {
	rectangles.push_back(Rectanglee(a, u, v, material));
}

void Scene::add_boxes(glm::vec3 min, glm::vec3 max, Material material) {
	boxes.push_back(Box(min, max, material));
}

void Scene::add_triangular_meshes(int n_triangles, float* triangle_list, Material material) {
	triangular_meshes.push_back(TriangularMesh(n_triangles, triangle_list, material));
}

void Scene::initialize() {
	Material gold = { glm::vec3{0.24725, 0.1995, 0.0745}, glm::vec3{0.75164, 0.60648, 0.22648}, 
		glm::vec3{0.628281, 0.555802, 0.366065}, 51.2, 0.0, 0.0, 1.0 };
	Material ruby = { glm::vec3{0.1745, 0.01175, 0.01175}, glm::vec3{0.61424, 0.04136, 0.04136},
	glm::vec3{0.727811, 0.626959, 0.626959}, 76.8, 0.0, 0.0, 1.0 };
	Material emerald = { glm::vec3{0.0215, 0.1745, 0.0215}, glm::vec3{0.07568, 0.61424, 0.07568},
		glm::vec3{0.633, 0.727811, 0.633}, 76.8, 0.0, 0.0, 1.0 };
	Material yellow_rubber = { glm::vec3{0.05, 0.05, 0.0}, glm::vec3{0.5, 0.5, 0.4},
	glm::vec3{0.7, 0.7, 0.04}, 10.0, 0.0, 0.0, 1.0 };
	Material silver = { glm::vec3{0.19225, 0.19225, 0.19225}, glm::vec3{0.50754, 0.50754, 0.50754},
		glm::vec3{0.508273, 0.508273, 0.508273}, 51.2, 0.0, 0.0, 1.0 };

	//add_sphere(glm::vec3(0.0f, 0.0f, 0.0f), 1.0f, gold);
	//add_sphere(glm::vec3(1.0f, 0.5f, 0.0f), 0.8f, ruby);
	add_rectangles(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f), glm::vec3(1.0f, 0.0f, -1.0f), emerald);
	add_boxes(glm::vec3(-2.0f, -2.0f, -0.5f), glm::vec3(0.0f, 0.0f, 0.0f), emerald);
	//add_triangular_meshes(n_triangles, tirangle_list, gold);
}

void Window::prepare(Camera camera) {
	pixmap = new unsigned char[width * height * 3];
	d = 1.0f;
	wc = camera.e - d * camera.n;
	wh = 2.0f * d * tan(camera.fovy * FromDegToRad / 2.0f);
	ww = camera.aspect * wh;

	a = wc - ww * camera.u / 2.0f + wh * camera.v / 2.0f;
	b = camera.u;
	c = -camera.v;
}

glm::vec3 Render::compute_ray_dir(int i, int j) {
	glm::vec3 tmp_vec3, dir;
	float length;

	for (int k = 0; k < 3; k++)
		tmp_vec3[k] = window.a[k] + ((float)i + 0.5f) * window.ww * window.b[k] / window.width
		+ ((float)j + 0.5f) * window.wh * window.c[k] / window.height;
	for (int k = 0; k < 3; k++)
		tmp_vec3[k] = window.a[k] 
		+ (i + 0.5) * window.ww * window.b[k] / window.width 
		+ (j + 0.5) * window.wh * window.c[k] / window.height - camera.e[k];
	length = sqrtf(glm::dot(tmp_vec3, tmp_vec3));  

	for (int k = 0; k < 3; k++)
		dir[k] = tmp_vec3[k] / length;
	return dir;
}
#define PIXMAP(i, j, k) *(window.pixmap + (j)*window.width*3 + (i)*3 + (k))
void Render::render() {
	Ray ray;
	ray.orig = camera.e;

	for (int j = 0; j < window.height; j++) {
		for (int i = 0; i < window.width; i++) {
			unsigned char* next_pixel = window.pixmap + 3 * j * window.width +  3 * i;
			glm::vec3 shaded_color;
			ray.dir = compute_ray_dir(i, j);
			trace_ray(ray, &shaded_color);

			*next_pixel = 255 * shaded_color[0] + 0.5;
			*(next_pixel + 1) = 255 * shaded_color[1] + 0.5;
			*(next_pixel + 2) = 255 * shaded_color[2] + 0.5;
		}
		// fprintf(stdout, "\nRAY_CASTER: DONE with scanline no. %d ...\n", j);
	}
	dump_pixmap();
}

void Render::trace_ray(Ray& ray, glm::vec3* shaded_color) {
	float t;
	glm::vec3 first_hit_pos, first_hit_norm; 
	Material* first_hit_material;
	first_hit_material = find_first_hit(ray, &t, &first_hit_pos, &first_hit_norm);
	if (first_hit_material == nullptr) {
		*shaded_color = glm::vec3(BG_COLOR); // background color
		return;
	}
	shade_first_hit(first_hit_pos, first_hit_norm, first_hit_material, shaded_color);
}

Material* Render::find_first_hit(Ray& ray, float* t, glm::vec3* position, glm::vec3* normal) {
	float cur_t = FLT_MAX, t_next;
	Material* cur_material = nullptr;
	glm::vec3 cur_position;
	glm::vec3 cur_normal;

	// Check Sphere Intersection
	for (auto iter = scene->spheres.begin(); iter != scene->spheres.end(); iter++) {
		t_next = intersect_sphere(ray, *iter);
		if (t_next < cur_t) {
			cur_t = t_next;
			cur_position = camera.e + cur_t * ray.dir;
			cur_normal = cur_position - (*iter).center;
			glm::normalize(cur_normal);
			cur_material = &(*iter).material;
		}
	}

	// Check Rectangle Intersection
	for (auto iter = scene->rectangles.begin(); iter != scene->rectangles.end(); iter++) {
		std::pair<float, glm::vec3> temp = intersect_rectangle(ray, *iter);
		t_next = temp.first;
		if (t_next < cur_t) {
			cur_t = t_next;
			cur_position = camera.e + cur_t * ray.dir;
			cur_normal = temp.second;
			cur_material = &(*iter).material;
		}
	}

	// Check AABB Box Intersection
	for (auto iter = scene->boxes.begin(); iter != scene->boxes.end(); iter++) {
		t_next = intersect_box(ray, *iter);
		if (t_next < cur_t) {
			cur_t = t_next;
			cur_position = camera.e + cur_t * ray.dir;

			// Compute the Normal of the intersected Box Face
			glm::vec3 epsilon(1e-4, 1e-4, 1e-4);  // Small offset to determine the box face
			if (glm::abs(cur_position.x - iter->min.x) < epsilon.x) cur_normal = glm::vec3(-1.0f, 0.0f, 0.0f);
			else if (glm::abs(cur_position.x - iter->max.x) < epsilon.x) cur_normal = glm::vec3(1.0f, 0.0f, 0.0f);
			else if (glm::abs(cur_position.y - iter->min.y) < epsilon.y) cur_normal = glm::vec3(0.0f, -1.0f, 0.0f);
			else if (glm::abs(cur_position.y - iter->max.y) < epsilon.y) cur_normal = glm::vec3(0.0f, 1.0f, 0.0f);
			else if (glm::abs(cur_position.z - iter->min.z) < epsilon.z) cur_normal = glm::vec3(0.0f, 0.0f, -1.0f);
			else if (glm::abs(cur_position.z - iter->max.z) < epsilon.z) cur_normal = glm::vec3(0.0f, 0.0f, 1.0f);

			cur_material = &(*iter).material;
		}
	}

	// Check Triangular Mesh Intersection
	// for () {}
	//*t = cur_t;
	//*position = cur_position;
	//*normal = cur_normal;
	return cur_material;
}

float Render::intersect_sphere(Ray& ray, Sphere& sphere) {
	float A, B, C, D, Dsqrt, t;
	glm::vec3 eminusg;

	A = glm::dot(ray.dir, ray.dir);
	eminusg = ray.orig - sphere.center;
	B = 2.0f * glm::dot(ray.dir, eminusg);
	C = glm::dot(eminusg, eminusg) - sphere.radius * sphere.radius;
	D = B * B - 4.0f * A * C;

	if (D < 0.0) return FLT_MAX;
	Dsqrt = sqrt(D);
	if ((t = 0.5 * (-B - Dsqrt) / A) >= 0.0f) return t;
	if ((t = 0.5 * (-B + Dsqrt) / A) >= 0.0f) return t;
	return FLT_MAX;
}

std::pair<float, glm::vec3> Render::intersect_rectangle(Ray& ray, Rectanglee& rectangle) {
	//Prepare Variables
	glm::vec3 e = ray.orig;
	glm::vec3 d = ray.dir;
	glm::vec3 a = rectangle.a;
	glm::vec3 u = rectangle.u;
	glm::vec3 v = rectangle.v;
	glm::vec3 n = glm::normalize(glm::cross(u, v));

	float t;
	if (glm::dot(n, d) == 0.0f) {
		return { FLT_MAX, glm::vec3(0.0f) };	// NO_INTERSECTION
	}
	t = glm::dot(n, (a - e)) / glm::dot(n, d);

	glm::vec3 p;
	if (t < 0.0f) {
		return { FLT_MAX, glm::vec3(0.0f) };	// NO_INTERSECTION
	}
	p = e + t * d;

	float alpha = glm::dot((p - a), u);
	float beta = glm::dot((p - a), v);
	if ((0.0f <= alpha && alpha <= glm::dot(u, u)) &&
		(0.0f <= beta && beta <= glm::dot(v, v))) {
		return { t, n };
	}

	return { FLT_MAX, glm::vec3(0.0f) };	// NO_INTERSECTION
}

float Render::intersect_box(Ray& ray, Box& box) {
	glm::vec3 t_near(- FLT_MAX);
	glm::vec3 t_far(FLT_MAX);
	glm::vec3 e = ray.orig;
	glm::vec3 d = ray.dir;

	// x-axis slab
	float x_min = box.min.x;
	float x_max = box.max.x;
	if (d.x == 0.0f) {
		if (e.x < x_min || e.x > x_max) {
			return FLT_MAX; // NO_INTERSECTION
		}
	}
	else {
		float t1 = (x_min - e.x) / d.x;
		float t2 = (x_max - e.x) / d.x;
		if (t1 > t2) std::swap(t1, t2);
		t_near.x = glm::max(t_near.x, t1);
		t_far.x = glm::min(t_far.x, t2);
		if (t_near.x > t_far.x || t_far.x < 0) return FLT_MAX; // NO_INTERSECTION
	}

	// y-axis slab
	float y_min = box.min.y;
	float y_max = box.max.y;
	if (d.y == 0.0f) {
		if (e.y < y_min || e.y > y_max) {
			return FLT_MAX; // NO_INTERSECTION
		}
	}
	else {
		float t1 = (y_min - e.y) / d.y;
		float t2 = (y_max - e.y) / d.y;
		if (t1 > t2) std::swap(t1, t2);
		t_near.y = glm::max(t_near.y, t1);
		t_far.y = glm::min(t_far.y, t2);
		if (t_near.y > t_far.y || t_far.y < 0) return FLT_MAX; // NO_INTERSECTION
	}

	// z-axis slab
	float z_min = box.min.z;
	float z_max = box.max.z;
	if (d.z == 0.0f) {
		if (e.z < z_min || e.z > z_max) {
			return FLT_MAX; // NO_INTERSECTION
		}
	}
	else {
		float t1 = (z_min - e.z) / d.z;
		float t2 = (z_max - e.z) / d.z;
		if (t1 > t2) std::swap(t1, t2);
		t_near.z = glm::max(t_near.z, t1);
		t_far.z = glm::min(t_far.z, t2);
		if (t_near.z > t_far.z || t_far.z < 0) return FLT_MAX; // NO_INTERSECTION
	}

	//glm::vec3 p = { e.x + t * d.x,  }

	return FLT_MAX;
}

//float Render::intersect_triangular_mesh(Ray& ray, TriangularMesh& box) {

//}

void Render::shade_first_hit(glm::vec3& pos, glm::vec3& norm, Material* material, glm::vec3* shaded_color) {
	//*shaded_color = 0.5f * norm + 0.5f;
	// Initialize the shaded color with the ambient component
	glm::vec3 ambient = material->ka;
	*shaded_color = ambient;

	// Light source properties (you can adjust these as needed)
	glm::vec3 light_pos = glm::vec3(5.0f, 10.0f, 5.0f);
	glm::vec3 light_color = glm::vec3(1.0f, 1.0f, 1.0f);

	// Calculate vectors needed for shading
	glm::vec3 light_dir = glm::normalize(light_pos - pos);  // Direction from the point to the light
	glm::vec3 view_dir = glm::normalize(camera.e - pos);     // Direction from the point to the camera
	glm::vec3 reflect_dir = glm::reflect(-light_dir, norm);  // Reflection of light direction around the normal

	// Diffuse shading (Lambertian reflection)
	float diff_intensity = glm::max(glm::dot(norm, light_dir), 0.0f);
	glm::vec3 diffuse = diff_intensity * material->kd * light_color;

	// Specular shading (Phong reflection)
	float spec_intensity = glm::pow(glm::max(glm::dot(view_dir, reflect_dir), 0.0f), material->n);
	glm::vec3 specular = spec_intensity * material->ks * light_color;

	// Sum the components to get the final color
	*shaded_color += diffuse + specular;

	// Clamp the color values to ensure they remain within [0,1]
	*shaded_color = glm::clamp(*shaded_color, 0.0f, 1.0f);
}

void Render::dump_pixmap() {
	FILE* fp;

	if ((fp = fopen(OUTPUT_FILE_NAME, "wb")) == NULL) {
		fprintf(stderr, "^^^Error: can not open the output file %s...\n", OUTPUT_FILE_NAME);
		exit(-1);
	}
	fprintf(fp, "P6\n%d %d\n255\n", window.width, window.height);
	fwrite(window.pixmap, window.width * window.height * 3, 1, fp);
	fclose(fp);
}