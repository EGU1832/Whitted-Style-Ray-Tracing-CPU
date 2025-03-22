//
// CSE6449 Real-Time Rendering
// Sogang Univiersity
//
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <opencv2/opencv.hpp>

#include "RayTracing.h"
#include "measure_cpu_time.h"

int main() {
	int image_width = 1200, image_height = 800;

	Camera camera(glm::vec3(0.0f, 0.0f, 3.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f),
		glm::vec3(-0.2f, 0.0f, 1.0f), 60.0f, (float) image_width / image_height);
	Window window(image_width, image_height, camera);
	Scene* scene = new Scene;
	scene->initialize();

	Render render(scene, camera, window);

	CHECK_TIME_START(_start, _freq);
	render.render();
	CHECK_TIME_END(_start, _end, _freq, _compute_time);
	fprintf(stdout, "*** Ray tracing time = %.3fms\n", _compute_time);

	// USER CODE (2025.03.21)
	// Read result.ppm redered Image
	cv::Mat image = cv::imread("result.ppm", cv::IMREAD_COLOR);
	if (image.empty()) {
		fprintf(stdout, "Error loading PPM file!\n");
		return -1;
	}

	// Print Image to Window
	cv::imshow("Rendered Image", image);
	cv::waitKey(0); // Wait for Key Input to Exit

}