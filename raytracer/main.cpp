/***********************************************************
Starter code for Assignment 3

This code was originally written by Jack Wang for
CSC418, SPRING 2005

***********************************************************/


#include "raytracer.h"
#include <cstdlib>
#include <iostream>
#include <chrono>
#include <ctime>




unsigned char* ReadBMP(char* filename)
{
	static unsigned char *texels;
	static int width, height;

	FILE *fd;
	fd = fopen(filename, "rb");
	if (fd == NULL)
	{
		printf("Error: fopen failed\n");
		return NULL;
	}

	unsigned char header[54];

	// Read header
	fread(header, sizeof(unsigned char), 54, fd);

	// Capture dimensions
	width = *(int*)&header[18];
	height = *(int*)&header[22];

	int padding = 0;

	// Calculate padding
	while ((width * 3 + padding) % 4 != 0)
	{
		padding++;
	}

	// Compute new width, which includes padding
	int widthnew = width * 3 + padding;

	// Allocate memory to store image data (non-padded)
	texels = (unsigned char *)malloc(width * height * 3 * sizeof(unsigned char));
	if (texels == NULL)
	{
		printf("Error: Malloc failed\n");
		return NULL;
	}

	// Allocate temporary memory to read widthnew size of data
	unsigned char* data = (unsigned char *)malloc(widthnew * sizeof(unsigned int));

	// Read row by row of data and remove padded data.
	for (int i = 0; i<height; i++)
	{
		// Read widthnew length of data
		fread(data, sizeof(unsigned char), widthnew, fd);

		// Retain width length of data, and swizzle RB component.
		// BMP stores in BGR format, my usecase needs RGB format
		for (int j = 0; j < width * 3; j += 3)
		{
			int index = (i * width * 3) + (j);
			texels[index + 0] = data[j + 2];
			texels[index + 1] = data[j + 1];
			texels[index + 2] = data[j + 0];
		}
	}

	free(data);
	fclose(fd);

	return texels;
}


int main(int argc, char* argv[])
{
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  

	auto start = std::chrono::system_clock::now();
	std::cout << "started";

	Raytracer raytracer;
	 int width = 1920;
	 int height = 1080;

	/*int width = 600;
	int height = 600;*/

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	//read texture
	unsigned char* texture = ReadBMP("tiles.bmp");

	// Defines a material for shading.
	Material gold(Colour(0.3, 0.3, 0.3), Colour(0.75164, 0.60648, 0.22648),
		Colour(0.628281, 0.555802, 0.366065),
		51.2, NULL);
	Material jade(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
		Colour(0.316228, 0.316228, 0.316228),
		12.8, NULL);

	Material orange(Colour(0, 0, 0), Colour(1.0, 0.6, 0.0),
		Colour(0.316228, 0.316228, 0.316228),
		30.0, NULL);

	Material yellow(Colour(0, 0, 0), Colour(0.0, 1.0, 1.0),
		Colour(0.316228, 0.316228, 0.316228),
		30.0, NULL);

	Material blue(Colour(0, 0, 0), Colour(0.0, 0.0, 1.0),
		Colour(0.316228, 0.316228, 0.316228),
		30.0, NULL);



	// Defines a point light source.

	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5),
		Colour(0.9, 0.9, 0.9)));

	/*raytracer.addLightSource(new PointLight(Point3D(0, 12, 5),
		Colour(0.9, 0.9, 0.9)));
*/
	/*raytracer.addLightSource(new PointLight(Point3D(-20, 15, 5),
		Colour(0.9, 0.9, 0.9)));*/

	// raytracer.addLightSource(new PointLight(Point3D(0, -12, 5),
	// 	Colour(0.9, 0.9, 0.9)));	

	



	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &gold);
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &jade);

	//SceneDagNode* sphere_2 = raytracer.addObject(new UnitSphere(), &blue);
	//SceneDagNode* sphere_3 = raytracer.addObject(new UnitSphere(), &yellow);


	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 2.0, 1.0 };
	double factor3[3] = { 0.7, 0.7, 0.7 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));
	raytracer.rotate(sphere, 'x', -45);
	raytracer.rotate(sphere, 'z', 45);
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);


	//raytracer.translate(sphere_2, Vector3D(-1, -1, -3));
	//raytracer.rotate(sphere_2, 'x', -45);
	//raytracer.rotate(sphere_2, 'z', 45);
	//raytracer.scale(sphere_2, Point3D(0, 0, 0), factor3);

	//raytracer.translate(sphere_3, Vector3D(1, 2, -7));
	//raytracer.rotate(sphere_3, 'x', -45);
	//raytracer.rotate(sphere_3, 'z', 45);
	//raytracer.scale(sphere_3, Point3D(0, 0, 0), factor3);

	raytracer.translate(plane, Vector3D(0, 0, -7));
	raytracer.rotate(plane, 'z', 45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");

	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	std::cout << "finished computation at " << std::ctime(&end_time)<< "elapsed time: " << elapsed_seconds.count() << "s\n";
	
 	return 0;
}