/***********************************************************
Starter code for Assignment 3

This code was originally written by Jack Wang for
CSC418, SPRING 2005

***********************************************************/


#include "raytracer.h"
#include <cstdlib>

int main1(int argc, char* argv[])
{
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	 int width = 600;
	 int height = 600;

	//int width = 600;
	//int height = 600;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, 0, 1);
	Vector3D view(0, 0, -1);
	Vector3D up(0, 1, 0);
	double fov = 60;

	// Defines a material for shading.
	Material glass(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.6, 0.6, 0.6), 80.0, 1, 1.33);
	Material jade(Colour(0, 0, 0), Colour(0.54, 0.89, 0.63),
		Colour(0.316228, 0.316228, 0.316228),12.8, 0.0, 0.0);

	// Defines a point light source.

	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5),
		Colour(0.9, 0.9, 0.9)));

	//raytracer.addLightSource(new PointLight(Point3D(0, 12, 5),
	//	Colour(0.9, 0.9, 0.9)));

	//raytracer.addLightSource(new PointLight(Point3D(-20, 15, 5),
	//	Colour(0.9, 0.9, 0.9)));

	// raytracer.addLightSource(new PointLight(Point3D(0, -12, 5),
	// 	Colour(0.9, 0.9, 0.9)));	

	



	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &glass);
	SceneDagNode* plane = raytracer.addObject(new UnitSquare(), &jade);

	// Apply some transformations to the unit square.
	double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(0, 0, -5));
	raytracer.rotate(sphere, 'x', -45);
	raytracer.rotate(sphere, 'z', 45);
	raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane, Vector3D(0, 0, -7));
	raytracer.rotate(plane, 'z', 45);
	raytracer.scale(plane, Point3D(0, 0, 0), factor2);

	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");

	// Render it from a different point of view.
	Point3D eye2(4, 2, 1);
	Vector3D view2(-4, -2, -6);
	//raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");

	return 0;
}


//another scene
int main(int argc, char* argv[])
{
	// Build your scene and setup your camera here, by calling 
	// functions from Raytracer.  The code here sets up an example
	// scene and renders it from two different view points, DO NOT
	// change this if you're just implementing part one of the 
	// assignment.  
	Raytracer raytracer;
	int width = 600;
	int height = 600;

	//int width = 600;
	//int height = 600;

	if (argc == 3) {
		width = atoi(argv[1]);
		height = atoi(argv[2]);
	}

	// Camera parameters.
	Point3D eye(0, -3, 2);
	Vector3D view(0, 1, 0);
	Vector3D up(0, 0, 1);
	double fov = 90;

	// Defines a material for shading.
	Material glass(Colour(0.0, 0.0, 0.0), Colour(0.0, 0.0, 0.0),
		Colour(0.6, 0.6, 0.6), 80.0, 1, 1.33);
	Material green(Colour(0, 0, 0), Colour(0, 0.89, 0),
		Colour(0,0,0), 12., 0.0, 0.0);
	Material blue(Colour(0, 0, 0), Colour(0, 0, 0.89),
		Colour(0, 0, 0), 12.8, 0.0, 0.0);
	Material red(Colour(0, 0, 0), Colour(0.89, 0, 0),
		Colour(0, 0, 0), 12.8, 0.0, 0.0);
	Material white(Colour(0, 0, 0), Colour(0.89, 0.89, 0.89),
		Colour(0, 0, 0), 12.8, 0.0, 0.0);
	// Defines a point light source.

	raytracer.addLightSource(new PointLight(Point3D(0, 0, 5),
		Colour(0.9, 0.9, 0.9)));

	//raytracer.addLightSource(new PointLight(Point3D(0, 12, 5),
	//	Colour(0.9, 0.9, 0.9)));

	//raytracer.addLightSource(new PointLight(Point3D(-20, 15, 5),
	//	Colour(0.9, 0.9, 0.9)));

	// raytracer.addLightSource(new PointLight(Point3D(0, -12, 5),
	// 	Colour(0.9, 0.9, 0.9)));	





	// Add a unit square into the scene with material mat.
	SceneDagNode* sphere = raytracer.addObject(new UnitSphere(), &glass);
	SceneDagNode* plane_back = raytracer.addObject(new UnitSquare(), &green);
	SceneDagNode* plane_left = raytracer.addObject(new UnitSquare(), &red);
	SceneDagNode* plane_right = raytracer.addObject(new UnitSquare(), &blue);
	SceneDagNode* plane_bot = raytracer.addObject(new UnitSquare(), &white);
	// Apply some transformations to the unit square.
	//double factor1[3] = { 1.0, 1.0, 1.0 };
	double factor2[3] = { 6.0, 6.0, 6.0 };
	raytracer.translate(sphere, Vector3D(-2, 2, 1));
	//raytracer.rotate(sphere, 'x', -45);
	//raytracer.rotate(sphere, 'z', 45);
	//raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

	raytracer.translate(plane_bot, Vector3D(0, 0, 0));
	//raytracer.rotate(plane_bot, 'z', 45);
	raytracer.scale(plane_bot, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_left, Vector3D(-3, 0, 3));
	raytracer.rotate(plane_left, 'x', 90);
	raytracer.rotate(plane_left, 'y', 90);
	raytracer.scale(plane_left, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_right, Vector3D(3, 0, 3));
	raytracer.rotate(plane_right, 'x', -90);
	raytracer.rotate(plane_right, 'y', -90);
	raytracer.scale(plane_right, Point3D(0, 0, 0), factor2);

	raytracer.translate(plane_back, Vector3D(0, 3, 1));
	raytracer.rotate(plane_back, 'x', 90);
	//raytracer.rotate(plane_back, 'y', 90);
	raytracer.scale(plane_back, Point3D(0, 0, 0), factor2);


	// Render the scene, feel free to make the image smaller for
	// testing purposes.	
	raytracer.render(width, height, eye, view, up, fov, "view1.bmp");

	// Render it from a different point of view.
	//Point3D eye2(4, 2, 1);
	//Vector3D view2(-4, -2, -6);
	//raytracer.render(width, height, eye2, view2, up, fov, "view2.bmp");

	return 0;
}