/***********************************************************
Starter code for Assignment 3

This code was originally written by Jack Wang for
CSC418, SPRING 2005

implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

Colour CalculatePhong(Ray3D& ray, Point3D& lightPos, Colour& ambient) {
	//Assume everything is in world coordinate 

	Point3D intersectPoint = ray.intersection.point;
	Vector3D normalVector = ray.intersection.normal;
	normalVector.normalize();
	Vector3D viewVector = -ray.dir;
	viewVector.normalize();
	Vector3D light_vector = (lightPos - intersectPoint);
	light_vector.normalize();
	Vector3D difference = (light_vector + viewVector);
	//calculate reflect
	Vector3D reflect = light_vector - 2 * (light_vector.dot(normalVector))*normalVector;
	reflect.normalize();
	difference.normalize();
	//double cos_specular = fmax(difference.dot(normalVector), 0.0);
	double cos_specular = fmax(reflect.dot(-viewVector), 0.0);
	//std::cout << cos_specular << "\n";
	double cos_diffuse = fmax(light_vector.dot(normalVector), 0.0);

	double ka = 1;
	double kd = 1;
	double ks = 1;
	//double ks = 3 - ka - kd;

	Colour current_phong_color = (ka*ray.intersection.mat->ambient + kd*cos_diffuse*ray.intersection.mat->diffuse + ks*pow(cos_specular, ray.intersection.mat->specular_exp)*ray.intersection.mat->specular);
	current_phong_color[0] = current_phong_color[0] * (ambient[0]);
	current_phong_color[1] = current_phong_color[1] * (ambient[1]);
	current_phong_color[2] = current_phong_color[2] * (ambient[2]);
	current_phong_color.clamp();
	return current_phong_color;
}



// TODO: implement this function to fill in values for ray.col 
// using phong shading.  Make sure your vectors are normalized, and
// clamp colour values to 1.0.

//
// It is assumed at this point that the intersection information in ray 
// is available.  So be sure that traverseScene() is called on the ray 
// before this function.  
void PointLight::shade(Ray3D& ray) {
	if (ray.intersection.mat->text_array != NULL) {
		unsigned char* texture = ray.intersection.mat->text_array;
		//normal of intersection point and sphere
		Vector3D nor = ray.intersection.normal;
		nor.normalize();

		double u = 0.5 + atan2(nor[2], nor[0]) / (2 * 3.141592653589793);
		double v = 0.5 - asin(nor[1]) / 3.141592653589793;

		int i = (int)(u * 256);//105 is image width
		int j = (int)(v * 256);
		unsigned int r;
		unsigned int g;
		unsigned int b;

		r = texture[j * 256 + i];
		g = texture[j * 256 + i + 1];
		b = texture[j * 256 + i + 2];

		double rf = r / 255.0;
		double gf = g / 255.0;
		double bf = b / 255.0;
		ray.col = Colour(rf, gf, bf);
	}
	else {
		//Assume everything is in world coordinate 
		Point3D lightPos = this->get_position();
		Colour ambient = this->get_ambient_light();
		Colour phongColor = CalculatePhong(ray, lightPos, ambient);
		ray.col = phongColor;
	}

}





void AreaLight::shade(Ray3D & ray)
{
}
