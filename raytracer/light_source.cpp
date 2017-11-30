/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

Colour PointLight::shade(Ray3D& ray, PointLight& light) {
	//Assume everything is in world coordinate 
	Point3D lightPos = light.get_position();

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
	current_phong_color[0] = current_phong_color[0] * (light._col_ambient[0]);
	current_phong_color[1] = current_phong_color[1] * (light._col_ambient[1]);
	current_phong_color[2] = current_phong_color[2] * (light._col_ambient[2]);
	current_phong_color.clamp();
	return current_phong_color;
}

void PointLight::shade( Ray3D& ray ) {
	//Assume everything is in world coordinate 
	Point3D lightPos = this->get_position();
	
	Point3D intersectPoint = ray.intersection.point;
	Vector3D normalVector = ray.intersection.normal;
	normalVector.normalize();
	Vector3D viewVector = -ray.dir;
	viewVector.normalize();
	Vector3D light_vector = (lightPos - intersectPoint);
	light_vector.normalize();
	Vector3D difference = (light_vector + viewVector);
	//calculate reflect
	Vector3D reflect = light_vector - 2*(light_vector.dot(normalVector))*normalVector;
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
	current_phong_color[0] = current_phong_color[0]*(this->_col_ambient[0]);
	current_phong_color[1] = current_phong_color[1]*(this->_col_ambient[1]);
	current_phong_color[2] = current_phong_color[2]*(this->_col_ambient[2]);
	
	//ray.col = 0.5*ray.col + 0.5*current_phong_color;
	ray.col = current_phong_color;

	//ray.col = ka*(0.8*ray.intersection.mat->ambient+0.2*this->_col_ambient) + kd*cos_diffuse*(0.8*ray.intersection.mat->diffuse+0.2*this->_col_diffuse) + ks*pow(cos_specular, ray.intersection.mat->specular_exp)*(0.8*ray.intersection.mat->specular+0.2*this->_col_specular);
	//clamp color to 1
	//std::cout << ray.col << "\n";
	ray.col.clamp();

	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

}

void AreaLight::shade(Ray3D & ray)
{
	//radomly generate many point lights in the area
	Point3D light_center = this->get_position();
	double rad = this->get_radius();
	double samples = 1 * rad;
	Colour col = Colour(0, 0, 0);

	for (int i = 0; i < samples; i++) {
		double rand_x = -rad/2 + ((double)rand() / (RAND_MAX)) * rad /2;
		double rand_y = -rad/2 + ((double)rand() / (RAND_MAX)) * rad / 2;
		Point3D light_pos = Point3D(light_center[0] + rand_x, light_center[1] + rand_y, light_center[2]);
		PointLight eachLight = PointLight(light_pos, this->_col_ambient);
		col = col +(1.0/samples) * eachLight.shade(ray, eachLight);
	}

	ray.col = col;

}
