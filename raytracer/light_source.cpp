/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {
	//Assume everything is in world coordinate 
	Point3D lightPos = this->get_position();
	Point3D intersectPoint = ray.intersection.point;
	Vector3D normalVector = ray.intersection.normal;
	Vector3D viewVector = ray.dir.normalize()*ray.dir;
	Vector3D light_vector = (lightPos - intersectPoint).normalize()*(lightPos - intersectPoint);
	Vector3D difference = (light_vector + viewVector).normalize()*(light_vector + viewVector);
	double cos_specular = fmax(difference.dot(normalVector), 0.0);
	double cos_diffuse = fmax(light_vector.dot(normalVector), 0.0);

	double ka = 0.3;
	double kd = 0.3;
	double ks = 1 - ka - kd;

	ray.col = ka*ray.intersection.mat->ambient + kd*cos_diffuse*ray.intersection.mat->diffuse + ks*pow(cos_specular, ray.intersection.mat->specular_exp)*ray.intersection.mat->specular;

	

	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

}

/////////sample
  //find the input light vector
  vec3 light_vector = normalize(lightPos - vertPos);
  //find the viewing_vector
  vec3 eye_vec = normalize(vec3(0,0,0) - vertPos);
  //find the reflect vector
  vec3 reflect_vec = normalize(reflect(-light_vector, normalInterp));

  //calculate crossproduct of eye vector and input light vector
  vec3 difference = normalize(light_vector + eye_vec);

  //calculate specular Color
  //I have implemented both way to find specular cos value
  //float cos_specular = max(dot(difference, normalInterp), 0.0);
  float cos_specular = max(dot(reflect_vec, eye_vec), 0.0);

  //calculate diffuse Color
  float cos_diffuse = max(dot(light_vector, normalInterp), 0.0);
  //float cos_diffuse = max(dot(light_vector, normalInterp), 0.0);



  gl_FragColor = vec4(Ka*ambientColor+Kd*diffuseColor*cos_diffuse+Ks*specularColor*(pow(cos_specular, shininessVal)), 1.0);