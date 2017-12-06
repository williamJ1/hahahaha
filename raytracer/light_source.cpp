/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade( Ray3D& ray ) {

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
		//std::cout << (int)texture[50] << (int)texture[51] << (int)texture[52];
		//std::cout << i << " " << j << "\n";
		r = texture[j * 256 + i];
		g = texture[j * 256 + i + 1];
		b = texture[j * 256 + i + 2];
		//std::cout << r << " " << g << " " << b << "\n";

		//for (int i = 0; i < 256; i++) {
		//	for (int j = 0; j < 256; j++) {
		//		std::cout << (int)texture[j * 256 + i] << " " << (int)texture[j * 256 + i + 1] << " " << (int)texture[j * 256 + i + 2] << "\n";
		//	}
		//}
		double rf = r / 255.0;
		double gf = g / 255.0;
		double bf = b / 255.0;



		ray.col = Colour(rf, gf, bf);
	}
	else {
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

		// if (ray.col[0] > 1) ray.col[0] = 1;
		// if (ray.col[1] > 1) ray.col[1] = 1;	
		// if (ray.col[2] > 1) ray.col[2] = 1;



		//double r = (double)rand() / (RAND_MAX) + 1;
		//double g = (double)rand() / (RAND_MAX)+1;
		//double b = (double)rand() / (RAND_MAX)+1;

		//ray.col = Colour(r, g, b);
	

		// TODO: implement this function to fill in values for ray.col 
		// using phong shading.  Make sure your vectors are normalized, and
		// clamp colour values to 1.0.
	
		//
		// It is assumed at this point that the intersection information in ray 
		// is available.  So be sure that traverseScene() is called on the ray 
		// before this function.  
	}


}

///////////sample
//  //find the input light vector
//  vec3 light_vector = normalize(lightPos - vertPos);
//  //find the viewing_vector
//  vec3 eye_vec = normalize(vec3(0,0,0) - vertPos);
//  //find the reflect vector
//  vec3 reflect_vec = normalize(reflect(-light_vector, normalInterp));
//
//  //calculate crossproduct of eye vector and input light vector
//  vec3 difference = normalize(light_vector + eye_vec);
//
//  //calculate specular Color
//  //I have implemented both way to find specular cos value
//  //float cos_specular = max(dot(difference, normalInterp), 0.0);
//  float cos_specular = max(dot(reflect_vec, eye_vec), 0.0);
//
//  //calculate diffuse Color
//  float cos_diffuse = max(dot(light_vector, normalInterp), 0.0);
//  //float cos_diffuse = max(dot(light_vector, normalInterp), 0.0);
//
//
//
//  gl_FragColor = vec4(Ka*ambientColor+Kd*diffuseColor*cos_diffuse+Ks*specularColor*(pow(cos_specular, shininessVal)), 1.0);