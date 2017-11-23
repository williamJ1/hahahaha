/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	//transform ray into object space
	Point3D ray_origin = worldToModel * (ray.origin);
	Vector3D ray_dir = worldToModel * (ray.dir);

	//init vertices of the unit squre
	Point3D vt[4];
	vt[0] = Point3D(0.5, 0.5, 0);
	vt[1] = Point3D(-0.5, 0.5, 0);
	vt[2] = Point3D(-0.5, -0.5, 0);
	vt[3] = Point3D(0.5, -0.5, 0);
	
	//get normal
	Vector3D vec1 = vt[1] - vt[0];
	Vector3D vec2 = vt[2] - vt[0];
	Vector3D n = vec1.cross(vec2);
	//square plain normal
	n.normalize();
	
	//check if dir dot n = 0
	if (ray_dir.dot(n) == 0){
		//ray.intersection.none = true;
		return false;
	}
	
	//compute t
	double t = (vt[0] - ray_origin).dot(n) / ray_dir.dot(n);

	//check if point lies inside unit square
	Point3D inter_p = ray_origin + t * ray_dir;
	if (inter_p[0] < -0.5 || inter_p[0] > 0.5 || inter_p[1] < -0.5 || inter_p[1] > 0.5){
		//outside unit square
		//ray.intersection.none = true;
		return false;
	}

	//valid intersection 
	Point3D world_inter_p = modelToWorld * inter_p;
	Vector3D world_n = modelToWorld * n;
	if (ray.intersection.none == true) {
		//no previous intersection
		ray.intersection.none = false;
		ray.intersection.normal = world_n;
		ray.intersection.point = world_inter_p;
		ray.intersection.t_value = t;
		//std::cout << ray.intersection.point << "\n";
		return true;
	}
	else if (ray.intersection.t_value >= t) {
		//found a closer intersection
		//replace
		ray.intersection.normal = world_n;
		ray.intersection.point = world_inter_p;
		ray.intersection.t_value = t;
		//std::cout << ray.intersection.point << "\n";
		return true;
	}
	else {
		//prev intersection is closer
		//do not replace
		return false;
	}


	//return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	//transform ray into model space
	Point3D ray_origin = worldToModel * (ray.origin);
	Vector3D ray_dir = worldToModel * (ray.dir);
	ray_dir.normalize();
	Point3D sphere_center = Point3D(0,0,0);

	
	double t0, t1; // solutions for t if the ray intersects 
	Vector3D L = sphere_center - ray_origin;
	double tca = L.dot(ray_dir);
	//if sphere lines on the other side of the origin
	if (tca < 0){
		//ray.intersection.none = true;
		return false;
	} 
	double d2 = L.dot(L) - tca * tca;
	//if d bigger than radius, ray does not intersect sphere
	if (d2 > 1){
		//ray.intersection.none = true;
		return false;
	}
	double thc = sqrt(1 - d2);
	t0 = tca - thc;
	t1 = tca + thc;
	
	//find the smaller value of t0 and t1, whic represents the cloest intersection point
	if (t0 > t1) std::swap(t0, t1); 
 
    if (t0 < 0) { 
        t0 = t1; // if t0 is negative, let's use t1 instead 
        if (t0 < 0) {
			//ray.intersection.none = true;
			return false; // both t0 and t1 are negative
		} 
    } 

	double previous_t_value = ray.intersection.t_value;
	//check if current object is in front of previous stored object
	if (previous_t_value > t0){
		ray.intersection.t_value = t0; 
		ray.intersection.none = false;
		Point3D inter_p = Point3D(ray_origin + t0*ray_dir);
		Vector3D n = inter_p - sphere_center;
		n.normalize();
		ray.intersection.normal = modelToWorld * n;
		ray.intersection.point = modelToWorld * inter_p;
	}

	std::cout << "orgin" << ray_origin << "\n";
	std::cout << ray.intersection.point << "\n";

	return true;
	
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	//ray.intersection.none = true;
	//return false;
}

