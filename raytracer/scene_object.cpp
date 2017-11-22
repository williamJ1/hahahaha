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
	n = n.normalize()*n;
	
	//check if dir.n = 0
	if (ray_dir.dot(n) == 0){
		ray.intersection.none = true;
		return false;
	}
	
	//compute t
	double t= (vt[0] - ray_origin).dot(n) / ray_dir.dot(n);

	//check if point lies inside unit square
	Point3D inter_p = ray_origin + t * ray_dir;
	if (inter_p[0] < -0.5 || inter_p[1] > 0.5 || inter_p[1] < -0.5 || inter_p[1] > 0.5){
		//outside unit square
		ray.intersection.none = true;
		return false;
	}

	//valid intersection 
	//transform back to model coord
	ray.intersection.normal =modelToWorld * n;
	ray.intersection.point = modelToWorld * inter_p;
	ray.intersection.t_value = t;

	
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	
	return false;
}

