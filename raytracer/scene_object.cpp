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
	t = std::abs(t);
	//check if point lies inside unit square
	Point3D inter_p = ray_origin + t * ray_dir;
	if (inter_p[0] < -0.5 || inter_p[0] > 0.5 || inter_p[1] < -0.5 || inter_p[1] > 0.5 || (inter_p[2]<-0.00001 || inter_p[2] > 0.00001)){
		//outside unit square
		//ray.intersection.none = true;
		return false;
	}

	//valid intersection 
	Point3D world_inter_p = modelToWorld * inter_p;
	Vector3D world_n = worldToModel.transpose() * n;
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
		return false;
	} 
	double d2 = L.dot(L) - tca * tca;
	//if d bigger than radius, ray does not intersect sphere
	if (d2 > 1){
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
			return false; // both t0 and t1 are negative
		} 
    } 

	double previous_t_value = ray.intersection.t_value;
	//check if current object is in front of previous stored object
	if (ray.intersection.none == true || previous_t_value > t0){
		ray.intersection.t_value = t0; 
		ray.intersection.none = false;
		Point3D inter_p = Point3D(ray_origin + t0*ray_dir);
		Vector3D n = inter_p - sphere_center;
		n.normalize();
		ray.intersection.normal = worldToModel.transpose() * n;
		ray.intersection.point = modelToWorld * inter_p;
		return true;
	}
	// std::cout << "orgin" << ray_origin << "\n";
	// std::cout << ray.intersection.point << "\n";
	return false;
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.
	//ray.intersection.none = true;
	//return false;
}

bool UnitCube::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld ) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
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
	t = std::abs(t);
	//check if point lies inside unit square
	Point3D inter_p = ray_origin + t * ray_dir;
	if (inter_p[0] < -0.5 || inter_p[0] > 0.5 || inter_p[1] < -0.5 || inter_p[1] > 0.5 || (inter_p[2]<-0.00001 || inter_p[2] > 0.00001)){
		//outside unit square
		//ray.intersection.none = true;
		return false;
	}

	//valid intersection 
	Point3D world_inter_p = modelToWorld * inter_p;
	Vector3D world_n = worldToModel.transpose() * n;
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

Point3D UnitSphere::BBmin(const Matrix4x4& modelToWorld){
	Matrix4x4 s = Matrix4x4();
	s.set_value(15, -1.0);
	//s = [ 1  0  0  0 ]
	//[ 0  1  0  0 ]
	//[ 0  0  1  0 ]
	//[ 0  0  0 -1 ]
	//Note: in this case inverse of s is the same
	Matrix4x4 m_transpose = modelToWorld.transpose();
	Matrix4x4 r = modelToWorld*s*m_transpose;
	// let (M S^-1 M^t) = R = [ r11 r12 r13 r14 ]  (note that R is symmetric: R=R^t)
	// 			              [ r12 r22 r23 r24 ]
	//    			          [ r13 r23 r33 r34 ]
	//     			          [ r14 r24 r34 r44 ]
	Point3D minxyz = Point3D();
	//z
	minxyz[2] = (r.get_value(11) + sqrt(r.get_value(11)*r.get_value(11) - r.get_value(15)*r.get_value(10)))/ r.get_value(15);
	//y
	minxyz[1] = (r.get_value(7) + sqrt(r.get_value(7)*r.get_value(7) - r.get_value(15)*r.get_value(5))) / r.get_value(15);
	//x
	minxyz[0] = (r.get_value(3) + sqrt(r.get_value(3)*r.get_value(3) - r.get_value(15)*r.get_value(0))) / r.get_value(15);
	// z =  (r34 +/- sqrt(r34^2 - r44 r33) ) / r44
	// y = (r24 +/- sqrt(r24^2 - r44 r22) ) / r44
	// x = (r14 +/- sqrt(r14^2 - r44 r11) ) / r44

	return minxyz;
}

Point3D UnitSphere::BBmax(const Matrix4x4& modelToWorld){
	Matrix4x4 s = Matrix4x4();
	s.set_value(15, -1.0);
	//s = [ 1  0  0  0 ]
	//[ 0  1  0  0 ]
	//[ 0  0  1  0 ]
	//[ 0  0  0 -1 ]
	//Note: in this case inverse of s is the same
	Matrix4x4 m_transpose = modelToWorld.transpose();
	Matrix4x4 r = modelToWorld*s*m_transpose;
	// let (M S^-1 M^t) = R = [ r11 r12 r13 r14 ]  (note that R is symmetric: R=R^t)
	// 			              [ r12 r22 r23 r24 ]
	//    			          [ r13 r23 r33 r34 ]
	//     			          [ r14 r24 r34 r44 ]
	Point3D maxxyz = Point3D();
	//z
	maxxyz[2] = (r.get_value(11) - sqrt(r.get_value(11)*r.get_value(11) - r.get_value(15)*r.get_value(10)))/ r.get_value(15);
	//y
	maxxyz[1] = (r.get_value(7) - sqrt(r.get_value(7)*r.get_value(7) - r.get_value(15)*r.get_value(5))) / r.get_value(15);
	//x
	maxxyz[0] = (r.get_value(3) - sqrt(r.get_value(3)*r.get_value(3) - r.get_value(15)*r.get_value(0))) / r.get_value(15);
	// z =  (r34 +/- sqrt(r34^2 - r44 r33) ) / r44
	// y = (r24 +/- sqrt(r24^2 - r44 r22) ) / r44
	// x = (r14 +/- sqrt(r14^2 - r44 r11) ) / r44

	return maxxyz;

}


Point3D UnitSquare::BBmin(const Matrix4x4& modelToWorld){
	Point3D p1 = modelToWorld * Point3D(-0.5, -0.5, 0);
	Point3D p2 = modelToWorld * Point3D(-0.5, 0.5, 0);
	Point3D p3 = modelToWorld * Point3D(0.5, -0.5, 0);
	Point3D p4 = modelToWorld * Point3D(0.5, 0.5, 0);

	Point3D minxyz = Point3D(0, 0, 0);
	//x
	minxyz[0] = fmin(p1[0], fmin(p2[0], fmin(p3[0], p4[0])));
	//y
	minxyz[1] = fmin(p1[1], fmin(p2[1], fmin(p3[1], p4[1])));
	//z
	minxyz[2] = fmin(p1[2], fmin(p2[2], fmin(p3[2], p4[2])));

	return minxyz;
}



Point3D UnitSquare::BBmax(const Matrix4x4& modelToWorld){
	Point3D p1 = modelToWorld * Point3D(-0.5, -0.5, 0);
	Point3D p2 = modelToWorld * Point3D(-0.5, 0.5, 0);
	Point3D p3 = modelToWorld * Point3D(0.5, -0.5, 0);
	Point3D p4 = modelToWorld * Point3D(0.5, 0.5, 0);
	// std::cout<< "p1" << p1 << "\n";
	// std::cout<< "p2" << p2 << "\n";
	// std::cout<< "p3" << p3 << "\n";
	// std::cout<< "p4" << p4 << "\n";

	Point3D maxxyz = Point3D(0, 0, 0);
	//x
	maxxyz[0] = fmax(p1[0], fmax(p2[0], fmax(p3[0], p4[0])));
	//y
	maxxyz[1] = fmax(p1[1], fmax(p2[1], fmax(p3[1], p4[1])));
	//z
	maxxyz[2] = fmax(p1[2], fmax(p2[2], fmax(p3[2], p4[2])));

	return maxxyz; 
}

Point3D UnitCube::BBmin(const Matrix4x4& modelToWorld){
	Point3D minxyz = Point3D(0, 0, 0);
	return minxyz;
}



Point3D UnitCube::BBmax(const Matrix4x4& modelToWorld){
	Point3D maxxyz = Point3D(0, 0, 0);
	return maxxyz; 
}