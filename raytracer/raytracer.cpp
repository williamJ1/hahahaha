/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		Implementations of functions in raytracer.h, 
		and the main function which specifies the 
		scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

Raytracer::Raytracer() : _lightSource(NULL) {
	_root = new SceneDagNode();
}

Raytracer::~Raytracer() {
	delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
		SceneObject* obj, Material* mat ) {
	SceneDagNode* node = new SceneDagNode( obj, mat );
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return node;;
}

void Raytracer::addObject_tree( SceneDagNode* parent, 
		SceneObject* obj, Material* mat, Matrix4x4 modelToWorld, Point3D BB_max, Point3D BB_min) {
	SceneDagNode* node = new SceneDagNode( obj, mat, modelToWorld, BB_max, BB_min);
	node->parent = parent;
	node->next = NULL;
	node->child = NULL;
	
	// Add the object to the parent's child list, this means
	// whatever transformation applied to the parent will also
	// be applied to the child.
	if (parent->child == NULL) {
		parent->child = node;
	}
	else {
		parent = parent->child;
		while (parent->next != NULL) {
			parent = parent->next;
		}
		parent->next = node;
	}
	
	return;
}




LightListNode* Raytracer::addLightSource( LightSource* light ) {
	LightListNode* tmp = _lightSource;
	_lightSource = new LightListNode( light, tmp );
	return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
		    node->trans = node->trans*rotation; 	
			angle = -angle;
		} 
		else {
			node->invtrans = rotation*node->invtrans; 
		}	
	}
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	node->trans = node->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	node->trans = node->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
		Vector3D up ) {
	Matrix4x4 mat; 
	Vector3D w;
	view.normalize();
	up = up - up.dot(view)*view;
	up.normalize();
	w = view.cross(up);

	mat[0][0] = w[0];
	mat[1][0] = w[1];
	mat[2][0] = w[2];
	mat[0][1] = up[0];
	mat[1][1] = up[1];
	mat[2][1] = up[2];
	mat[0][2] = -view[0];
	mat[1][2] = -view[1];
	mat[2][2] = -view[2];
	mat[0][3] = eye[0];
	mat[1][3] = eye[1];
	mat[2][3] = eye[2];

	return mat; 
}


void Raytracer::computeTransforms( SceneDagNode* node )
{
    SceneDagNode *childPtr;
    if (node->parent != NULL )
    {
		node->modelToWorld = node->parent->modelToWorld*node->trans; 
		//modelToWorld is the transform from current child to world, trans is the transform from child to parent
        node->worldToModel = node->invtrans*node->parent->worldToModel; 
    }
    else
    {
        node->modelToWorld = node->trans;
        node->worldToModel = node->invtrans; 
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        computeTransforms(childPtr);
        childPtr = childPtr->next;
    }



}

//check if current ray intersect object in scene
void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    SceneDagNode *childPtr;

    // Applies transformation of the current node to the global
    // transformation matrices.
    if (node->obj) {
        // Perform intersection.
        if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
            ray.intersection.mat = node->mat;
        }
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        traverseScene(childPtr, ray);
        childPtr = childPtr->next;
    }

}

void Raytracer::traverseScene_BSP( AABB_node* tree_root, Ray3D& ray ) {
    SceneDagNode *childPtr;
	Matrix4x4 BB_worldToModel = Matrix4x4();
	Matrix4x4 BB_modelToWorld = Matrix4x4();


	bool ray_intersect_BB = tree_root->cube_obj->intersect(ray, BB_worldToModel, BB_modelToWorld);
	//check with current bounding box for intersection

	//if intersect
		//check if current tree node is leaf 
		//if its leaf, calculate intersection
		//else call traverseScene_BSP with left/right children

	//else return
}

void Raytracer::computeShading( Ray3D& ray, int* phong_count) {
    LightListNode* curLight = _lightSource;
	LightListNode* curLight_temp = _lightSource;

	double light_num = 0;

	for (;;){
		if(curLight_temp == NULL) break;
		light_num ++;
   		curLight_temp = curLight_temp->next;
	}
	//std::cout << light_num << "\n";
	double light_ratio = 1/light_num;
	//std::cout <<"ratio   " << light_ratio << "\n";

    for (;;) {
        if (curLight == NULL) break;
        // Each lightSource provides its own shading function.
        // Implement shadows here if needed.
		//====================================
		//shoot a new ray
		Vector3D ray_dir = curLight->light->get_position() - ray.intersection.point;
		ray_dir.normalize();
		Ray3D shadow_ray = Ray3D((ray.intersection.point + 0.01 * ray.intersection.normal), ray_dir);
		traverseScene(_root, shadow_ray);
		//double t_light = (curLight->light->get_position() - ray.intersection.point)/ray_dir;
		Colour temp = Colour(0, 0, 0);
		if (shadow_ray.intersection.none){
			//TODO:check if the object is between light and start point
			temp = ray.col;
			curLight->light->shade(ray);
			*phong_count = *phong_count + 1;
			ray.col = light_ratio*ray.col + (1-light_ratio)*temp;
			//ray.col = 0.6*ray.col + 0.4*temp;
			//std::cout << ray.col << "\n";
		}
		else {
			ray.col = ray.col + Colour(0, 0, 0);
		}
		//====================================
		// Colour temp = ray.col;
		// for (int i = -1; i < 1; i++){
		// 	Vector3D ray_dir = curLight->light->get_position() - ray.intersection.point;
		// 	ray_dir[0] = ray_dir[0] + i;
		// 	ray_dir[1] = ray_dir[1] + i;
		// 	ray_dir[2] = ray_dir[2] + i;
		// 	ray_dir.normalize();
		// 	Ray3D shadow_ray = Ray3D((ray.intersection.point + 0.01 * ray_dir), ray_dir);
		// 	traverseScene(_root, shadow_ray);
		// 	curLight->light->shade(ray);
		// 	if(shadow_ray.intersection.none){
		// 		temp = temp + 0.2*ray.col;
		// 	}
			
		// }
		// ray.col = temp;
		//std::cout << temp << "\n";
		//====================================

		//curLight->light->shade(ray);
        curLight = curLight->next;
    }
}

void Raytracer::initPixelBuffer() {
    int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
    _rbuffer = new unsigned char[numbytes];
    std::fill_n(_rbuffer, numbytes,0);
    _gbuffer = new unsigned char[numbytes];
    std::fill_n(_gbuffer, numbytes,0);
    _bbuffer = new unsigned char[numbytes];
    std::fill_n(_bbuffer, numbytes,0);
}

void Raytracer::flushPixelBuffer( std::string file_name ) {
    bmp_write( file_name.c_str(), _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
    delete _rbuffer;
    delete _gbuffer;
    delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray , int depth, int d_end, Colour col) {
	if (depth >= d_end){
		//std::cout << "combine color" << col << "\n";
		return col;
	}
	//_root is scene graph
    traverseScene(_root, ray); 
    // Don't bother shading if the ray didn't hit 
    // anything.
	//ray.col = Colour(0,0,0);
	int phong_count = 0;

    if (!ray.intersection.none) {
		// if (depth == 1) {
		// 	//std::cout << "wrong" << "\n";
		// 	//std::cout << "ray_origin" << ray.origin << "\n";
		// 	//std::cout << "ray direction" << ray.dir << "\n";
		// 	//traverseScene(_root, ray);
		// 	//return Colour(0,0,0.9);
		// }
        computeShading(ray, &phong_count); 
		//std::cout << "count" << phong_count<<"\n";
		// ray.col[0] = ray.col[0] / phong_count;
		// ray.col[1] = ray.col[1] / phong_count;
		// ray.col[2] = ray.col[2] / phong_count;
		//if (depth == 1) {
		//	if (ray.intersection.point[2] > -6.9) {
		//		if (col[0] != 0 || col[1] != 0 || col[2] != 0){
		//			if (ray.col[0] != 0 || ray.col[1] != 0 || ray.col[2] != 0) {
		//				std::cout << "prev color" << 255 * col << "\n";
		//				std::cout << "cur color" << 255 * ray.col << "\n";
		//			}
		//		}
		//	}
		//}

		Vector3D ray_dir = ray.dir;
		ray_dir.normalize();
		Vector3D ray_n = ray.intersection.normal;
		ray_n.normalize();
		Vector3D ref_dir = ray_dir - 2 * (ray_dir.dot(ray_n))*ray_n;
		ref_dir.normalize();
		
		Ray3D ref_ray = Ray3D((ray.intersection.point + 0.0001 * ray_n), ref_dir);
		//std::cout << "combine color" << col << "\n";
		if (depth == 0) {
			return shadeRay(ref_ray, depth + 1, d_end, col + ray.col);
		}
		else {
			return shadeRay(ref_ray, depth + 1, d_end, col + 0.1 * ray.col);
		}
    }
	else {
		//std::cout << "combine color" << col << "\n";
		return col;
	}

    // You'll want to call shadeRay recursively (with a different ray, 
    // of course) here to implement reflection/refraction effects.  
	//calculate reflection
	//Vector3D ray_dir = ray.dir;
	//ray_dir.normalize();
	//Vector3D ray_n = ray.intersection.normal;
	//ray_n.normalize();
	//Vector3D ref_dir = ray_dir - 2*(ray.dir.dot(ray_n))*ray_n;
	//Ray3D ref_ray = Ray3D(ray.intersection.point, ref_dir);
	//return shadeRay(ref_ray, depth, d_end, col);
    
}	

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
        Vector3D up, double fov, std::string fileName ) {
    computeTransforms(_root);
	computeBB(_root);
	AABB_node *tree_root =new AABB_node();
	tree_root->cube_obj = new UnitCube();
	buildAABBtree(_root, tree_root);
    Matrix4x4 viewToWorld;
    _scrWidth = width;
    _scrHeight = height;
    double factor = (double(height)/2)/tan(fov*M_PI/360.0);

    initPixelBuffer();
    viewToWorld = initInvViewMatrix(eye, view, up);

    // Construct a ray for each pixel.
	// double largest_color_value = 0.0;
    for (int i = 0; i < _scrHeight; i++) {
        for (int j = 0; j < _scrWidth; j++) {
            // Sets up ray origin and direction in view space, 
            // image plane is at z = -1.
            Point3D origin(0, 0, 0);
			Point3D imagePlane;
			imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
			imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
			imagePlane[2] = -1;

			// TODO: Convert ray to world space and call 
			// shadeRay(ray) to generate pixel colour. 	
			Colour init_color = Colour(0, 0, 0);
			Colour col = Colour(0, 0, 0);
			int SA_times = 1;
			double bias = (1.0 / factor) / SA_times;
			for (int k = 0; k < SA_times; k++) {
				for (int l = 0; l < SA_times; l++) {
					Vector3D RayDirection = Vector3D(imagePlane[0] + bias*k, imagePlane[1] + bias*l, imagePlane[2]);
					Ray3D ray = Ray3D(viewToWorld*origin, viewToWorld*RayDirection);
					int d_end = 1;
					col = col + (1.0/std::pow((double)SA_times,2)) * shadeRay(ray, 0, d_end, init_color);
				}
			}

			col.clamp();
			//std::cout << "after" << col << "\n"; 
			// if (col[0] >1 || col[1]>1 || col[2]>1){
			// 	col = (1.0 / (double)d_end)*col;
			// }
			// double current_max = fmax(col[0], col[1]);
			// current_max = fmax(current_max, col[2]);
			// largest_color_value = fmax(largest_color_value, current_max);

				
			_rbuffer[i*width+j] = int(col[0]*255);
			_gbuffer[i*width+j] = int(col[1]*255);
			_bbuffer[i*width+j] = int(col[2]*255);

			
		}
	}




	flushPixelBuffer(fileName);
}


void Raytracer::computeBB( SceneDagNode* node)
{
	//traverse all scene objects and find bounding box for all
	if (node->obj != NULL){
		node->BB_max = node->obj->BBmax(node->modelToWorld);
		node->BB_min = node->obj->BBmin(node->modelToWorld);
	}
	// std::cout << "BB_MAX:  " << node->BB_max << "\n";
	// std::cout << "BB_MIN:  " << node->BB_min << "\n";

	SceneDagNode * childPtr = node->child;
	while (childPtr != NULL){
		computeBB(childPtr);
		childPtr = childPtr->next;
	}
}


//build AABB tree
//find current bounding box 
void Raytracer::buildAABBtree( SceneDagNode* node_list, AABB_node* tree_node)
{
	//assume all SceneDagNode has no first term, first term is node_list->next
	//Case1: node_list has 1 item, no need to continue partition
	SceneDagNode* child;
	if (node_list->child != NULL) {
		child = node_list->child;
	}
	else {
		return;
	}

	if (child->next == NULL){
		tree_node->scene_obj = child;
		// std::cout << "inhere3" << "\n";
		return;
	}
	//Cse2: node_list has several items, do BSP 
	//assume can go through all nodes using node_list->next
	//find current bounding box
	Point3D cur_min = Point3D(999.0, 999.0, 999.0);
	Point3D cur_max = Point3D(-999.0, -999.0, -999.0);
	SceneDagNode* temp_node = node_list->child;
	while(temp_node != NULL){
		cur_min[0] = fmin(temp_node->BB_min[0], cur_min[0]);
		cur_min[1] = fmin(temp_node->BB_min[1], cur_min[1]);
		cur_min[2] = fmin(temp_node->BB_min[2], cur_min[2]);

		cur_max[0] = fmax(temp_node->BB_max[0], cur_max[0]);
		cur_max[1] = fmax(temp_node->BB_max[1], cur_max[1]);
		cur_max[2] = fmax(temp_node->BB_max[2], cur_max[2]);

		temp_node = temp_node->next;
	}

	//set current bounding box size 
	tree_node->BB_max = cur_max;
	tree_node->BB_min = cur_min;

	//partition current bounding box and seperate objects according to the plain
	//always split along x axis
	//find normal of the plain, and a point on plain where the normal shoot out:
	double delta_z = (cur_max[2] - cur_min[2])/2;
	Vector3D normal_vec = Vector3D(0, 0, -delta_z);
	Point3D int_point = Point3D(cur_min[0], cur_max[1], ((cur_max[2]) - delta_z));
    SceneDagNode *left = new SceneDagNode();
	SceneDagNode *right = new SceneDagNode();
	//traverse through list of nodes and seperate into left and right
	SceneDagNode *current_node = node_list->child;
	//SceneDagNode o = *current_node;
	while(current_node != NULL){
		//find object center
		Point3D object_center = current_node->modelToWorld * Point3D(0, 0, 0);
		Vector3D object_dir = object_center - int_point;
		double cos_value = normal_vec.dot(object_dir);
		if (cos_value <= 0){
			addObject_tree(left, current_node->obj, current_node->mat, current_node->modelToWorld, current_node->BB_max, current_node->BB_min);
			
		}
		else{
			addObject_tree(right, current_node->obj, current_node->mat, current_node->modelToWorld, current_node->BB_max, current_node->BB_min);
		}

		current_node = current_node->next;
	}

	AABB_node* left_tree_node = new AABB_node();
	AABB_node* right_tree_node = new AABB_node();
	left_tree_node->parent = tree_node;
	right_tree_node->parent = tree_node;

	left_tree_node->cube_obj = new UnitCube();
	right_tree_node->cube_obj = new UnitCube();

	tree_node->left = left_tree_node;
	tree_node->right = right_tree_node;

	buildAABBtree(left, left_tree_node);
	buildAABBtree(right, right_tree_node);

	//should not reach here
	return;
}