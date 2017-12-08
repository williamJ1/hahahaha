/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		   light source classes

***********************************************************/

#include "util.h"

Colour CalculatePhong(Ray3D& ray, Point3D& lightPos, Colour& ambient);

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
	virtual void shade( Ray3D& ) = 0;
	virtual Point3D get_position() const = 0; 
	virtual double get_radius() const = 0;
	virtual Colour get_ambient_light() const = 0;
	virtual int get_light_type() const = 0;
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
	PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
	_col_diffuse(col), _col_specular(col) {}
	PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
	: _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
	_col_specular(specular) {}
	void shade( Ray3D& ray );
	Point3D get_position() const { return _pos; }
	double get_radius() const { return 0; }
	Colour get_ambient_light() const {return _col_ambient; }
	int get_light_type() const {return 0;}

private:
	Point3D _pos;
	Colour _col_ambient;
	Colour _col_diffuse; 
	Colour _col_specular; 
};



// A area light is defined by its center radius and color
class AreaLight : public LightSource {
public:
	AreaLight(Point3D center, double radius, Colour col) : _pos(center),_radius(radius), _col_ambient(col),
		_col_diffuse(col), _col_specular(col) {}
	AreaLight(Point3D center, double radius, Colour ambient, Colour diffuse, Colour specular)
		: _pos(center), _radius(radius), _col_ambient(ambient), _col_diffuse(diffuse),
		_col_specular(specular) {}
	void shade(Ray3D& ray);
	Point3D get_position() const { return _pos; }
	double get_radius() const { return _radius; }
	Colour get_ambient_light() const {return _col_ambient; }
	int get_light_type() const {return 1;}
private:
	Point3D _pos;
	double _radius;
	Colour _col_ambient;
	Colour _col_diffuse;
	Colour _col_specular;
};
