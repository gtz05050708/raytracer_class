#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <utility>
#include <vector>
#include <stack>    
#include <deque>
#include <fstream>
#include <sstream>
#include <algorithm> 
#include <FreeImage.h>
//#include "readfile.h"
#include "Transform.h"
using namespace std;

class pointLight {
public:
	glm::vec3 position;
	glm::vec3 color;
	glm::vec3 attenuation = glm::vec3(0.0f, 0.0f, 1.0f);

	pointLight(glm::vec3 position, glm::vec3 color) {
		this->position = position;
		this->color = color;
	}
};

class directLight {
public:
	glm::vec3 position;
	glm::vec3 color;
	glm::vec3 attenuation = glm::vec3(1.0f, 0.0f, 0.0f);

	directLight(glm::vec3 position, glm::vec3 color) {
		this->position = position;
		this->color = color;
	}
};


class Sphere {
public:
	glm::vec3 center; //center
	float r; //radius

	glm::vec3 emission;
	glm::vec3 ambient = glm::vec3(0.2f, 0.2f, 0.2f);
	glm::vec3 specular;
	glm::vec3 diffuse;
	float shininess;
	glm::mat4 transform;
	// attenuation

	Sphere(glm::vec3 center, float radius) {
		this->center = center;
		this->r = radius;
	}

};

class Ray {
public:
	glm::vec3 p0;
	glm::vec3 d;
	Ray(const glm::vec3 p0, float fovy, float width, float height, float curI, float curJ, glm::vec3 u, glm::vec3 v, glm::vec3 w) {
		float tanfovy = glm::tan(glm::radians(fovy / 2.0f));
		float tanfovx = tanfovy * width / height;
		float alpha = tanfovx * (((curJ + 0.5f) - width / 2.0f) / (width / 2.0f));
		float beta = tanfovy * (((curI + 0.5f) - height / 2.0f) / (height / 2.0f));
		this->p0 = p0;
		d = glm::normalize(alpha * u + beta * v - w);
	}

	Ray(glm::vec3 p0, glm::vec3 d) {
		this->p0 = p0;
		this->d = d;
	}

	glm::vec3 getPoint(float t) {
		return p0 + d * t;
	}
};

class Triangle {
public:
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;

	glm::vec3 emission;
	glm::vec3 ambient = glm::vec3(0.2f, 0.2f, 0.2f);
	glm::vec3 specular;
	glm::vec3 diffuse;
	float shininess;
	glm::mat4 transform;
	// attenuation

	Triangle(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
	}
};

class Trinormal {
public:
	glm::vec3 p1;
	glm::vec3 p2;
	glm::vec3 p3;

	glm::vec3 n1;
	glm::vec3 n2;
	glm::vec3 n3;

	glm::vec3 emission;
	glm::vec3 ambient = glm::vec3(0.2f, 0.2f, 0.2f);
	glm::vec3 specular;
	glm::vec3 diffuse;
	float shininess;
	glm::mat4 transform;
	// attenuation

	Trinormal(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec3 n1, glm::vec3 n2, glm::vec3 n3) {
		this->p1 = p1;
		this->p2 = p2;
		this->p3 = p3;
		this->n1 = n1;
		this->n2 = n2;
		this->n3 = n3;
	}
};

class Camera {
public:
	float fovy;
	glm::vec3 eye = glm::vec3(0, 0, 1);
	glm::vec3 center = glm::vec3(0, 0, 0);
	glm::vec3 up = glm::vec3(0, 1, 0);
	glm::vec3 u;
	glm::vec3 v;
	glm::vec3 w;

	Camera(glm::vec3 eye, glm::vec3 center, glm::vec3 up, float fovy) {
		this->eye = eye;
		this->center = center;
		this->up = up;
		this->w = glm::normalize(eye - center);
		this->u = glm::normalize(glm::cross(up, w));
		this->v = glm::normalize(glm::cross(w, u));
		this->fovy = fovy;
	}

	void getUVW() {
		this->w = glm::normalize(eye - center);
		this->u = glm::normalize(glm::cross(up, w));
		this->v = glm::cross(w, u);
	}
};

class Scene {
public:
	int numobjects = 0;
	int maxobjects = 10;
	string outputName;

	// General
	// size
	float width;
	float height;
	// maxdepth
	int maxdepth = 5;

	// Camera
	vector<Camera> cameras;

	// Geometry
	vector<Sphere> allSpheres;
	vector<Triangle> allTriangles;
	vector<Trinormal> allTrinormals;
	// implement triwithnormals

	// light
	static const int numLights = 10;
	int numPoint = 0;                     // How many point light are used 
	int numDirection = 0; // How many directional light are used 

	//vector<glm::vec3> pointPos;
	//vector<glm::vec3> pointColor;
	vector<pointLight> allPointLights;
	vector<directLight> allDirectLights;
	//vector<glm::vec3> directionaltPos;
	//vector<glm::vec3> directionalColor;

	float c0; // att
	float c1; // att
	float c2; // att

	int maxverts;
	int maxvertnorms;

	int numverts = 0;
	int numvertnorms = 0;

	vector<glm::vec3> vertices;
	vector<glm::vec3> vertexnormals;
	vector<glm::vec3> vertexnoramlsn;
};


void matransform(stack<glm::mat4> &transfstack, float *values)
{
	glm::mat4 transform = transfstack.top();
	glm::vec4 valvec = glm::vec4(values[0], values[1], values[2], values[3]);
	glm::vec4 newval = transform * valvec;
	for (int i = 0; i < 4; i++) values[i] = newval[i];
}

void rightmultiply(const glm::mat4 & M, stack<glm::mat4> &transfstack)
{
	glm::mat4 &T = transfstack.top();
	T = T * M;
}

// Function to read the input data values
// Use is optional, but should be very helpful in parsing.  
bool readvals(stringstream &s, const int numvals, float* values)
{
	for (int i = 0; i < numvals; i++) {
		s >> values[i];
		if (s.fail()) {
			cout << "Failed reading value " << i << " will skip\n";
			return false;
		}
	}
	return true;
}

void readfile(const char* filename, Scene &scene)
{
	string str, cmd;
	ifstream in;
	in.open(filename);

	if (in.is_open()) {

		// I need to implement a matrix stack to store transforms.  
		// This is done using standard STL Templates 
		stack <glm::mat4> transfstack;
		transfstack.push(glm::mat4(1.0));  // identity

		getline(in, str);

		float ambient[3] = { 0.0f,0.0f,0.0f};// local
		float diffuse[3] = { 0.0f,0.0f,0.0f};// local
		float specular[3] = { 0.0f,0.0f,0.0f};// local
		float emission[3] = {0.0f,0.0f,0.0f};// local
		float attenuation[3] = { 1.0f, 0.0f, 0.0f }; // local
		float shininess = 0.0f;// local

		while (in) {
			if ((str.find_first_not_of(" \t\r\n") != string::npos) && (str[0] != '#')) {
				// Ruled out comment and blank lines 

				stringstream s(str);
				s >> cmd;
				int i;
				float values[10]; 
				bool validinput; // Validity of input 

				// Process the light, add it to database.
				// Lighting Command
				if (cmd == "point") {
					if (scene.numPoint + scene.numDirection == scene.numLights) { // No more Lights 
						cerr << "Reached Maximum Number of Lights " << (scene.numPoint + scene.numDirection) << " Will ignore further lights\n";
					}
					else {
						validinput = readvals(s, 6, values); // Position/color for lts.
						if (validinput) {

							// YOUR CODE FOR HW 2 HERE. 
							// Note that values[0...7] shows the read in values 
							// Make use of lightposn[] and lightcolor[] arrays in variables.h
							// Those arrays can then be used in display too.  
							glm::vec3 lightpos = glm::vec3(values[0], values[1], values[2]);
							glm::vec3 lightcolor = glm::vec3(values[3], values[4], values[5]);
							pointLight point(lightpos, lightcolor);
							point.attenuation = glm::vec3(attenuation[0], attenuation[1], attenuation[2]);
							scene.allPointLights.push_back(point);
							++scene.numPoint;
						}
					}
				}

				else if (cmd == "directional") {
					if (scene.numPoint + scene.numDirection == scene.numLights) { // No more Lights 
						cerr << "Reached Maximum Number of Lights " << (scene.numPoint + scene.numDirection) << " Will ignore further lights\n";
					}
					else {
						validinput = readvals(s, 6, values); // Position/color for lts.
						if (validinput) {

							// YOUR CODE FOR HW 2 HERE. 
							// Note that values[0...7] shows the read in values 
							// Make use of lightposn[] and lightcolor[] arrays in variables.h
							// Those arrays can then be used in display too.  
							glm::vec3 lightpos = glm::vec3(values[0], values[1], values[2]);
							glm::vec3 lightcolor = glm::vec3(values[3], values[4], values[5]);
							directLight direct(lightpos, lightcolor);
							direct.attenuation = glm::vec3(1.0f, 0.0f, 0.0f);
							scene.allDirectLights.push_back(direct);
							++scene.numDirection;
						}
					}
				}

				else if (cmd == "attenuation") {
					validinput = readvals(s, 3, values); // attenuation 
					if (validinput) {
						for (i = 0; i < 3; i++) {
							attenuation[i] = values[i];
						}
					}
				}

				// Material Commands 
				// Ambient, diffuse, specular, shininess properties for each object.
				// Filling this in is pretty straightforward, so I've left it in 
				// the skeleton, also as a hint of how to do the more complex ones.
				// Note that no transforms/stacks are applied to the colors. 

				else if (cmd == "ambient") {
					validinput = readvals(s, 3, values); // colors 
					if (validinput) {
						for (i = 0; i < 3; i++) {
							ambient[i] = values[i];
						}
					}
				}
				else if (cmd == "diffuse") {
					validinput = readvals(s, 3, values);
					if (validinput) {
						for (i = 0; i < 3; i++) {
							diffuse[i] = values[i];
						}
					}
				}
				else if (cmd == "specular") {
					validinput = readvals(s, 3, values);
					if (validinput) {
						for (i = 0; i < 3; i++) {
							specular[i] = values[i];
						}
					}
				}
				else if (cmd == "emission") {
					validinput = readvals(s, 3, values);
					if (validinput) {
						for (i = 0; i < 3; i++) {
							emission[i] = values[i];
						}
					}
				}
				else if (cmd == "shininess") {
					validinput = readvals(s, 1, values);
					if (validinput) {
						shininess = values[0];
					}
				}
				else if (cmd == "size") {
					validinput = readvals(s, 2, values);
					if (validinput) {
						scene.width = (int)values[0]; scene.height = (int)values[1];
					}
				}
				else if (cmd == "output") {
					s >> scene.outputName;
				}
				else if (cmd == "maxdepth") {
					validinput = readvals(s, 1, values);
					if (validinput) {
						scene.maxdepth = (int)values[0];
					}
				}
				else if (cmd == "camera") {
					validinput = readvals(s, 10, values); // 10 values eye cen up fov
					if (validinput) {

						// YOUR CODE FOR HW 2 HERE
						// Use all of values[0...9]
						// You may need to use the upvector fn in Transform.cpp
						// to set up correctly. 
						// Set eyeinit upinit center fovy in variables.h 
						float fovy = values[9];
						glm::vec3 eyeinit = glm::vec3(values[0], values[1], values[2]);
						glm::vec3 center = glm::vec3(values[3], values[4], values[5]);
						glm::vec3 upinit = glm::vec3(values[6], values[7], values[8]);
						Camera cam(eyeinit, center, upinit, fovy);
						scene.cameras.push_back(cam);
					}
				}

				else if (cmd == "sphere") {
					validinput = readvals(s, 4, values);
					if (validinput) {
						glm::vec3 center = glm::vec3(values[0], values[1], values[2]);
						Sphere sphere(center, values[3]);
						for (int i = 0; i < 3; i++) {
							(sphere.ambient)[i] = ambient[i];
							(sphere.diffuse)[i] = diffuse[i];
							(sphere.specular)[i] = specular[i];
							(sphere.emission)[i] = emission[i];
						}
						sphere.shininess = shininess;
						sphere.transform = transfstack.top();
						scene.allSpheres.push_back(sphere);
					}
				}

				else if (cmd == "maxverts") {
					validinput = readvals(s, 1, values);
					scene.maxverts = values[0];
				}

				else if (cmd == "maxvertnorms") {
					validinput = readvals(s, 1, values);
					scene.maxvertnorms = values[0];
				}

				else if (cmd == "vertex") {
					validinput = readvals(s, 3, values);
					scene.vertices.push_back(glm::vec3(values[0], values[1], values[2]));
					scene.numverts++;
				}

				else if (cmd == "vertexnormal") {
					validinput = readvals(s, 6, values);
					scene.vertexnormals.push_back(glm::vec3(values[0], values[1], values[2]));
					scene.vertexnoramlsn.push_back(glm::vec3(values[3], values[4], values[5]));
					scene.numvertnorms++;
				}

				else if (cmd == "tri") {
					validinput = readvals(s, 3, values);
					float v1 = values[0];
					float v2 = values[1];
					float v3 = values[2];
					Triangle tri(scene.vertices[v1], scene.vertices[v2], scene.vertices[v3]);
					for (int i = 0; i < 3; i++) {
						(tri.ambient)[i] = ambient[i];
						(tri.diffuse)[i] = diffuse[i];
						(tri.specular)[i] = specular[i];
						(tri.emission)[i] = emission[i];
					}
					tri.shininess = shininess;
					tri.transform = transfstack.top();
					scene.allTriangles.push_back(tri);
				}

				else if (cmd == "trinormal") {
					validinput = readvals(s, 3, values);
					float v1 = values[0];
					float v2 = values[1];
					float v3 = values[2];
					Trinormal tri(scene.vertexnormals[v1], scene.vertexnormals[v2], scene.vertexnormals[v3],
						scene.vertexnoramlsn[v1], scene.vertexnoramlsn[v2], scene.vertexnoramlsn[v3]);
					for (int i = 0; i < 3; i++) {
						(tri.ambient)[i] = ambient[i];
						(tri.diffuse)[i] = diffuse[i];
						(tri.specular)[i] = specular[i];
						(tri.emission)[i] = emission[i];
					}
					tri.shininess = shininess;
					tri.transform = transfstack.top();
					scene.allTrinormals.push_back(tri);
				}
				/*
				// I've left the code for loading objects in the skeleton, so
				// you can get a sense of how this works.
				// Also look at demo.txt to get a sense of why things are done this way.
				else if (cmd == "sphere" || cmd == "cube" || cmd == "teapot") {
					if (scene.numobjects == scene.maxobjects) { // No more objects
						cerr << "Reached Maximum Number of Objects " << scene.numobjects << " Will ignore further objects\n";
					} else {
						validinput = readvals(s, 1, values);
						if (validinput) {
							object * obj = &(objects[numobjects]);
							obj->size = values[0];

							// Set the object's light properties
							for (i = 0; i < 4; i++) {
								(obj->ambient)[i] = ambient[i];
								(obj->diffuse)[i] = diffuse[i];
								(obj->specular)[i] = specular[i];
								(obj->emission)[i] = emission[i];
							}
							obj->shininess = shininess;

							// Set the object's transform
							obj->transform = transfstack.top();

							// Set the object's type
							if (cmd == "sphere") {
								obj->type = sphere;
							} else if (cmd == "cube") {
								obj->type = cube;
							} else if (cmd == "teapot") {
								obj->type = teapot;
							}
						}
						++numobjects;
					}
				}
				*/
				else if (cmd == "translate") {
					validinput = readvals(s, 3, values);
					if (validinput) {
						// YOUR CODE FOR HW 2 HERE.  
						// Think about how the transformation stack is affected
						// You might want to use helper functions on top of file. 
						// Also keep in mind what order your matrix is!
						glm::mat4 ret = Transform::translate(values[0], values[1], values[2]);
						rightmultiply(ret, transfstack);
					}
				}
				else if (cmd == "scale") {
					validinput = readvals(s, 3, values);
					if (validinput) {

						// YOUR CODE FOR HW 2 HERE.  
						// Think about how the transformation stack is affected
						// You might want to use helper functions on top of file.  
						// Also keep in mind what order your matrix is!
						mat4 ret = Transform::scale(values[0], values[1], values[2]);
						rightmultiply(ret, transfstack);
					}
				}
				else if (cmd == "rotate") {
					validinput = readvals(s, 4, values);
					if (validinput) {

						// YOUR CODE FOR HW 2 HERE. 
						// values[0..2] are the axis, values[3] is the angle.  
						// You may want to normalize the axis (or in Transform::rotate)
						// See how the stack is affected, as above.  
						// Note that rotate returns a mat3. 
						// Also keep in mind what order your matrix is!
						vec3 axis = vec3(values[0], values[1], values[2]);
						mat4 ret = mat4(Transform::rotate(values[3], axis));
						rightmultiply(ret, transfstack);
					}
				}

				// I include the basic push/pop code for matrix stacks
				else if (cmd == "pushTransform") {
					transfstack.push(transfstack.top());
				}
				else if (cmd == "popTransform") {
					if (transfstack.size() <= 1) {
						cerr << "Stack has no elements.  Cannot Pop\n";
					}
					else {
						transfstack.pop();
					}
				}

				else {
					cerr << "Unknown Command: " << cmd << " Skipping \n";
				}
			}
			getline(in, str);
		}

		// Set up initial position for eye, up and amount
		// As well as booleans 

	}
	else {
		cerr << "Unable to Open Input Data File " << filename << "\n";
		throw 2;
	}
}
