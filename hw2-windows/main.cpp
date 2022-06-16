/******************************************************************************/
/* This is the program skeleton for homework 2 in CSE167 by Ravi Ramamoorthi  */
/* Extends HW 1 to deal with shading, more transforms and multiple objects    */
/******************************************************************************/

// You shouldn't have to edit this file at all!

#include "hw4.h"
#include <omp.h>

float interceptSphereTest(Ray& ray, Sphere& sphere) {
	float b = glm::dot(ray.d, (ray.p0 - sphere.center)) * 2.0f;
	float bsquare = glm::pow(b, 2);
	float ac4 = 4.0f * glm::pow(glm::length(ray.d), 2) * (glm::pow(glm::length(ray.p0 - sphere.center), 2) - glm::pow(sphere.r, 2));
	float discriminate = bsquare - ac4;
	if (discriminate < 0) {
		return -99.0f;
	}
	float t1 = (-1.0f * b + glm::sqrt(discriminate)) / (2.0f * glm::pow(glm::length(ray.d), 2));
	float t2 = (-1.0f * b - glm::sqrt(discriminate)) / (2.0f * glm::pow(glm::length(ray.d), 2));
	float t = t1;
	if (t1 >= 0 && t2 >= 0) {
		t = t1 < t2 ? t1 : t2;
	}
	else if (t1 >= 0 && t2 < 0) {
		t = t1;
	}
	else if (t2 >= 0 && t1 < 0) {
		t = t2;
	}
	else {
		return -99.0f;
	}

	return t;
}

float interceptTriangleTest(Ray& ray, Triangle& tri) {
	glm::vec3 p1 = tri.p1;
	glm::vec3 p2 = tri.p2;
	glm::vec3 p3 = tri.p3;
	glm::vec3 d = ray.d * -1.0f;
	glm::vec4 p0Vec = glm::vec4(ray.p0[0], ray.p0[1], ray.p0[2], 1.0f);
	glm::mat4 inverse = glm::mat4(p1.x, p1.y, p1.z, 1.0f,
		p2.x, p2.y, p2.z, 1.0f,
		p3.x, p3.y, p3.z, 1.0f,
		d.x, d.y, d.z, 0.0f);
	inverse = glm::inverse(inverse);
	glm::vec4 result = inverse * p0Vec;
	if (result[0] >= 0.0f && result[1] >= 0.0f && result[2] >= 0.0f && result[3] >= 0.0f) {
		return result[3];
	}
	else {
		return -99.0f;
	}
}

pair<glm::vec3, int> getInterceptSphere(Ray& ray, Scene& scene) {

	glm::vec3 point(-99.0f, -99.0f, -99.0f);
	glm::vec3 normal(-99.0f, -99.0f, -99.0f);
	float smalldepth = -99.0f;
	bool foundIntersect = false;
	int index = 0;
	int num = -99;
	for (Sphere& sphere : scene.allSpheres) {
		glm::vec4 transp0 = glm::vec4(ray.p0, 1.0f);
		glm::vec4 transd = glm::vec4(ray.d, 0.0f);
		glm::mat4 inverseMat = glm::inverse(sphere.transform);
		transp0 = inverseMat * transp0;
		transp0 = transp0 / transp0.w;
		transd = inverseMat * transd;
		glm::vec3 newp0 = glm::vec3(transp0);
		glm::vec3 newd = glm::vec3(glm::normalize(transd));
		Ray transRay(newp0, newd);
		float t = interceptSphereTest(transRay, sphere);
		if (t != -99.0f) {
			glm::vec3 intersectionLocal = transRay.d * (t - 0.01f) + transRay.p0;
			glm::vec4 intersection4DLocal = glm::vec4(intersectionLocal, 1.0f);
			glm::vec4 intersection4DWorld = sphere.transform * intersection4DLocal;
			intersection4DWorld = intersection4DWorld / intersection4DWorld.w;
			glm::vec3 intersetcionWorld = glm::vec3(intersection4DWorld);
			float worldZ = glm::length(intersetcionWorld - ray.p0);
			if (!foundIntersect || worldZ < smalldepth) {
				point = intersetcionWorld;
				foundIntersect = true;
				smalldepth = worldZ;
				num = index;
			}
		}
		index++;
	}
	pair<glm::vec3, int> thePair = make_pair(point, num);
	return thePair;
}

pair<glm::vec3, int> getInterceptTriangle(Ray& ray, Scene& scene) {

	glm::vec3 output(-99.0f, -99.0f, -99.0f);
	float smallT = 0;
	bool foundIntersect = false;
	int num = -99;
	int index = 0;
	for (Triangle& tri : scene.allTriangles) {
		float t = interceptTriangleTest(ray, tri);
		if (t != -99.0f) {
			if (!foundIntersect || t < smallT) {
				glm::vec3 intersection3D = ray.d * (t - 0.01f) + ray.p0;
				smallT = t;
				foundIntersect = true;
				output = intersection3D;
				num = index;
			}
		}
		index++;
	}
	pair<glm::vec3, int> thePair = make_pair(output, num);
	return thePair;
}

void saveScreenshot(string fname, int width, int height, BYTE* pixels) {

	FIBITMAP* img = FreeImage_ConvertFromRawBits(pixels, width, height, width * 3, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);

	std::cout << "Saving screenshot: " << fname << "\n";

	FreeImage_Save(FIF_PNG, img, fname.c_str(), 0);
}

glm::vec3 ComputeLight(const glm::vec3 direction, const glm::vec3 lightcolor, const glm::vec3  normal, const glm::vec3 halfvec,
	const glm::vec3 mydiffuse, const glm::vec3 myspecular, const float myshininess) {

	float nDotL = glm::dot(normal, direction);
	glm::vec3 lambert = mydiffuse * lightcolor * max(nDotL, 0.0f);

	float nDotH = glm::dot(normal, halfvec);
	glm::vec3 phong = myspecular * lightcolor * pow(max(nDotH, 0.0f), myshininess);

	glm::vec3 retval = lambert + phong;
	return retval;
}

glm::vec3 getColor(Scene& scene, int type, int index, int depth, glm::vec3 intersect, glm::vec3 eye, Ray& ray) {
	if (depth == 0) {
		return glm::vec3(0, 0, 0);
	}
	glm::vec3 normal(0, 0, 0);
	glm::vec3 emission(0, 0, 0);
	glm::vec3 ambient(0, 0, 0);
	glm::vec3 specular(0, 0, 0);
	glm::vec3 diffuse(0, 0, 0);
	float shininess = 0;
	// object intersected is a sphere
	if (type == 0) {
		glm::vec4 intersectionWorld4D = glm::vec4(intersect, 1.0f);
		glm::vec4 intersectionLocal4D = glm::inverse(scene.allSpheres[index].transform) * intersectionWorld4D;
		intersectionLocal4D = intersectionLocal4D / intersectionLocal4D.w;
		glm::vec3 intersectionLocal = glm::vec3(intersectionLocal4D);

		glm::vec3 normalLocal = glm::normalize(intersectionLocal - scene.allSpheres[index].center);
		glm::mat3 transformMat = glm::mat3(scene.allSpheres[index].transform);
		transformMat = glm::inverse(glm::transpose(transformMat));
		normal = glm::normalize(transformMat * normalLocal);

		emission = scene.allSpheres[index].emission;
		ambient = scene.allSpheres[index].ambient;
		specular = scene.allSpheres[index].specular;
		diffuse = scene.allSpheres[index].diffuse;
		shininess = scene.allSpheres[index].shininess;
	}
	else {
		glm::vec3 side1 = scene.allTriangles[index].p3 - scene.allTriangles[index].p1;
		glm::vec3 side2 = scene.allTriangles[index].p3 - scene.allTriangles[index].p2;
		normal = glm::normalize(glm::cross(side1, side2));
		emission = scene.allTriangles[index].emission;
		ambient = scene.allTriangles[index].ambient;
		specular = scene.allTriangles[index].specular;
		diffuse = scene.allTriangles[index].diffuse;
		shininess = scene.allTriangles[index].shininess;
	}


	// now check whether the light is visible 
	glm::vec3 colorlvl = emission + ambient;

	for (pointLight& point : scene.allPointLights) {
		glm::vec3 lightDirn = glm::normalize(point.position - intersect);
		Ray lightRay(intersect, lightDirn);
		pair<glm::vec3, int> spherePair = getInterceptSphere(lightRay, scene);
		pair<glm::vec3, int> trianglePair = getInterceptTriangle(lightRay, scene);
		// if the light is visible, compute the color
		bool visible = false;
		float sphereZ = glm::length(spherePair.first - intersect);
		float triZ = glm::length(trianglePair.first - intersect);
		float lightZ = glm::length(point.position - intersect);
		if (spherePair.second == -99 && trianglePair.second == -99) {
			visible = true;
		}
		else if (spherePair.second != -99 && trianglePair.second == -99) {
			if (sphereZ > lightZ) {
				visible = true;
			}
		}
		else if (spherePair.second == -99 && trianglePair.second != -99) {
			if (triZ > lightZ) {
				visible = true;
			}
		}
		else {
			if (sphereZ > lightZ && triZ > lightZ) {
				visible = true;
			}
		}

		if (visible) {
			// get half angle for point light
			glm::vec3 eyedirn = glm::normalize(eye - intersect);
			glm::vec3 halfAngle = glm::normalize(lightDirn + eyedirn);
			glm::vec3 color = ComputeLight(lightDirn, point.color, normal, halfAngle, diffuse, specular, shininess);
			float distance = glm::length(lightDirn);
			distance = point.attenuation.x + point.attenuation.y * distance + point.attenuation.z * distance * distance;
			color = color / distance;
			colorlvl += color;
		}
	}

	for (directLight& direct : scene.allDirectLights) {
		glm::vec3 lightDirn = glm::normalize(direct.position);

		Ray lightRay(intersect, lightDirn);
		pair<glm::vec3, int> spherePair = getInterceptSphere(lightRay, scene);
		pair<glm::vec3, int> trianglePair = getInterceptTriangle(lightRay, scene);
		// if the light is visible, compute the color
		bool visible = false;
		if (spherePair.second == -99 && trianglePair.second == -99) {
			visible = true;
		}
		if (visible) {
			// get half angle for point light
			glm::vec3 eyedirn = glm::normalize(eye - intersect);
			glm::vec3 halfAngle = glm::normalize(lightDirn + eyedirn);
			glm::vec3 color = ComputeLight(lightDirn, direct.color, normal, halfAngle, diffuse, specular, shininess);
			float distance = direct.attenuation.x;
			color /= distance;
			colorlvl += color;
		}
	}

	glm::vec3 reflected = glm::normalize(ray.d + 2.0f * glm::dot(-1.0f * ray.d, normal) * normal);
	Ray reflectRay(intersect, reflected);
	pair<glm::vec3, int> sphereIntersect = getInterceptSphere(reflectRay, scene);
	pair<glm::vec3, int> triIntersect = getInterceptTriangle(reflectRay, scene);

	if (sphereIntersect.second != -99 || triIntersect.second != -99) {
		int newDepth = depth - 1;
		int type = -1;
		int theIndex = -1;
		float sphereDist = glm::length(sphereIntersect.first - reflectRay.p0);
		float triangleDist = glm::length(triIntersect.first - reflectRay.p0);
		glm::vec3 theIntersect(0, 0, 0);
		if (sphereIntersect.second != -99 && triIntersect.second == -99) {
			type = 0;
			theIndex = sphereIntersect.second;
			theIntersect = sphereIntersect.first;
		}
		else if (sphereIntersect.second == -99 && triIntersect.second != -99) {
			type = 1;
			theIndex = triIntersect.second;
			theIntersect = triIntersect.first;
		}
		else {
			if (triangleDist > sphereDist) {
				type = 0;
				theIndex = sphereIntersect.second;
				theIntersect = sphereIntersect.first;
			}
			else {
				type = 1;
				theIndex = triIntersect.second;
				theIntersect = triIntersect.first;
			}
		}
		colorlvl += getColor(scene, type, theIndex, newDepth, theIntersect, intersect, reflectRay) * specular;
	}

	return colorlvl;
}

int main(int argc, char* argv[]) {
	FreeImage_Initialise();
	if (argc < 2 || argc > 2) {
		cerr << "argument: scene.test\n";
		exit(-1);
	}

	Scene theScene;
	readfile(argv[1], theScene);
	for (Triangle& tri : theScene.allTriangles) {
		glm::vec4 a1 = glm::vec4(tri.p1.x, tri.p1.y, tri.p1.z, 1.0f);
		glm::vec4 a2 = glm::vec4(tri.p2.x, tri.p2.y, tri.p2.z, 1.0f);
		glm::vec4 a3 = glm::vec4(tri.p3.x, tri.p3.y, tri.p3.z, 1.0f);
		a1 = tri.transform * a1;
		a2 = tri.transform * a2;
		a3 = tri.transform * a3;
		a1 = a1 / a1.w;
		a2 = a2 / a2.w;
		a3 = a3 / a3.w;
		tri.p1 = glm::vec3(a1);
		tri.p2 = glm::vec3(a2);
		tri.p3 = glm::vec3(a3);
	}

	BYTE* pixels = new BYTE[3 * theScene.width * theScene.height];
	int num = 0;
	/*
	for (int x = 0; x < theScene.cameras.size(); x++) {
	*/
	Camera cam = theScene.cameras[0];

	int height = theScene.height;
	int width = theScene.width;
	int i = 0;
	int j = 0;
	omp_set_dynamic(0);     // disable dynamic teams
	omp_set_num_threads(12); // manually set number of threads (cpu dependent)
	#pragma omp parallel for collapse(2) private(i) private(j)
	for (i = 0; i < height; i++) {
		cout << i << endl;
		for (j = 0; j < width; j++) {
			Ray curRay(cam.eye, cam.fovy, theScene.width, theScene.height, (float)i, (float)j, cam.u, cam.v, cam.w);
			pair<glm::vec3, int> sphereIntersect = getInterceptSphere(curRay, theScene);
			pair<glm::vec3, int> triIntersect = getInterceptTriangle(curRay, theScene);
			float sphereDist = glm::length(sphereIntersect.first - curRay.p0);
			float triangleDist = glm::length(triIntersect.first - curRay.p0);
			if (sphereIntersect.second != -99 || triIntersect.second != -99) {
				int type = -1;
				int theIndex = -1;
				glm::vec3 theIntersect(0, 0, 0);
				if (sphereIntersect.second != -99 && triIntersect.second == -99) {
					type = 0;
					theIndex = sphereIntersect.second;
					theIntersect = sphereIntersect.first;
				}
				else if (sphereIntersect.second == -99 && triIntersect.second != -99) {
					type = 1;
					theIndex = triIntersect.second;
					theIntersect = triIntersect.first;
				}
				else {
					if (triangleDist > sphereDist) {
						type = 0;
						theIndex = sphereIntersect.second;
						theIntersect = sphereIntersect.first;
					}
					else {
						type = 1;
						theIndex = triIntersect.second;
						theIntersect = triIntersect.first;
					}
				}
				glm::vec3 theColor = getColor(theScene, type, theIndex, theScene.maxdepth, theIntersect, cam.eye, curRay);
				theColor *= 255.0f;
				pixels[i * width * 3 + j * 3] = (BYTE)theColor.z;
				pixels[i * width * 3 + j * 3 + 1] = (BYTE)theColor.y;
				pixels[i * width * 3 + j * 3 + 2] = (BYTE)theColor.x;
			}
			else {
				pixels[i * width * 3 + j * 3] = (BYTE)0;
				pixels[i * width * 3 + j * 3 + 1] = (BYTE)0;
				pixels[i * width * 3 + j * 3 + 2] = (BYTE)0;
			}
		}
	}
	saveScreenshot(theScene.outputName, theScene.width, theScene.height, pixels);
	FreeImage_DeInitialise;
	delete pixels;
	return 0;

}


