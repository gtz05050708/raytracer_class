// Transform.cpp: implementation of the Transform class.

// Note: when you construct a matrix using mat4() or mat3(), it will be COLUMN-MAJOR
// Keep this in mind in readfile.cpp and display.cpp
// See FAQ for more details or if you're having problems.

#include "Transform.h"

// Helper rotation function.  Please implement this.  
mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
	vec3 newaxis = glm::normalize(axis);
	float angle = glm::radians(degrees);
	glm::mat3 A = glm::mat3(0.0f, newaxis.z, -1 * newaxis.y, -1 * newaxis.z, 0.0f, newaxis.x, newaxis.y, -1 * newaxis.x, 0.0f);
	glm::mat3 aaT = glm::outerProduct(newaxis, newaxis);
	glm::mat3 R = glm::cos(angle) * glm::mat3(1.0f) + (1 - glm::cos(angle)) * aaT + glm::sin(angle) * A;
	return R;
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{
	glm::mat3 R = Transform::rotate(degrees, up);
	eye = R * eye;
	up = R * up;
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
	glm::mat3 R = Transform::rotate(degrees, glm::normalize(glm::cross(eye, up)));
	eye = R * eye;
	up = R * up;
}

mat4 Transform::lookAt(const vec3 &eye, const vec3 &center, const vec3 &up) 
{
	glm::vec3 b = glm::normalize(up);
	glm::vec3 c = glm::normalize(eye);
	glm::vec3 a = glm::cross(b, c);
	glm::mat4 viewMat = glm::mat4(a.x, b.x, c.x, 0, a.y, b.y, c.y, 0, a.z, b.z, c.z, 0, glm::dot((-1.0f * a), eye),
		glm::dot((-1.0f * b), eye), glm::dot((-1.0f * c), eye), 1);
	return viewMat;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{
	float tanFovy = glm::tan(glm::radians(fovy) / 2);
    mat4 ret = mat4(
		1 / (aspect * tanFovy), 0, 0, 0,
		0, 1 / tanFovy, 0, 0,
		0, 0, -1 * (zFar + zNear) / (zFar - zNear), -1,
		0, 0, -1 * 2 * zFar * zNear / (zFar - zNear), 0
	);
    return ret;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
	mat4 ret = mat4(
		sx, 0, 0, 0,
		0, sy, 0, 0,
		0, 0, sz, 0,
		0, 0, 0, 1
	);
    return ret;
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
	mat4 ret = mat4(
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		tx, ty, tz, 1);
    return ret;
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
    vec3 x = glm::cross(up,zvec); 
    vec3 y = glm::cross(zvec,x); 
    vec3 ret = glm::normalize(y); 
    return ret; 
}


Transform::Transform()
{

}

Transform::~Transform()
{

}
