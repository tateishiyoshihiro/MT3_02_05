#pragma once

struct Vector3 final {
	float x;
	float y;
	float z;
};

struct Matrix4x4 final {
	float m[4][4];
};

struct Sphere final {
	Vector3 center;
	float radius;
};

struct Line final {
	Vector3 origin;
	Vector3 diff;
};

struct Ray final {
	Vector3 origin;
	Vector3 diff;
};

struct Segment final {
	Vector3 origin;
	Vector3 diff;
};

struct Plane final {
	Vector3 normal;
	float distance;
};

struct Triangle {
	Vector3 vertices[3];
};

struct AABB {
	Vector3 min;
	Vector3 max;
};