#pragma once
#include "MatrixStruct.h"

class MatrixMath
{
public:

	static float Dot(const Vector3& v1, const Vector3& v2);

	static float Length(const Vector3& v);

	static Vector3 Add(const Vector3& v1, const Vector3& v2);

	static Vector3 Cross(const Vector3& v1, const Vector3& v2);

	static Vector3 Subtract(const Vector3& v1, const Vector3& v2);

	static Vector3 Multiply(const float& v1, const Vector3& v2);

	static Vector3 TransformCoord(Vector3 vector, Matrix4x4 matrix);

	static Vector3 Project(const Vector3& v1, const Vector3& v2);

	static Vector3 ClosestPoint(const Vector3& point, const Segment& segment);

	static Vector3 Normalize(const Vector3& v);

	static Vector3 Perpendicular(const Vector3& vector);

	static Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2);

	static Matrix4x4 MakeTranslateMatrix(const Vector3 translate);

	static Matrix4x4 MakeScaleMatrix(const Vector3 scale);

	static Matrix4x4 MakeRotateXMatrix(float radian);

	static Matrix4x4 MakeRotateYMatrix(float radian);

	static Matrix4x4 MakeRotateZMatrix(float radian);

	static Matrix4x4 MakeAffineMatrix(const Vector3 scale, const Vector3 rotate, const Vector3 translate);

	static Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRetio, float nearClip, float farClip);

	static Matrix4x4 MakeOrthographicMatrix(float left, float top, float right, float bottom, float nearClip, float farClip);

	static Matrix4x4 MakeViewPortMatrix(float left, float top, float width, float height, float minDepth, float maxDepth);

	static Matrix4x4 Inverse(const Matrix4x4& m);

	static Matrix4x4 Transpose(const Matrix4x4& m);

	static Matrix4x4 MakeIdentity4x4();

	static bool IsCollision(const Sphere& s1, const Sphere& s2);

	static bool IsCollision(const Sphere& s1, const Plane& plane);

	static bool IsCollision(const Segment& line, const Plane& plane);

	static bool IsCollision(const Ray& line, const Plane& plane);

	static bool IsCollision(const Line& line, const Plane& plane);

	static bool IsCollision(const Triangle& triangle, const Segment& segment);

	static bool IsCollision(const AABB& aabb1, const AABB& aabb2);

	static bool IsCollision(const AABB& aabb, const Sphere& sphere);
};