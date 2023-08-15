//
// Created by hugo on 22/04/2020.
//

#ifndef MN_BORIS_VECTOR2_H
#define MN_BORIS_VECTOR2_H


class Vector2 {
private:
    double x{};
    double y{};
public:
    double getY();

    double getX();

    double getNorm();

    void setX(double x);

    void setY(double y);

    void setComponents(double xValue, double yValue);

    Vector2 normalize();

    explicit Vector2(double value);

    Vector2(double x, double y);

    Vector2();

    [[maybe_unused]] void display() const;

};

Vector2 projectOntoVector(Vector2 vec, Vector2 vector);

double getDistanceBetweenVectors(Vector2 vector1, Vector2 vector2);

Vector2 crossProduct(Vector2 vector1, Vector2 vector2);

Vector2 operator+(Vector2 vector1, Vector2 vector2);

Vector2 operator-(Vector2 vector1, Vector2 vector2);

Vector2 operator*(double coefficient, Vector2 vector);

#endif //MN_BORIS_VECTOR2_H

