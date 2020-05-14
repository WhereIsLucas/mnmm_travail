#include "Vector2.h"

class Plan {

private:
    Vector2 position;
public:
    const Vector2 &getPosition() const;

    void setPosition(const Vector2 &position);

    const Vector2 &getNormal() const;

    void setNormal(const Vector2 &normal);

private:
    Vector2 normal;

public:
    Plan();

    ~Plan();

    void initPlan(Vector2, Vector2);

    void initPlanFromCoordinates(Vector2 a, Vector2 b);
};
