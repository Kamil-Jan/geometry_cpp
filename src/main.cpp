#include <iostream>
#include "geometry.h"

int main() {
    Point a(0, 0);
    Point b(0, 1);
    Point c(1, 2);
    Triangle triangle(a, b, c);

    Circle circle = triangle.circumscribedCircle();
    std::cout << "Circumscribed circle's radius = " << circle.radius() << "\n";
    circle = triangle.inscribedCircle();
    std::cout << "Inscribed circle's radius = " << circle.radius() << "\n";
    circle = triangle.ninePointsCircle();
    std::cout << "Nine points circle's radius = " << circle.radius() << "\n";

    Line l = triangle.EulerLine();
    triangle.reflect(l);

    return 0;
}
