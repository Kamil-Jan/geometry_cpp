#include "geometry.h"

Point Point::middle(const Point& p1, const Point& p2) {
    return {(p1.x + p2.x) / 2, (p1.y + p2.y) / 2};
}

Point Point::rotate(double angle) const {
    double cos = std::cos(angle);
    double sin = std::sin(angle);
    return {x * cos - y * sin, x * sin + y * cos};
}

Point Point::rotate(const Point& p, double angle) const {
    return shift(-p).rotate(angle).shift(p);
}

Point Point::shift(const Point& p) const {
    return {x + p.x, y + p.y};
}

Point Point::reflect(const Point& center) const {
    double dx = center.x - x;
    double dy = center.y - y;
    return shift({dx, dy}).shift({dx, dy});
}

Point Point::scale(const Point& center, double coefficient) const {
    double dx = x - center.x;
    double dy = y - center.y;
    dx = coefficient * dx;
    dy = coefficient * dy;
    return Point{dx, dy}.shift(center);
}

Point Point::operator-() const {
    return {x * -1, y * -1};
}

bool operator<(const Point& p1, const Point& p2) {
    return (p1.x < p2.x) || (math_geometry::equal(p1.x, p2.y) && p1.y < p2.y);
}

bool operator==(const Point& p1, const Point& p2) {
    return math_geometry::equal(p1.x, p2.x) && math_geometry::equal(p1.y, p2.y);
}

bool operator!=(const Point& p1, const Point& p2) {
    return !(p1 == p2);
}

Line::Line(double _a, double _b, double _c) : a(_a), b(_b), c(_c) {}

Line::Line(const Point& p1, const Point& p2) {
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = p1.x * p2.y - p2.x * p1.y;
    convertToCanonForm();
}

Line::Line(double k, double shift) {
    a = -k;
    b = 1;
    c = -shift;
}

Line::Line(const Point& p, double k) {
    a = -k;
    b = 1;
    c = k * p.x - p.y;
}

Point Line::intersect(const Line& line) const {
    double x = (b * line.c - line.b * c) / (a * line.b - line.a * b);
    double y = (line.a * c - a * line.c) / (a * line.b - line.a * b);
    return {x, y};
}

Point Line::reflect(const Point& p) const {
    Point proj = project(p);
    Point shi = {proj.x - p.x, proj.y - p.y};
    return p.shift(shi).shift(shi);
}

Point Line::project(const Point& p) const {
    double new_x = (p.x - a * p.y - a * c) / (1 + a * a);
    double new_y = (a * a * p.y - a * p.x - c) / (1 + a * a);
    return {new_x, new_y};
}

Line Line::normal(const Point& p) const {
    Point proj = project(p);
    double n_a = -b;
    double n_b = a;
    double n_c = (b - a) * proj.x - (a + b) * proj.y - c;
    Line n(n_a, n_b, n_c);
    n.convertToCanonForm();
    return n;
}

void Line::convertToCanonForm() {
    if (math_geometry::equal(b, 0)) {
        c /= a;
        a = 1;
    } else {
        a /= b;
        c /= b;
        b = 1;
    }
}

bool operator==(const Line& l1, const Line& l2) {
    return math_geometry::equal(l1.a, l2.a) &&
           math_geometry::equal(l1.b, l2.b) && math_geometry::equal(l1.c, l2.c);
}

bool operator!=(const Line& l1, const Line& l2) {
    return !(l1 == l2);
}

bool operator==(const Shape& shape1, const Shape& shape2) {
    return shape1.equal(shape2);
}

bool operator!=(const Shape& shape1, const Shape& shape2) {
    return !shape1.equal(shape2);
}

math_geometry::NumSign math_geometry::getNumSign(double x) {
    if (x > 0) {
        return math_geometry::NumSign::POSITIVE;
    }
    if (math_geometry::equal(x, 0)) {
        return math_geometry::NumSign::ZERO;
    }
    return math_geometry::NumSign::NEGATIVE;
}

bool math_geometry::equal(double x, double y) {
    return std::abs(x - y) < math_geometry::EPS;
}

double math_geometry::Vec2::length() const {
    return std::sqrt(x * x + y * y);
}

double math_geometry::Vec2::angle(const Vec2& v1, const Vec2& v2) {
    return std::acos(dotProduct(v1, v2) / (v1.length() * v2.length()));
}

double math_geometry::Vec2::dotProduct(const Vec2& v1, const Vec2& v2) {
    return v1.x * v2.x + v1.y * v2.y;
}

double math_geometry::Vec2::zCrossProduct(const Vec2& v1, const Vec2& v2) {
    return v1.x * v2.y - v1.y * v2.x;
}

Ellipse::Ellipse(const Point& p1, const Point& p2, double k) {
    if (p1 < p2) {
        f1 = p1;
        f2 = p2;
    } else {
        f1 = p2;
        f2 = p1;
    }
    a = k / 2;
}

std::pair<Point, Point> Ellipse::focuses() const {
    return {f1, f2};
}

std::pair<Line, Line> Ellipse::directrices() const {
    Line axis(f1, f2);
    return {axis.normal(f1), axis.normal(f2)};
}

double Ellipse::eccentricity() const {
    return getC() / a;
}

Point Ellipse::center() const {
    return Point::middle(f1, f2);
}

double Ellipse::perimeter() const {
    double b = getB();
    return M_PI * (3 * (a + b) - std::sqrt((3 * a + b) * (a + 3 * b)));
}

double Ellipse::area() const {
    return M_PI * a * getB();
}

bool Ellipse::equal(const Shape& another) const {
    if (const Ellipse* elli = dynamic_cast<const Ellipse*>(&another)) {
        return (f1 == elli->f1 && f2 == elli->f2 && a == elli->a);
    }
    return false;
}

bool Ellipse::isCongruentTo(const Shape& another) const {
    using math_geometry::Vec2;
    if (const Ellipse* elli = dynamic_cast<const Ellipse*>(&another)) {
        Ellipse canon1 = getCanonForm();
        Ellipse canon2 = elli->getCanonForm();
        double angle1 =
            Vec2::angle(Vec2({0, 0}, {1, 0}), Vec2(canon1.f1, canon1.f2));
        double angle2 =
            Vec2::angle(Vec2({0, 0}, {1, 0}), Vec2(canon2.f1, canon2.f2));
        canon1.rotate({0, 0}, -angle1);
        canon2.rotate({0, 0}, -angle2);
        return canon1 == canon2;
    }
    return false;
}

bool Ellipse::isSimilarTo(const Shape& another) const {
    using math_geometry::Vec2;
    const Ellipse* elli = dynamic_cast<const Ellipse*>(&another);
    if (elli == nullptr) {
        return false;
    }
    Ellipse canon1 = getCanonForm();
    Ellipse canon2 = elli->getCanonForm();
    double k1 = 1 / Vec2(canon1.f1, canon1.f2).length();
    double k2 = 1 / Vec2(canon2.f1, canon2.f2).length();
    canon1.scale({0, 0}, k1);
    canon2.scale({0, 0}, k2);
    return canon1.isCongruentTo(canon2);
}

bool Ellipse::containsPoint(const Point& p) const {
    Point cen = center();
    double b = getB();
    double dx = (p.x - cen.x);
    double dy = (p.y - cen.y);
    return dx * dx / (a * a) + dy * dy / (b * b) - 1 < math_geometry::EPS;
}

void Ellipse::rotate(const Point& center, double angle) {
    f1 = f1.rotate(center, angle);
    f2 = f2.rotate(center, angle);
}

void Ellipse::reflect(const Point& center) {
    f1 = f1.reflect(center);
    f2 = f2.reflect(center);
}

void Ellipse::reflect(const Line& axis) {
    f1 = axis.reflect(f1);
    f2 = axis.reflect(f2);
}

void Ellipse::scale(const Point& center, double coefficient) {
    f1 = f1.scale(center, coefficient);
    f2 = f2.scale(center, coefficient);
    a *= coefficient;
}

double Ellipse::getB() const {
    double c = getC();
    return std::sqrt(a * a - c * c);
}

double Ellipse::getC() const {
    return math_geometry::Vec2(f1, center()).length();
}

Ellipse Ellipse::getCanonForm() const {
    return Ellipse({0, 0}, f2.shift(-f1), 2 * a);
}

Circle::Circle(const Point& center, double rad)
    : Ellipse(center, center, 2 * rad) {}

double Circle::radius() const {
    return a;
}

Polygon::Polygon(const std::vector<Point>& points) : vertices_(points) {}

size_t Polygon::verticesCount() const {
    return vertices_.size();
}

const std::vector<Point>& Polygon::getVertices() const {
    return vertices_;
}

bool Polygon::isConvex() const {
    using math_geometry::getNumSign;
    using math_geometry::NumSign;
    using math_geometry::Vec2;

    size_t sz = vertices_.size();
    if (sz <= 3) {
        return true;
    }

    NumSign sign = NumSign::ZERO;
    for (size_t i = 0; i < sz; ++i) {
        double z = Vec2::zCrossProduct(
            Vec2(vertices_[i], vertices_[(i + 1) % sz]),
            Vec2(vertices_[(i + 1) % sz], vertices_[(i + 2) % sz]));
        NumSign z_sign = getNumSign(z);
        if (sign == NumSign::ZERO) {
            sign = z_sign;
        } else if (z_sign != sign) {
            return false;
        }
    }
    return true;
}

double Polygon::perimeter() const {
    double p = 0;
    size_t sz = vertices_.size();
    for (size_t i = 0; i < sz; ++i) {
        p +=
            math_geometry::Vec2(vertices_[i], vertices_[(i + 1) % sz]).length();
    }
    return p;
}

double Polygon::area() const {
    double a = 0;
    size_t sz = vertices_.size();
    for (size_t i = 1; i < sz - 1; ++i) {
        a += math_geometry::Vec2::zCrossProduct(
                 math_geometry::Vec2(vertices_[0], vertices_[i]),
                 math_geometry::Vec2(vertices_[0], vertices_[(i + 1) % sz])) /
             2;
    }
    return std::abs(a);
}

bool Polygon::equal(const Shape& another) const {
    if (const Polygon* poly = dynamic_cast<const Polygon*>(&another)) {
        size_t n = vertices_.size();
        size_t m = poly->verticesCount();
        if (n != m) {
            return false;
        }

        const std::vector<Point>& poly_vert = poly->getVertices();
        size_t i = 0;
        size_t j = 0;
        while (j < m && vertices_[i] != poly_vert[j]) {
            ++j;
        }
        if (j == m) {
            return false;
        }
        if (vertices_[i + 1] == poly_vert[(j + 1) % m]) {
            ++i;
            size_t start = j;
            for (size_t k = (start + 1) % m; k != start; k = (k + 1) % m) {
                if (vertices_[i++] != poly_vert[k]) {
                    return false;
                }
            }
            return true;
        }

        ++i;
        size_t start = j;
        for (size_t k = (start - 1 + m) % m; k != start; k = (k - 1 + m) % m) {
            if (vertices_[i++] != poly_vert[k]) {
                return false;
            }
        }
        return true;
    }
    return false;
}

bool Polygon::isCongruentTo(const Shape& another) const {
    if (const Polygon* poly = dynamic_cast<const Polygon*>(&another)) {
        if (vertices_.size() != poly->verticesCount()) {
            return false;
        }
        std::vector<Point> another_vertices = poly->getVertices();
        if (checkCongruency(vertices_, another_vertices)) {
            return true;
        }
        std::reverse(another_vertices.begin(), another_vertices.end());
        return (checkCongruency(vertices_, another_vertices));
    }
    return false;
}

bool Polygon::isSimilarTo(const Shape& another) const {
    if (const Polygon* poly = dynamic_cast<const Polygon*>(&another)) {
        if (vertices_.size() != poly->verticesCount()) {
            return false;
        }
        std::vector<Point> another_vertices = poly->getVertices();
        if (checkSimilarity(vertices_, another_vertices)) {
            return true;
        }
        std::reverse(another_vertices.begin(), another_vertices.end());
        return (checkSimilarity(vertices_, another_vertices));
    }
    return false;
}

bool Polygon::containsPoint(const Point& point) const {
    using math_geometry::getNumSign;
    using math_geometry::NumSign;
    using math_geometry::Vec2;

    NumSign sign = NumSign::ZERO;
    size_t sz = vertices_.size();
    for (size_t i = 0; i < sz; ++i) {
        double z =
            Vec2::zCrossProduct(Vec2(vertices_[i], point),
                                Vec2(vertices_[i], vertices_[(i + 1) % sz]));
        NumSign z_sign = getNumSign(z);
        if (sign == NumSign::ZERO) {
            sign = z_sign;
        } else if (z_sign != NumSign::ZERO && z_sign != sign) {
            return false;
        }
    }
    return true;
}

void Polygon::rotate(const Point& center, double angle) {
    for (Point& p : vertices_) {
        p = p.rotate(center, angle);
    }
}

void Polygon::reflect(const Point& center) {
    for (Point& p : vertices_) {
        p = p.reflect(center);
    }
}

void Polygon::reflect(const Line& axis) {
    for (Point& p : vertices_) {
        p = axis.reflect(p);
    }
}

void Polygon::scale(const Point& center, double coefficient) {
    for (Point& p : vertices_) {
        p = p.scale(center, coefficient);
    }
}

Polygon Polygon::getCanonForm() const {
    std::vector<Point> points = vertices_;
    for (Point& p : points) {
        p = p.shift(-vertices_[0]);
    }
    return Polygon(points);
}

bool Polygon::checkCongruency(const std::vector<Point>& vertices1,
                              const std::vector<Point>& vertices2) {
    using math_geometry::equal;
    using math_geometry::Vec2;

    size_t n = vertices1.size();
    for (size_t i = 0; i < n; ++i) {
        bool is_congruent = true;
        for (size_t j = 0; j < n; ++j) {
            Vec2 side1(vertices1[(i + j) % n], vertices1[(i + j + 1) % n]);
            Vec2 side2(vertices2[j], vertices2[(j + 1) % n]);
            if (!equal(side1.length(), side2.length())) {
                is_congruent = false;
                break;
            }
            double area1 = Vec2::zCrossProduct(
                side1,
                Vec2(vertices1[(i + j) % n], vertices1[(i + j - 1 + n) % n]));
            double area2 = Vec2::zCrossProduct(
                side2, Vec2(vertices2[j], vertices2[(j - 1 + n) % n]));
            if (!equal(std::abs(area1), std::abs(area2))) {
                is_congruent = false;
                break;
            }
        }
        if (is_congruent) {
            return true;
        }
    }
    return false;
}

bool Polygon::checkSimilarity(const std::vector<Point>& vertices1,
                              const std::vector<Point>& vertices2) {
    using math_geometry::equal;
    using math_geometry::Vec2;

    size_t n = vertices1.size();
    for (size_t i = 0; i < n; ++i) {
        bool is_similar = true;
        double k = -1;
        for (size_t j = 0; j < n; ++j) {
            Vec2 side1(vertices1[(i + j) % n], vertices1[(i + j + 1) % n]);
            Vec2 side2(vertices2[j], vertices2[(j + 1) % n]);
            double ratio = side1.length() / side2.length();
            if (k == -1) {
                k = ratio;
            } else if (!equal(k, ratio)) {
                is_similar = false;
                break;
            }
            double area1 = Vec2::zCrossProduct(
                side1,
                Vec2(vertices1[(i + j) % n], vertices1[(i + j - 1 + n) % n]));
            double area2 = Vec2::zCrossProduct(
                side2, Vec2(vertices2[j], vertices2[(j - 1 + n) % n]));
            if (!equal(std::abs(area1), k * k * std::abs(area2))) {
                is_similar = false;
                break;
            }
        }
        if (is_similar) {
            return true;
        }
    }
    return false;
}

Rectangle::Rectangle(const Point& p1, const Point& p3, double ratio) {
    Point center = Point::middle(p1, p3);
    if (ratio < 1) {
        ratio = 1 / ratio;
    }
    Point p2 = p1.rotate(center, M_PI - 2 * std::atan(ratio));
    Point p4 = p2.reflect(center);
    vertices_ = {p1, p2, p3, p4};
}

double Rectangle::perimeter() const {
    using math_geometry::Vec2;
    return 2 * (Vec2(vertices_[0], vertices_[1]).length() +
                Vec2(vertices_[1], vertices_[2]).length());
}

double Rectangle::area() const {
    using math_geometry::Vec2;
    return Vec2(vertices_[0], vertices_[1]).length() *
           Vec2(vertices_[1], vertices_[2]).length();
}

Point Rectangle::center() const {
    return Point::middle(vertices_[0], vertices_[2]);
}

std::pair<Line, Line> Rectangle::diagonals() const {
    return {Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3])};
}

Square::Square(const Point& p1, const Point& p3) {
    Point middle = Point::middle(p1, p3);
    Point shift = p1.shift(-middle);
    Point p2 = {-shift.y + middle.x, shift.x + middle.y};
    Point p4 = {shift.y + middle.x, -shift.x + middle.y};
    vertices_ = {p1, p2, p3, p4};
}

Circle Square::circumscribedCircle() const {
    return Circle(center(),
                  math_geometry::Vec2(vertices_[0], center()).length());
}

Circle Square::inscribedCircle() const {
    return Circle(center(),
                  math_geometry::Vec2(vertices_[0], vertices_[1]).length() / 2);
}

Triangle::Triangle(const Point& p1, const Point& p2, const Point& p3) {
    vertices_ = {p1, p2, p3};
}

Circle Triangle::circumscribedCircle() const {
    Point m1 = Point::middle(vertices_[0], vertices_[1]);
    Point m2 = Point::middle(vertices_[1], vertices_[2]);
    Line bisector1 = Line(vertices_[0], vertices_[1]).normal(m1);
    Line bisector2 = Line(vertices_[1], vertices_[2]).normal(m2);
    Point inter = bisector1.intersect(bisector2);
    return Circle(inter, math_geometry::Vec2(inter, vertices_[0]).length());
}

Circle Triangle::inscribedCircle() const {
    using math_geometry::Vec2;
    Vec2 side1(vertices_[0], vertices_[1]);
    Vec2 side2(vertices_[1], vertices_[2]);
    Vec2 side3(vertices_[2], vertices_[0]);
    double k1 = side1.length() / side2.length();
    double k2 = side2.length() / side3.length();
    Point p1 = {(vertices_[0].x + k1 * vertices_[2].x) / (1 + k1),
                (vertices_[0].y + k1 * vertices_[2].y) / (1 + k1)};
    Point p2 = {(vertices_[1].x + k2 * vertices_[0].x) / (1 + k2),
                (vertices_[1].y + k2 * vertices_[0].y) / (1 + k2)};
    Line median1 = Line(vertices_[1], p1);
    Line median2 = Line(vertices_[2], p2);
    Point inter = median1.intersect(median2);
    Point proj = Line(vertices_[0], vertices_[1]).project(inter);
    return Circle(inter, math_geometry::Vec2(inter, proj).length());
}

Point Triangle::centroid() const {
    double x = (vertices_[0].x + vertices_[1].x + vertices_[2].x) / 3;
    double y = (vertices_[0].y + vertices_[1].y + vertices_[2].y) / 3;
    return {x, y};
}

Point Triangle::orthocenter() const {
    Line side1 = Line(vertices_[0], vertices_[1]);
    Line side2 = Line(vertices_[1], vertices_[2]);
    return side1.normal(vertices_[2]).intersect(side2.normal(vertices_[0]));
}

Line Triangle::EulerLine() const {
    return Line(centroid(), orthocenter());
}

Circle Triangle::ninePointsCircle() const {
    Point m1 = Point::middle(vertices_[0], vertices_[1]);
    Point m2 = Point::middle(vertices_[1], vertices_[2]);
    Point m3 = Point::middle(vertices_[0], vertices_[2]);
    return Triangle(m1, m2, m3).circumscribedCircle();
}
