#pragma once

#include <algorithm>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <vector>

struct Point {
    double x, y;

    Point() : x(0), y(0) {}
    Point(double _x, double _y) : x(_x), y(_y) {}

    Point rotate(double angle) const;
    Point rotate(const Point& p, double angle) const;
    Point shift(const Point& p) const;
    Point reflect(const Point& center) const;
    Point scale(const Point& center, double coefficient) const;
    Point operator-() const;
    static Point middle(const Point& p1, const Point& p2);
};

bool operator<(const Point& p1, const Point& p2);
bool operator==(const Point& p1, const Point& p2);
bool operator!=(const Point& p1, const Point& p2);

class Line {
  public:
    double a, b, c;
    Line(double _a, double _b, double _c);
    Line(const Point& p1, const Point& p2);
    Line(double k, double shift);
    Line(const Point& p, double k);

    Point intersect(const Line& line) const;
    Point reflect(const Point& p) const;
    Point project(const Point& p) const;
    Line normal(const Point& p) const;

  private:
    void convertToCanonForm();
};

bool operator==(const Line& l1, const Line& l2);
bool operator!=(const Line& l1, const Line& l2);

class Shape {
  public:
    virtual ~Shape() = default;
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool equal(const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;
    virtual void rotate(const Point& center, double angle) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis) = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
};

bool operator==(const Shape& shape1, const Shape& shape2);
bool operator!=(const Shape& shape1, const Shape& shape2);

namespace math_geometry {
const double EPS = 1e-5;
bool equal(double x, double y);

enum class NumSign {
    ZERO,
    POSITIVE,
    NEGATIVE,
};

NumSign getNumSign(double x);

class Vec2 {
  public:
    double x, y;

    Vec2(const Point& p1, const Point& p2) {
        x = p2.x - p1.x;
        y = p2.y - p1.y;
    }

    double length() const;
    static double angle(const Vec2& v1, const Vec2& v2);
    static double dotProduct(const Vec2& v1, const Vec2& v2);
    static double zCrossProduct(const Vec2& v1, const Vec2& v2);
};
};  // namespace math_geometry

class Ellipse : public Shape {
  public:
    Ellipse(const Point& p1, const Point& p2, double k);

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    double perimeter() const final;
    double area() const final;
    bool equal(const Shape& another) const final;
    bool isCongruentTo(const Shape& another) const final;
    bool isSimilarTo(const Shape& another) const final;
    bool containsPoint(const Point& p) const final;
    void rotate(const Point& center, double angle) final;
    void reflect(const Point& center) final;
    void reflect(const Line& axis) final;
    void scale(const Point& center, double coefficient) final;

  protected:
    Point f1, f2;
    double a;

    double getB() const;
    double getC() const;
    Ellipse getCanonForm() const;
};

class Circle : public Ellipse {
  public:
    Circle(const Point& center, double rad);
    double radius() const;
};

class Polygon : public Shape {
  public:
    Polygon(){};
    Polygon(const std::vector<Point>& points);

    template <typename... Args>
    Polygon(Args... args) {
        (vertices_.push_back(args), ...);
    }

    size_t verticesCount() const;
    const std::vector<Point>& getVertices() const;
    bool isConvex() const;

    double perimeter() const override;
    double area() const override;
    bool equal(const Shape& another) const final;
    bool isCongruentTo(const Shape& another) const final;
    bool isSimilarTo(const Shape& another) const final;
    bool containsPoint(const Point& point) const final;
    void rotate(const Point& center, double angle) final;
    void reflect(const Point& center) final;
    void reflect(const Line& axis) final;
    void scale(const Point& center, double coefficient) final;

  protected:
    std::vector<Point> vertices_;
    Polygon getCanonForm() const;
    static bool checkCongruency(const std::vector<Point>& vertices1,
                                const std::vector<Point>& vertices2);
    static bool checkSimilarity(const std::vector<Point>& vertices1,
                                const std::vector<Point>& vertices2);
};

class Rectangle : public Polygon {
  public:
    Rectangle(){};
    Rectangle(const Point& p1, const Point& p3, double ratio);
    double perimeter() const final;
    double area() const final;
    Point center() const;
    std::pair<Line, Line> diagonals() const;
};

class Square : public Rectangle {
  public:
    Square(const Point& p1, const Point& p3);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
};

class Triangle : public Polygon {
  public:
    Triangle(const Point& p1, const Point& p2, const Point& p3);
    Circle circumscribedCircle() const;
    Circle inscribedCircle() const;
    Point centroid() const;
    Point orthocenter() const;
    Line EulerLine() const;
    Circle ninePointsCircle() const;
};
