#include <iostream>
#include <string>
#include <math.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <stack>
#include <limits>
#include <queue>
#include <list>
#include <cstring>
#include <bitset>

#define PI 3.14159265

using namespace std;

// Brand new comp geom template!

const long double EPS = 1e-6;

struct point {
    long double x, y;
    
    point() {}
    point(long double a, long double b) : x(a), y(b) {}
    
    point operator-(const point & o) {
        return point(x - o.x, y - o.y);
    }
    point operator+(const point & o) {
        return point(x + o.x, y + o.y);
    }
    bool operator==(const point & o) {
        return equals(x, o.x) && equals(y, o.y);
    }
    bool operator!=(const point & o) {
        return !equals(x, o.x) || !equals(y, o.y);
    }
    bool equals(long double a, long double b) {
        return fabs(a - b) < EPS;
    }
    point operator*(const long double m) {
        return point(m * x, m * y);
    }
    point operator/(const long double m) {
        return point(x / m, y / m);
    }
    bool operator<(const point & o) const {
        if (x == o.x) return y > o.y;
        else return x < o.x;
    }
    long double dist(const point & o) const {
        return sqrt((x - o.x) * (x - o.x) + (y - o.y) * (y - o.y));
    }
    long double dot(const point & o) const {
        return x * o.x + y * o.y;
    }
    long double cross(const point & o) const {
        return x * o.y - y * o.x;
    }
    void print(int p) const {
        cout << fixed << setprecision(p) << x << " " << y << endl;
    }
    
    
    // Rotates point b around point by alpha radians counterclockwise
    point rotate(const point & b, long double alpha) {
        point a = point(x, y);
        point x = b;
        point c = x - a;
        
        long double newx = c.x * cos(alpha) - c.y * sin(alpha);
        long double newy = c.x * sin(alpha) + c.y * cos(alpha);
        
        return point(newx, newy) + a;
        
    }
    
};

struct line {
    
    long double a, b, c;
    
    point p1, p2;
    
    line() {}
    
    line(point A, point B) {
        a = B.y - A.y;
        b = A.x - B.x;
        c = a * A.x + b * A.y;
        p1.x = A.x;
        p1.y = A.y;
        p2.x = B.x;
        p2.y = B.y;
    }
    
    line perpendicular(point p) {
        if (a == 0) {
            return line(p, point(p.x, p.y + 1));
        }
        long double d = -b * p.x + a * p.y;
        long double newx = p.x + 1;
        long double newy = (d + newx * b) / a;
        point newpoint = point(newx, newy);
        return line(p, newpoint);
    }
    
    line perpendicularBisector() {
        point midpoint = (p1 + p2) / 2;
        return perpendicular(midpoint);
    }
    
    point intersect(const line & y) {
        long double det = a * y.b - y.a * b;
        
        if (det == 0) {
            // PARALLEL
            // GIVE BETTER INPUTS PLS
            return point(-1e9, -1e9);
            
        } else {
            long double xx = (y.b * c - b * y.c) / det;
            long double yy = (a * y.c - y.a * c) / det;
            return point(xx, yy);
        }
    }
    
    bool equals(long double a, long double b) {
        return fabs(a - b) < EPS;
    }
    
    bool intersectsSegment(const line & y) {
        point p = intersect(y);
        // if (p == y.p1 || p == y.p2 || p == p1 || p == p2) return false;
        bool b1 = equals(p.dist(y.p1) + p.dist(y.p2), y.p1.dist(y.p2));
        bool b2 = equals(p.dist(p1) + p.dist(p2), p1.dist(p2));
        return b1 && b2;
    }
    
    point shift(const point & p, const long double & d, long double multi) {
        
        if (a == 0) {
            point np = point(0, multi * d) + p;
            return (np);
        }
        
        point x = point(b * sqrt(d * d / (a * a + b * b)), a * sqrt(d * d / (a * a + b * b)));
        return (x * multi) + p;
        
    }
    
    point reflect(point p) {
        line l = perpendicular(p);
        point y = intersect(l);
        return (y - (p - y));
    }
    
    // Distance from point C to line ab (or segment if isSegment is true)
    long double pointDist(point C, bool isSegment) {
        point A = p1;
        point B = p2;
        long double dist = (B - A).cross(C - A) / sqrt((B - A).dot(B - A));
        
        if (isSegment) {
            // Check the two endpoints
            long double dot1 = (C - B).dot(B - A);
            if (dot1 > 0) return sqrt((B - C).dot(B - C));
            long double dot2 = (C - A).dot(A - B);
            if (dot2 > 0) return sqrt((A - C).dot(A - C));
        }
        
        return fabs(dist);
    }
    
    bool isLeft(point a) {
        // equal to 0 means a is on the line
        // < 0 means a is to the right of the line
        return ((p2.x - p1.x) * (a.y - p1.y) - (p2.y - p1.y) * (a.x - p1.x)) > 0;
    }
    
    bool isParallel(line y) {
        long double det = a * y.b - y.a * b;
        return det == 0;
    }
    
    void print() {
        cout << p1.x << " " << p1.y << " " << p2.x << " " << p2.y << endl;
    }
    
};

struct circle {
    point c;
    long double r;
    
    circle() {}
    circle(point p, long double a) : c(p), r(a) {}
    
    circle(point a, point b, point p) {
        line ab = line(a, b);
        line bp = line(b, p);
        line d = ab.perpendicularBisector();
        line e = bp.perpendicularBisector();
        c = d.intersect(e);
        r = c.dist(a);
    }
    
};

long double calcArea(vector<point> p) {
    long double ret = 0;
    for (int i = 1; i < p.size() - 1; i++) {
        point np1 = p[i] - p[0];
        point np2 = p[i + 1] - p[0];
        ret += np1.cross(np2);
    }
    return fabs(ret / 2.0);
}



const long double INF = 1e18;

vector<point> convexHull(vector<point> x, bool onEdge) {
    
    vector<bool> used(x.size(), false);
    int p = 0;
    for (int i = 1; i < x.size(); i++) {
        if (x[i] < x[p]) p = i;
    }
    
    int start = p;
    vector<point> ret;
    do {
        int n = -1;
        long double dist = onEdge ? INF : 0;
        for (int i = 0; i < x.size(); i++) {
            if (i == p) continue;
            if (used[i]) continue;
            if (n == -1) n = i;
            
            long double cross = (x[i] - x[p]).cross(x[n] - x[p]);
            
            long double d = (x[i] - x[p]).dot(x[n] - x[p]);
            
            if (cross < 0) {
                n = i;
                dist = d;
            } else if (cross == 0) {
                if (onEdge && d < dist) {
                    dist = d;
                    n = i;
                } else if (!onEdge && d > dist) {
                    dist = d;
                    n = i;
                }
            }
        }
        p = n;
        used[p] = true;
        ret.push_back(x[p]);
    } while (start != p);
    
    return ret;
}

struct rect {
    point a, b, c, d;
    
    rect(point x, point y, point z) {
        if (x.dist(y) > x.dist(z) && x.dist(y) > y.dist(z)) {
            a = x;
            b = z;
            c = y;
        } else if (x.dist(z) > x.dist(y) && x.dist(z) > y.dist(z)) {
            a = x;
            b = y;
            c = z;
        } else {
            a = y;
            b = x;
            c = z;
        }
        point m = (a + c) / 2;
        d = b + (m - b) * 2;
    }
    bool equals(long double a, long double b) {
        return fabs(a - b) < EPS;
    }
    
};

bool isPointInPolygon(vector<point> v, point p) {
    int i, j, c = 0;
    for (i = 0, j = v.size() - 1; i < v.size(); j = i++) {
        if (((v[i].y + EPS > p.y) != (v[j].y + EPS > p.y)) && (p.x < (v[j].x - v[i].x) * (p.y - v[i].y) / (v[j].y - v[i].y) + v[i].x + EPS)) c = !c;
    }
    return c;
}

// FIGURE OUT WHY THIS DOESN'T WORK LATER
/*
bool isPointInPolygon(vector<point> v, point p) {
    const point outside = point(-5, p.y);
    line l = line(outside, p);

    vector<line> e;
    for (int i = 0; i < v.size(); i++) {
        e.push_back(line(v[i], v[(i + 1) % v.size()]));
    }
    int cnt = 0;
    for (int i = 0; i < v.size(); i++) {
        point intersect = l.intersect(e[i]);
        point maxi = point(max(e[i].p1.x, e[i].p2.x), max(e[i].p1.y, e[i].p2.y));
        point mini = point(min(e[i].p1.x, e[i].p2.x), min(e[i].p1.y, e[i].p2.y));
        if (intersect.x > mini.x && intersect.x < maxi.x && intersect.y > mini.y && intersect.y < maxi.y) {
            if (intersect.x < p.x) cnt++;
        }
    }
    return cnt % 2 == 1;
}*/

struct loc {
    point p;
    long double d;
    bool operator<(const loc & o) const {
        return o.d < d;
    }
    
    loc(point a) : p(a), d(0) {}
    loc(point a, long double dist) : p(a), d(dist) {}
};

struct point3 {
    long double x, y, z;
    point3() {}
    point3(long double a, long double b, long double c) : x(a), y(b), z(c) {}
    
    point3 operator-(const point3 & o) const {
        return point3(x - o.x, y - o.y, z - o.z);
    }
    point3 operator+(const point3 & o) const {
        return point3(x + o.x, y + o.y, z + o.z);
    }
    long double dot(const point3 & o) const {
        return x * o.x + y * o.y + z * o.z;
    }
    bool operator==(const point3 & o) {
        return equals(x, o.x) && equals(y, o.y);
    }
    bool operator!=(const point3 & o) {
        return !equals(x, o.x) || !equals(y, o.y) || !equals(z, o.z);
    }
    bool equals(long double a, long double b) {
        return fabs(a - b) < EPS;
    }
    point3 operator*(const long double m) {
        return point3(m * x, m * y, m * z);
    }
    point3 operator/(const long double m) {
        return point3(x / m, y / m, z / m);
    }
    
};

struct line3 {
    point3 s, e;
    
    
    
    line3() {}
    line3(point3 a, point3 b) : s(a), e(b) {}
    
    
    // test if r >= 0 for ray, and 0 <= r <= 1 for segment
    point3 param(long double r) {
        return s + (e - s) * r;
    }
};

bool equals(long double a, long double b) {
    return fabs(a - b) < EPS;
}

int main() {
    
    
    
}
