#include <vector>
#include <fstream>
#include <numeric>
#include <cmath>
#include <type_traits>
#include <iostream>
#include <algorithm>
#include <random>

struct DelaunayPoint {
    double x;
    double y;
    bool operator==(const DelaunayPoint& rhs);
};


bool DelaunayPoint::operator==(const DelaunayPoint& rhs) {
    if (fabs(x - rhs.x) <= 100 * std::numeric_limits<double>::epsilon()) {
        if (fabs(y - rhs.y) <= 100 * std::numeric_limits<double>::epsilon()) {
            return true;
        }
    }
    return false;
}


struct Rectangle {
    DelaunayPoint A;
    DelaunayPoint B;
    DelaunayPoint C;
    DelaunayPoint D;
};

double _super_size = 100000;

double _len(double x0, double y0, double x1, double y1) {
    return sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));
}

std::vector<DelaunayPoint> _jarvis(std::vector<DelaunayPoint> S) {
    // https://iq.opengenus.org/gift-wrap-jarvis-march-algorithm-convex-hull/
    std::vector<DelaunayPoint> P;
    if (S.size() < 3) return P;

    auto orientation = [](DelaunayPoint p, DelaunayPoint q, DelaunayPoint r) {
        int d = ((q.y - p.y) * (r.x - q.x)) - ((q.x - p.x) * (r.y - q.y));
        if (d == 0) return 0; // Colinear
        if (d > 0) return 1; // Clockwise
        else return 2; // Counter clockwise
    };

    int l = 0;
    int p = 0;
    int q = 0;
    DelaunayPoint point_on_hull{10000, 10000};
    int i = 0;
    for (auto& s : S) {
        if (s.x < point_on_hull.x) {
            point_on_hull = s;
            l = i;
            p = i;
        }
        ++i;
    }

    bool stop = false;
    while (!stop) {
        P.push_back(S.at(p));
        q = (p + 1) % S.size();
        for (int i = 0; i < S.size(); i++) {
           if (orientation(S.at(p), S.at(i), S.at(q)) == 2)
               q = i;
        }
        p = q;
        if (p == l) stop = true;
    }
    int smallest_y_index = 0;
    double smallest_y = _super_size;
    for (int i = 0; i < P.size(); ++i) {
        if (P.at(i).y < smallest_y) {
            smallest_y = P.at(i).y;
            smallest_y_index = i;
        }
    }
    // We need P on standard form - The vertices are listed in counterclockwise order with smallest y first
    std::rotate(P.begin(), P.begin() + smallest_y_index, P.end()); 

    std::ofstream _jarvis_file;
    _jarvis_file.open("csv/jarvis.csv", std::ios::trunc);
    for (size_t i = 0; i < P.size(); ++i) {
         _jarvis_file << P.at(i).x << ", " << P.at(i).y << "\n";
    }
    _jarvis_file.close();
 
    return P;
}

Rectangle _get_rectangle(double x0, double y0, double x1, double y1, double l1, double l2) {
    auto resized = [](double x0, double y0, double x1, double y1, double d) {
        double len = _len(x0, y0, x1, y1);
        double ux = (x1 - x0) / len;
        double uy = (y1 - y0) / len;
        return DelaunayPoint{x0 + ux * d, y0 + uy * d};
    };
    auto perp_line = [&resized](double x0, double y0, double x1, double y1, double d) {
        double len = _len(x0, y0, x1, y1);
        double ux = (y0 - y1) / len;
        double uy = (x1 - x0) / len;
        d = d + (y1 - y0);
        return DelaunayPoint{x1 + ux * d, y1 + uy * d};
    };

    DelaunayPoint A{x0, y0};
    DelaunayPoint B{resized(x0, y0, x1, y1, l2)};
    x1 = B.x;
    y1 = B.y;
    DelaunayPoint C{perp_line(x0, y0, x1, y1, l1)};
    DelaunayPoint D{perp_line(x1, y1, x0, y0, -l1)};
    
    Rectangle rectangle{A, B, C, D};
    return rectangle;
}

// 0 Point inside, 1 Move up, 2 Move Left, 3 Move Down, 4 Move Right
int _pip(Rectangle P, DelaunayPoint p) { // Input: Convex counterclockwise polygon and a point
    double x0 = P.A.x;
    double y0 = P.A.y;
    double x1 = P.B.x;
    double y1 = P.B.y;
    double x2 = P.C.x;
    double y2 = P.C.y;
    double x3 = P.D.x;
    double y3 = P.D.y;

    double A = -(y1 - y0);
    double B = x1 - x0;
    double C = -(A * x0 + B * y0);
    double D = A * p.x + B * p.y + C;
    if (D < 0) {
        return 1;
    }
    A = -(y2 - y1);
    B = x2 - x1;
    C = -(A * x1 + B * y1);
    D = A * p.x + B * p.y + C;
    if (D < 0) {
        return 2;
    }
    A = -(y3 - y2);
    B = x3 - x2;
    C = -(A * x2 + B * y2);
    D = A * p.x + B * p.y + C;
    if (D < 0) {
        return 3;
    }
    A = -(y0 - y3);
    B = x0 - x3;
    C = -(A * x3 + B * y3);
    D = A * p.x + B * p.y + C;
    if (D < 0) {
        return 4;
    }
    return 0; 
}

void _move_rectangle(Rectangle& R, double distance, int direction) { // 1 Move up, 2 Move Left, 3 Move Down, 4 Move Right
    if (direction == 1) {
        R.A.y += distance;
        R.B.y += distance;
        R.C.y += distance;
        R.D.y += distance;
    } else if (direction == 2) {
        R.A.x -= distance;
        R.B.x -= distance;
        R.C.x -= distance;
        R.D.x -= distance;
    } else if (direction == 3) {
        R.A.y -= distance;
        R.B.y -= distance;
        R.C.y -= distance;
        R.D.y -= distance;
    } else if (direction == 4) {
        R.A.x += distance;
        R.B.x += distance;
        R.C.x += distance;
        R.D.x += distance;
    }
}

void rotating_calipers(std::vector<DelaunayPoint> P) {
    //https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1381&context=cstech
    auto perp_dist = [](double x1, double y1, double x2, double y2, double m) {
        if (m == 0) {
            return (fabs(x1 - x2));
        } else {
            double b1 = y1 - m * x1;
            double b2 = y2 - m * x2;
            return fabs(b2 - b1) / sqrt(pow(m, 2) + 1);
        }
    };

    auto dot_pr = [](int cur, int nxt, const std::vector<DelaunayPoint>& polygon) {
        // Returns the dot product of the vectors from vertices cur to cur + 1 and nxt to nxt + 1
        int cur1 = cur + 1;
        if (cur1 == polygon.size()) cur1 = 0;
        int nxt1 = nxt + 1;
        if (nxt1 == polygon.size()) nxt1 = 0;
        return ((polygon.at(cur1).x - polygon.at(cur).x) * (polygon.at(nxt1).x - polygon.at(nxt).x)) +
            ((polygon.at(cur1).y - polygon.at(cur).y) * (polygon.at(nxt1).y - polygon.at(nxt).y));
    };

    auto cross_pr = [](int cur, int nxt, const std::vector<DelaunayPoint>& polygon) {
        // Returns the cross product of the vectors from vertices cur to cur + 1 and nxt to nxt + 1
        int cur1 = cur + 1;
        if (cur1 == polygon.size()) cur1 = 0;
        int nxt1 = nxt + 1;
        if (nxt1 == polygon.size()) nxt1 = 0;
        return ((polygon.at(cur1).x - polygon.at(cur).x) * (polygon.at(nxt1).y - polygon.at(nxt).y)) -
            ((polygon.at(cur1).y - polygon.at(cur).y) * (polygon.at(nxt1).x - polygon.at(nxt).x));
    };

    auto inc = [](int index, int n) {
        ++index;
        return index % n;
    };

    double area = _super_size * _super_size;
    double l1 = 0;
    double l2 = 0;
    int vertex = 0;

    int j = 1;
    int n = P.size();
    for (int i = 1; i < n; ++i) {
        while (dot_pr(i, j, P) > 0.0) {
            j = inc(j, n);
        }
        int k = 0;
        if (i == 0) k = j;
        while (cross_pr(i, k, P) > 0.0) {
            k = inc(k, n);
        }
        int m = 0;
        if (i == 0) m = k;
        while (dot_pr(i, m, P) < 0.0) {
            m = inc(m, n);
        }

        double d1 = 0;
        double d2 = 0;
        if (P.at(i).x - P.at(inc(i, n)).x <= 100 * std::numeric_limits<double>::epsilon()) { // Does not seem to work
            d1 = fabs(P.at(k).x - P.at(i).x);
            d2 = fabs(P.at(m).y - P.at(j).y);
        } else if (P.at(i).y - P.at(inc(i, n)).y <= 100 * std::numeric_limits<double>::epsilon()) {
            d1 = fabs(P.at(k).y - P.at(i).y);
            d2 = fabs(P.at(m).x - P.at(j).x);
        } else {
            double slope = (P.at(inc(i, n)).y - P.at(i).y) / (P.at(inc(i, n)).x - P.at(i).x);
            d1 = perp_dist(P.at(i).x, P.at(i).y, P.at(k).x, P.at(k).y, slope);
            d2 = perp_dist(P.at(j).x, P.at(j).y, P.at(m).x, P.at(m).y, -1 / slope);
        }
        if (d1 * d2 < area && d1 != 0 && d2 != 0) {
            vertex = i;
            l1 = d1;
            l2 = d2;
            area = l1 * l2;
        }
    }
    
    double x0 = P.at(vertex).x;
    double y0 = P.at(vertex).y;
    double x1 = P.at(inc(vertex, n)).x;
    double y1 = P.at(inc(vertex, n)).y;
    Rectangle rectangle = std::move(_get_rectangle(x0, y0, x1, y1, l1, l2));

    // double distance = 0.2;
    // bool fail = false;
    // for (auto& p : P) {
    //     int direction = _pip(rectangle, p);
    //     int counter = 0;
    //     while (direction != 0 && !fail) {
    //         _move_rectangle(rectangle, distance, direction);
    //         int prev_direction = direction;
    //         direction = _pip(rectangle, p);
    //         if (direction == 1 && prev_direction == 3) direction = 0;
    //         if (direction == 2 && prev_direction == 4) direction = 0;
    //         if (counter > 1000) {
    //             fail = true;
    //             std::cout << "Fail. Some points probably have the same x coordiante\n";
    //         }
    //         ++counter;
    //     }
    // }

    std::ofstream _rotating_calipers_file;
    _rotating_calipers_file.open("csv/rotating_calipers.csv", std::ios::trunc);
    _rotating_calipers_file << rectangle.A.x << ", " << rectangle.A.y << "\n";
    _rotating_calipers_file << rectangle.B.x << ", " << rectangle.B.y << "\n";
    _rotating_calipers_file << rectangle.C.x << ", " << rectangle.C.y << "\n";
    _rotating_calipers_file << rectangle.D.x << ", " << rectangle.D.y << "\n";
    _rotating_calipers_file.close();

}

int main() {
    double lower_bound = -10;
    double upper_bound = 10;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    std::vector<DelaunayPoint> P;
    for (size_t i = 0; i < 5; ++i) {
        P.push_back(DelaunayPoint{unif(re), unif(re)});
    }
    rotating_calipers(_jarvis(P));
}
