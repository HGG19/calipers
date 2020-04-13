#include <vector>
#include <fstream>
#include <numeric>
#include <cmath>
#include <type_traits>
#include <iostream>

struct DelaunayPoint {
    double x;
    double y;
    bool operator==(const DelaunayPoint& rhs);
};

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

bool DelaunayPoint::operator==(const DelaunayPoint& rhs) {
    if (fabs(x - rhs.x) <= 100 * std::numeric_limits<double>::epsilon()) {
        if (fabs(y - rhs.y) <= 100 * std::numeric_limits<double>::epsilon()) {
            return true;
        }
    }
    return false;
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

void rotating_calipers() {
    std::ofstream _jarvis_file;
    std::ofstream _rotating_calipers_file;
    
    std::vector<DelaunayPoint> P;
    P.push_back(DelaunayPoint{-1, 0});
    P.push_back(DelaunayPoint{4, 0});
    P.push_back(DelaunayPoint{3, 3});
    P.push_back(DelaunayPoint{1, 3});
    P.push_back(DelaunayPoint{-1, 0});

    //https://docs.lib.purdue.edu/cgi/viewcontent.cgi?article=1381&context=cstech
    auto perp_dist = [](double x1, double y1, double x2, double y2, double m) {
        if (m == 0) {
            std::cout << "M == 0\n";
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
    double s = 0;
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
        double slope = 0;
        if (P.at(i).x - P.at(inc(i, n)).x <= 100 * std::numeric_limits<double>::epsilon()) {
            d1 = fabs(P.at(k).x - P.at(i).x);
            d2 = fabs(P.at(m).y - P.at(j).y);
            slope = std::numeric_limits<double>::infinity();
        } else if (P.at(i).y - P.at(inc(i, n)).y <= 100 * std::numeric_limits<double>::epsilon()) {
            d1 = fabs(P.at(k).y - P.at(i).y);
            d2 = fabs(P.at(m).x - P.at(j).x);
            slope = 0;
        } else {
            slope = (P.at(inc(i, n)).y - P.at(i).y) / (P.at(inc(i, n)).x - P.at(i).x);
            d1 = perp_dist(P.at(i).x, P.at(i).y, P.at(k).x, P.at(k).y, slope);
            d2 = perp_dist(P.at(j).x, P.at(j).y, P.at(m).x, P.at(m).y, -1 / slope);
        }
        if (d1 * d2 < area && d1 != 0 && d2 != 0) {
            vertex = i;
            l1 = d1;
            l2 = d2;
            area = l1 * l2;
            s = slope;
        }
    }
    
    double x0 = P.at(vertex).x;
    double y0 = P.at(vertex).y;
    double x1 = P.at(inc(vertex, n)).x;
    double y1 = P.at(inc(vertex, n)).y;
    Rectangle rectangle = std::move(_get_rectangle(x0, y0, x1, y1, l1, l2));

    _rotating_calipers_file.open("rotating_calipers.csv", std::ios::trunc);
    _rotating_calipers_file << rectangle.A.x << ", " << rectangle.A.y << "\n";
    _rotating_calipers_file << rectangle.B.x << ", " << rectangle.B.y << "\n";
    _rotating_calipers_file << rectangle.C.x << ", " << rectangle.C.y << "\n";
    _rotating_calipers_file << rectangle.D.x << ", " << rectangle.D.y << "\n";
    _rotating_calipers_file.close();

    _jarvis_file.open("jarvis.csv", std::ios::trunc);
    for (int i = 0; i < P.size(); ++i) {
         _jarvis_file << P.at(i).x << ", " << P.at(i).y << "\n";
    }
    _jarvis_file.close();

}

int main() {
    rotating_calipers();
}
