// Small numerical helpers: unit parsing, quadrature nodes, special functions.
#include "stacktherm/stacktherm.hpp"

#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace st {

namespace {

// lower-case, drop separators, map the micro sign to 'u'
std::string unit_key(const std::string& unit) {
    std::string out;
    for (size_t i = 0; i < unit.size(); ++i) {
        unsigned char c = static_cast<unsigned char>(unit[i]);
        if (c == 0xC2 && i + 1 < unit.size()) {            // U+00B5 micro, U+00B7 middle dot
            unsigned char d = static_cast<unsigned char>(unit[++i]);
            if (d == 0xB5) out += 'u';
            continue;
        }
        if (c == 0xCE && i + 1 < unit.size()) {            // U+03BC greek mu
            if (static_cast<unsigned char>(unit[++i]) == 0xBC) out += 'u';
            continue;
        }
        if (c == ' ' || c == '^' || c == '*' || c == '.') continue;
        out += static_cast<char>(std::tolower(c));
    }
    return out;
}

const std::map<std::string, double>& scales() {
    static const std::map<std::string, double> s = {
        {"m", 1.0}, {"mm", 1e-3}, {"um", 1e-6}, {"nm", 1e-9}, {"a", 1e-10},
        {"w/m2", 1.0}, {"w/cm2", 1e4}, {"w/mm2", 1e6}, {"kw/cm2", 1e7},
        {"w/m2k", 1.0}, {"mw/m2k", 1e6}, {"gw/m2k", 1e9}, {"kw/m2k", 1e3},
        {"m2k/w", 1.0}, {"mm2k/w", 1e-6}, {"cm2k/w", 1e-4}, {"m2k/gw", 1e-9}, {"m2k/mw", 1e-6},
        {"k", 1.0}};
    return s;
}

}  // namespace

double parse_quantity(const std::string& text, const std::string& default_unit) {
    const char* begin = text.c_str();
    char* end = nullptr;
    double number = std::strtod(begin, &end);
    if (end == begin) throw std::invalid_argument("cannot parse quantity '" + text + "'");
    std::string unit(end);
    unit.erase(0, unit.find_first_not_of(" \t"));
    unit.erase(unit.find_last_not_of(" \t") + 1);
    if (unit.empty()) {
        if (default_unit.empty()) return number;
        return number * scales().at(unit_key(default_unit));
    }
    if (unit == "C" || unit == "degC" || unit == "\xC2\xB0""C") return number + 273.15;
    auto it = scales().find(unit_key(unit));
    if (it == scales().end())
        throw std::invalid_argument("unknown unit '" + unit + "' in '" + text + "'");
    return number * it->second;
}

double parse_interface(const std::string& text) {
    std::string key = unit_key(text);
    const std::string suffix = "w/m2k";
    if (key.size() >= suffix.size() && key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0)
        return 1.0 / parse_quantity(text);
    return parse_quantity(text);
}

void gauss_legendre(int n, double a, double b, Vec& x, Vec& w) {
    x.assign(n, 0.0);
    w.assign(n, 0.0);
    for (int i = 0; i < (n + 1) / 2; ++i) {
        double z = std::cos(kPi * (i + 0.75) / (n + 0.5)), pp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p1 = 1.0, p2 = 0.0;
            for (int j = 1; j <= n; ++j) {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
            }
            pp = n * (z * p1 - p2) / (z * z - 1.0);
            double dz = p1 / pp;
            z -= dz;
            if (std::abs(dz) < 1e-15) break;
        }
        double xm = 0.5 * (b + a), xl = 0.5 * (b - a);
        x[i] = xm - xl * z;
        x[n - 1 - i] = xm + xl * z;
        w[i] = w[n - 1 - i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    }
}

double expint(int n, double x) {
    const double euler = 0.5772156649015329, eps = 1e-15, tiny = 1e-300;
    if (x > 700.0) return 0.0;
    if (x == 0.0) return n > 1 ? 1.0 / (n - 1) : kInf;
    if (x > 1.0) {                                   // Lentz continued fraction
        double b = x + n, c = 1.0 / tiny, d = 1.0 / b, h = d;
        for (int i = 1; i <= 10000; ++i) {
            double an = -double(i) * (n - 1 + i);
            b += 2.0;
            d = 1.0 / (an * d + b);
            c = b + an / c;
            double del = c * d;
            h *= del;
            if (std::abs(del - 1.0) < eps) break;
        }
        return h * std::exp(-x);
    }
    double ans = (n - 1) != 0 ? 1.0 / (n - 1) : -std::log(x) - euler, fact = 1.0;
    for (int i = 1; i <= 10000; ++i) {
        fact *= -x / i;
        double del;
        if (i != n - 1) {
            del = -fact / (i - (n - 1));
        } else {
            double psi = -euler;
            for (int k = 1; k <= n - 1; ++k) psi += 1.0 / k;
            del = fact * (-std::log(x) + psi);
        }
        ans += del;
        if (std::abs(del) < std::abs(ans) * eps) break;
    }
    return ans;
}

namespace {

// Faddeeva function w(z) for Im z >= 0: Weideman's rational approximation
// (SIAM J. Numer. Anal. 31, 1497 (1994)), N = 32.
std::complex<double> faddeeva(std::complex<double> z) {
    static const int N = 32;
    static const double L = std::sqrt(N / std::sqrt(2.0));
    static const std::vector<double> a = [] {
        const int M = 2 * N, M2 = 2 * M;
        std::vector<double> f(M2, 0.0);
        for (int k = -M + 1; k <= M - 1; ++k) {
            double t = L * std::tan(k * kPi / M / 2.0);
            f[k + M] = std::exp(-t * t) * (L * L + t * t);
        }
        std::vector<double> coef(N);
        for (int j = 1; j <= N; ++j) {               // real part of the FFT of fftshift(f)
            double re = 0.0;
            for (int m = 0; m < M2; ++m) re += f[(m + M) % M2] * std::cos(2.0 * kPi * j * m / M2);
            coef[N - j] = re / M2;                   // flipped: highest degree first
        }
        return coef;
    }();
    const std::complex<double> I(0.0, 1.0);
    std::complex<double> lz = L - I * z, Z = (L + I * z) / lz, p = 0.0;
    for (double c : a) p = p * Z + c;
    return 2.0 * p / (lz * lz) + (1.0 / std::sqrt(kPi)) / lz;
}

}  // namespace

double voigt(double x, double sigma, double gamma) {
    if (sigma <= 0.0) return gamma / kPi / (x * x + gamma * gamma);
    if (gamma <= 0.0) return std::exp(-0.5 * x * x / (sigma * sigma)) / (sigma * std::sqrt(2.0 * kPi));
    std::complex<double> z(x / (sigma * std::sqrt(2.0)), gamma / (sigma * std::sqrt(2.0)));
    return faddeeva(z).real() / (sigma * std::sqrt(2.0 * kPi));
}

std::string data_dir() {
    if (const char* env = std::getenv("STACKTHERM_DATA")) return env;
#ifdef STACKTHERM_DATA_DIR
    return STACKTHERM_DATA_DIR;
#else
    return "data";
#endif
}

std::string json_escape(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            default: out += c;
        }
    }
    return out;
}

std::string fmt(double v, int significant) {
    if (!std::isfinite(v)) return "null";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*g", significant, v);
    return buf;
}

}  // namespace st
