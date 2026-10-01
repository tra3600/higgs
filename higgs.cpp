// higgs.cpp — Mécanisme de Higgs et brisure spontanée de symétrie.
//
// Modèle : champ scalaire à N composantes réelles (N=1 : Z2, N=2 : U(1) "chapeau
// mexicain", N=4 : doublet du Modèle Standard) avec le potentiel
//
//     V(phi) = -1/2 mu^2 |phi|^2 + 1/4 lambda |phi|^4 + mu^4/(4 lambda)
//
// (la constante est choisie pour que V = 0 dans le vide). Résultats analytiques :
//     vev            v      = mu / sqrt(lambda)
//     masse du Higgs m_h    = sqrt(2) mu = sqrt(2 lambda) v
//     Goldstones     N-1 modes de masse nulle
//     transition     T_c^2  = 24 mu^2 / ((N+2) lambda)    (correction thermique)
//     mur de domaine phi(x) = v tanh(mu x / sqrt 2),  énergie (2 sqrt2 / 3) mu^3/lambda
//
// Compilation :  g++ -O2 -std=c++17 higgs.cpp -o higgs      (ou : make)
// Aide        :  ./higgs --help

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using Vec = std::vector<double>;

// ───────────────────────────── Modèle ─────────────────────────────

struct Params {
    double mu = 1.0;       // paramètre de masse (dimension d'une énergie)
    double lambda = 0.5;   // auto-couplage quartique
    int N = 1;             // nombre de composantes réelles du champ
    int points = 61;       // nombre de points des tables / courbes
    double range = 2.0;    // le champ est balayé sur [-range*v, +range*v]
    std::string csv;       // fichier CSV de sortie (optionnel)

    double v() const { return mu / std::sqrt(lambda); }
    double mh() const { return std::sqrt(2.0) * mu; }
    double Tc() const { return std::sqrt(24.0 * mu * mu / ((N + 2) * lambda)); }
    double thermal_c() const { return (N + 2) * lambda / 24.0; }  // m^2(T) = -mu^2 + c T^2
};

// V en fonction de la norme r = |phi|
double V_radial(const Params& p, double r) {
    double r2 = r * r;
    return -0.5 * p.mu * p.mu * r2 + 0.25 * p.lambda * r2 * r2 + std::pow(p.mu, 4) / (4 * p.lambda);
}
double dV(const Params& p, double r) { return -p.mu * p.mu * r + p.lambda * r * r * r; }
double d2V(const Params& p, double r) { return -p.mu * p.mu + 3 * p.lambda * r * r; }

// V(phi) pour un vecteur complet de composantes
double V_vec(const Params& p, const Vec& phi) {
    double r2 = 0;
    for (double x : phi) r2 += x * x;
    return V_radial(p, std::sqrt(r2));
}

// Potentiel à température T (approximation haute température, ordre T^2)
double V_thermal(const Params& p, double r, double T) {
    return V_radial(p, r) + 0.5 * p.thermal_c() * T * T * r * r;
}

// ───────────────────────────── Utilitaires ─────────────────────────────

Vec linspace(double a, double b, int n) {
    Vec x(n);
    for (int i = 0; i < n; ++i) x[i] = n > 1 ? a + (b - a) * i / (n - 1) : a;
    return x;
}

void section(const std::string& title) {
    std::cout << "\n=== " << title << " ===\n";
}

struct Check {
    int passed = 0, failed = 0;
    void operator()(bool ok, const std::string& what) {
        (ok ? passed : failed)++;
        std::cout << (ok ? "  [OK]   " : "  [ECHEC] ") << what << "\n";
    }
};

void write_csv(const std::string& path, const std::string& header, const std::vector<Vec>& cols) {
    if (path.empty()) return;
    std::ofstream f(path);
    if (!f) throw std::runtime_error("impossible d'ouvrir " + path);
    f << std::setprecision(10) << header << "\n";
    for (size_t i = 0; i < cols[0].size(); ++i) {
        for (size_t j = 0; j < cols.size(); ++j) f << (j ? "," : "") << cols[j][i];
        f << "\n";
    }
    std::cout << "  -> données écrites dans " << path << "\n";
}

// Courbe ASCII y(x). Les x sont uniformément répartis.
void ascii_plot(const Vec& x, const Vec& y, const std::string& xlabel, int W = 70, int H = 18) {
    double ymin = *std::min_element(y.begin(), y.end());
    double ymax = *std::max_element(y.begin(), y.end());
    if (ymax - ymin < 1e-300) ymax = ymin + 1;
    std::vector<std::string> g(H, std::string(W, ' '));
    // axe y=0 éventuel
    if (ymin < 0 && ymax > 0) {
        int r0 = (int)std::lround((ymax - 0) / (ymax - ymin) * (H - 1));
        for (int c = 0; c < W; ++c) g[r0][c] = '-';
    }
    // axe x=0 (phi = 0)
    if (x.front() < 0 && x.back() > 0) {
        int c0 = (int)std::lround((0 - x.front()) / (x.back() - x.front()) * (W - 1));
        for (int r = 0; r < H; ++r) g[r][c0] = (g[r][c0] == '-') ? '+' : '|';
    }
    for (int c = 0; c < W; ++c) {
        double xi = x.front() + (x.back() - x.front()) * c / (W - 1);
        // interpolation linéaire
        double t = (xi - x.front()) / (x.back() - x.front()) * (x.size() - 1);
        size_t i = std::min<size_t>((size_t)t, x.size() - 2);
        double yi = y[i] + (y[i + 1] - y[i]) * (t - i);
        int r = (int)std::lround((ymax - yi) / (ymax - ymin) * (H - 1));
        r = std::clamp(r, 0, H - 1);
        g[r][c] = '*';
    }
    std::cout << std::fixed << std::setprecision(3);
    for (int r = 0; r < H; ++r) {
        double yl = ymax - (ymax - ymin) * r / (H - 1);
        if (r % 4 == 0 || r == H - 1) std::cout << std::setw(9) << yl << " ";
        else std::cout << std::string(10, ' ');
        std::cout << g[r] << "\n";
    }
    std::cout << std::string(10, ' ') << std::left << std::setw(W / 2) << x.front() << std::right
              << std::setw(W - W / 2) << x.back() << "\n"
              << std::string(10, ' ') << xlabel << "\n";
    std::cout << std::defaultfloat;
}

// ───────────────────────────── Modules ─────────────────────────────

void show_parameters(const Params& p) {
    section("Paramètres et prédictions analytiques");
    std::cout << std::setprecision(6)
              << "  N (composantes)        = " << p.N << "\n"
              << "  mu                     = " << p.mu << "\n"
              << "  lambda                 = " << p.lambda << "\n"
              << "  vev  v = mu/sqrt(lam)  = " << p.v() << "\n"
              << "  m_Higgs = sqrt(2) mu   = " << p.mh() << "\n"
              << "  m_Goldstone            = 0   (" << p.N - 1 << " mode(s))\n"
              << "  T_c                    = " << p.Tc() << "\n"
              << "  V(0) - V(v)            = " << V_radial(p, 0) << "  (profondeur du puits)\n";
    if (p.N == 1)
        std::cout << "  (symétrie Z2 discrète : pas de boson de Goldstone, 2 vides +v / -v)\n";
    else
        std::cout << "  (symétrie continue O(" << p.N << ") -> O(" << p.N - 1
                  << ") : le vide est une sphère de rayon v)\n";
}

void potential_table(const Params& p) {
    section("Potentiel V(phi) le long d'un axe");
    double v = p.v();
    Vec phi = linspace(-p.range * v, p.range * v, p.points), V(p.points);
    for (int i = 0; i < p.points; ++i) V[i] = V_radial(p, std::fabs(phi[i]));
    ascii_plot(phi, V, "phi  (minima en +/- v = " + std::to_string(v).substr(0, 6) + ")");
    std::cout << "\n" << std::setw(12) << "phi" << std::setw(16) << "V(phi)" << std::setw(16) << "V'(phi)"
              << std::setw(16) << "V''(phi)" << "\n";
    int step = std::max(1, p.points / 20);
    for (int i = 0; i < p.points; i += step) {
        double r = phi[i];
        double s = r < 0 ? -1 : 1;
        std::cout << std::setw(12) << r << std::setw(16) << V[i] << std::setw(16) << s * dV(p, std::fabs(r))
                  << std::setw(16) << d2V(p, std::fabs(r)) << "\n";
    }
    write_csv(p.csv, "phi,V", {phi, V});
}

// Vue en dessus du chapeau mexicain (plan phi1-phi2)
void mexican_hat(const Params& p) {
    section("Chapeau mexicain (plan phi1-phi2, vue de dessus)");
    const char* shades = " .:-=+*#%@";
    double v = p.v(), R = 1.8 * v;
    double vmax = V_radial(p, R);
    int W = 61, H = 29;
    for (int j = 0; j < H; ++j) {
        double y = R - 2 * R * j / (H - 1);
        std::string line;
        for (int i = 0; i < W; ++i) {
            double x = -R + 2 * R * i / (W - 1) * 1.0;
            double r = std::hypot(x, y * 1.0);
            double t = std::clamp(V_radial(p, r) / vmax, 0.0, 1.0);
            char c = shades[(int)(t * 9.999)];
            if (std::fabs(r - v) < 0.06 * v) c = 'O';  // cercle des vides
            line += c;
        }
        std::cout << "  " << line << "\n";
    }
    std::cout << "  'O' = ensemble des vides |phi| = v ; centre = maximum local instable (V = "
              << V_radial(p, 0) << ")\n";
}

// Minimisation numérique de V par Newton + vérification des masses par Hessien discret
void vacuum_and_masses(const Params& p, Check& chk) {
    section("Vide, masses et spectre (vérification numérique)");
    // Newton sur V'(r) = 0
    double r = 1.5 * p.v();  // au-delà du point d'inflexion v/sqrt3 : converge vers le minimum
    int it = 0;
    for (; it < 100; ++it) {
        double step = dV(p, r) / d2V(p, r);
        r -= step;
        if (std::fabs(step) < 1e-14) break;
    }
    std::cout << std::setprecision(10) << "  Newton : v_num = " << r << "   (analytique " << p.v()
              << ", " << it << " itérations)\n";
    chk(std::fabs(r - p.v()) < 1e-9 * p.v(), "vev numérique = mu/sqrt(lambda)");

    // Hessien par différences finies centrées dans toutes les directions
    Vec phi(p.N, 0.0);
    phi[0] = r;
    double h = 1e-4 * p.v();
    Vec masses2(p.N);
    for (int i = 0; i < p.N; ++i) {
        Vec a = phi, b = phi;
        a[i] += h;
        b[i] -= h;
        masses2[i] = (V_vec(p, a) - 2 * V_vec(p, phi) + V_vec(p, b)) / (h * h);
    }
    std::cout << "  Masses (valeurs propres du Hessien de V au vide) :\n";
    for (int i = 0; i < p.N; ++i)
        std::cout << "    mode " << i << " : m^2 = " << std::setw(14) << masses2[i]
                  << (i == 0 ? "   <- Higgs radial" : "   <- Goldstone") << "\n";
    chk(std::fabs(std::sqrt(std::max(masses2[0], 0.0)) - p.mh()) < 1e-5 * p.mh(), "m_Higgs = sqrt(2) mu");
    for (int i = 1; i < p.N; ++i) chk(std::fabs(masses2[i]) < 1e-6 * p.mu * p.mu, "Goldstone sans masse (mode " + std::to_string(i) + ")");
    chk(d2V(p, 0) < 0, "V''(0) < 0 : l'origine est instable (brisure spontanée)");
    std::cout << std::setprecision(6);
}

// Masse des bosons W, Z et du fermion via le couplage au Higgs (cas Modèle Standard)
void gauge_masses(const Params& p, double g, double gp, double yt) {
    section("Masses engendrées par le mécanisme de Higgs (Modèle Standard)");
    double v = p.v();
    std::cout << "  m_W = g v / 2              = " << g * v / 2 << "\n"
              << "  m_Z = sqrt(g^2+g'^2) v / 2 = " << std::sqrt(g * g + gp * gp) * v / 2 << "\n"
              << "  m_photon                   = 0\n"
              << "  m_fermion = y v / sqrt(2)  = " << yt * v / std::sqrt(2.0) << "   (y = " << yt << ")\n"
              << "  rho = m_W^2/(m_Z^2 cos^2 th_W) = 1 à l'arbre\n";
}

// Oscillation du champ homogène autour du vide : phi'' = -V'(phi). Intégrateur RK4.
void oscillation(const Params& p, Check& chk) {
    section("Dynamique : oscillation autour du vide (le quantum = boson de Higgs)");
    double v = p.v();
    double r = v * 1.02, w = 0.0, t = 0.0;  // petite excitation
    double dt = 2e-3 / p.mu;
    auto acc = [&](double x) { return -dV(p, x); };
    auto rk4 = [&](double& x, double& u) {
        double k1x = u, k1u = acc(x);
        double k2x = u + 0.5 * dt * k1u, k2u = acc(x + 0.5 * dt * k1x);
        double k3x = u + 0.5 * dt * k2u, k3u = acc(x + 0.5 * dt * k2x);
        double k4x = u + dt * k3u, k4u = acc(x + dt * k3x);
        x += dt / 6 * (k1x + 2 * k2x + 2 * k3x + k4x);
        u += dt / 6 * (k1u + 2 * k2u + 2 * k3u + k4u);
    };
    double E0 = 0.5 * w * w + V_radial(p, r);
    double prev = r - v, last_cross = -1;
    Vec periods;
    Vec ts, rs;
    for (int n = 0; n < 2000000 && periods.size() < 12; ++n) {
        rk4(r, w);
        t += dt;
        double cur = r - v;
        if (prev < 0 && cur >= 0) {  // passage ascendant par v
            double tc = t - dt * cur / (cur - prev);
            if (last_cross >= 0) periods.push_back(tc - last_cross);
            last_cross = tc;
        }
        prev = cur;
        if (n % 50 == 0) { ts.push_back(t); rs.push_back(r); }
    }
    double T = 0;
    for (double q : periods) T += q;
    T /= std::max<size_t>(1, periods.size());
    double omega = 2 * M_PI / T;
    double E1 = 0.5 * w * w + V_radial(p, r);
    std::cout << "  amplitude initiale = 2% de v\n"
              << "  période mesurée   T = " << T << "  ->  omega = " << omega << "\n"
              << "  m_Higgs analytique   = " << p.mh() << "  (écart relatif "
              << std::fabs(omega - p.mh()) / p.mh() << ", dû à l'anharmonicité)\n"
              << "  dérive de l'énergie  = " << std::fabs(E1 - E0) / E0 << "\n";
    chk(std::fabs(omega - p.mh()) / p.mh() < 0.01, "fréquence d'oscillation = m_Higgs (à 1%)");
    chk(std::fabs(E1 - E0) / E0 < 1e-8, "conservation de l'énergie (RK4)");
    ascii_plot(ts, rs, "temps", 70, 12);
    write_csv(p.csv, "t,phi", {ts, rs});
}

// Mur de domaine 1+1D : champ sur un réseau, évolution leapfrog, kink boosté.
void domain_wall(const Params& p, Check& chk) {
    section("Mur de domaine (kink) en 1+1 dimensions");
    double mu = p.mu, v = p.v();
    double Mexact = 2 * std::sqrt(2.0) / 3 * std::pow(mu, 3) / p.lambda;
    const double beta = 0.5, gamma = 1 / std::sqrt(1 - beta * beta);
    const double L = 80.0 / mu;
    const int n = 4000;
    const double dx = L / n, dt = 0.25 * dx;
    Vec x(n + 1), f(n + 1), pi(n + 1), force(n + 1);
    double x0 = -0.2 * L;
    for (int i = 0; i <= n; ++i) {
        x[i] = -L / 2 + i * dx;
        double u = gamma * mu * (x[i] - x0) / std::sqrt(2.0);
        double th = std::tanh(u), sech2 = 1 - th * th;
        f[i] = v * th;
        pi[i] = -v * beta * gamma * mu / std::sqrt(2.0) * sech2;  // d phi/dt du kink boosté
    }
    auto energy = [&]() {
        double E = 0;
        for (int i = 0; i < n; ++i) {
            double g = (f[i + 1] - f[i]) / dx;
            double fm = 0.5 * (f[i] + f[i + 1]), pm = 0.5 * (pi[i] + pi[i + 1]);
            E += (0.5 * pm * pm + 0.5 * g * g + V_radial(p, std::fabs(fm))) * dx;
        }
        return E;
    };
    auto compute_force = [&]() {
        for (int i = 1; i < n; ++i)
            force[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (dx * dx) - (-mu * mu * f[i] + p.lambda * f[i] * f[i] * f[i]);
    };
    auto position = [&]() {
        for (int i = 0; i < n; ++i)
            if (f[i] < 0 && f[i + 1] >= 0) return x[i] - f[i] * dx / (f[i + 1] - f[i]);
        return 0.0;
    };
    double E0 = energy(), xs = position();
    int steps = (int)(0.4 * L / (beta) / dt * 0.5);  // le kink parcourt ~0.2 L
    compute_force();
    for (int s = 0; s < steps; ++s) {  // vitesse-Verlet
        for (int i = 1; i < n; ++i) pi[i] += 0.5 * dt * force[i];
        for (int i = 1; i < n; ++i) f[i] += dt * pi[i];
        compute_force();
        for (int i = 1; i < n; ++i) pi[i] += 0.5 * dt * force[i];
    }
    double E1 = energy(), xe = position(), tt = steps * dt;
    double vel = (xe - xs) / tt;
    std::cout << std::setprecision(6)
              << "  énergie au repos analytique (2sqrt2/3) mu^3/lambda = " << Mexact << "\n"
              << "  énergie initiale (boost beta=" << beta << ")  = " << E0 << "   (gamma*M = " << gamma * Mexact << ")\n"
              << "  énergie finale                         = " << E1 << "\n"
              << "  vitesse mesurée du mur                = " << vel << "   (attendue " << beta << ")\n";
    chk(std::fabs(E0 - gamma * Mexact) / (gamma * Mexact) < 1e-3, "énergie du kink = gamma * (2sqrt2/3) mu^3/lambda");
    chk(std::fabs(E1 - E0) / E0 < 1e-3, "conservation de l'énergie sur le réseau");
    chk(std::fabs(vel - beta) < 1e-2, "le kink se propage à la vitesse du boost (invariance de Lorentz)");
    Vec xs_plot, fs_plot;
    for (int i = 0; i <= n; i += n / 100) { xs_plot.push_back(x[i]); fs_plot.push_back(f[i]); }
    std::cout << "  profil final phi(x) (transition entre les deux vides -v et +v) :\n";
    ascii_plot(xs_plot, fs_plot, "x", 70, 10);
}

// Transition de phase thermique : restauration de la symétrie à haute température
void thermal(const Params& p, Check& chk) {
    section("Restauration de la symétrie à haute température");
    double Tc = p.Tc(), c = p.thermal_c();
    std::cout << "  m^2(T) = -mu^2 + c T^2,  c = (N+2) lambda/24 = " << c << "   ->  T_c = " << Tc << "\n\n"
              << std::setw(10) << "T/Tc" << std::setw(14) << "v(T)" << std::setw(14) << "v(T)/v" << "   phase\n";
    Vec Ts, vs;
    int n = std::max(p.points, 21);
    for (int i = 0; i < n; ++i) {
        double T = 1.4 * Tc * i / (n - 1);
        double m2 = -p.mu * p.mu + c * T * T;
        double vt = m2 < 0 ? std::sqrt(-m2 / p.lambda) : 0.0;
        Ts.push_back(T);
        vs.push_back(vt);
        if (i % std::max(1, n / 14) == 0)
            std::cout << std::setw(10) << T / Tc << std::setw(14) << vt << std::setw(14) << vt / p.v()
                      << (vt > 0 ? "   brisée" : "   symétrique") << "\n";
    }
    ascii_plot(Ts, vs, "T   (v(T) s'annule en T_c : transition du second ordre)", 70, 10);
    // vérification numérique : minimum de V_T trouvé par balayage
    auto argmin = [&](double T) {
        double best = 0, bv = V_thermal(p, 0, T);
        for (int i = 1; i <= 4000; ++i) {
            double r = 2 * p.v() * i / 4000;
            double val = V_thermal(p, r, T);
            if (val < bv) { bv = val; best = r; }
        }
        return best;
    };
    chk(std::fabs(argmin(0.5 * Tc) - p.v() * std::sqrt(1 - 0.25)) < 1e-3 * p.v(), "T = 0.5 Tc : minimum en v*sqrt(1-(T/Tc)^2)");
    chk(argmin(1.1 * Tc) == 0.0, "T = 1.1 Tc : minimum à l'origine (symétrie restaurée)");
    write_csv(p.csv, "T,vev", {Ts, vs});
}

// ───────────────────────────── CLI ─────────────────────────────

void usage() {
    std::cout <<
R"(Usage : higgs [commande] [options]

Commandes :
  all          (défaut) toutes les sections ci-dessous
  info         paramètres et prédictions analytiques
  potential    table et courbe de V(phi)            (--csv fichier pour exporter)
  hat          chapeau mexicain en ASCII (plan phi1-phi2)
  vacuum       minimisation de V, masses du Higgs et des Goldstones
  oscillation  dynamique du champ autour du vide, mesure de m_Higgs
  kink         mur de domaine 1+1D sur réseau (kink boosté)
  thermal      transition de phase thermique, T_c
  test         exécute toutes les vérifications ; code de sortie 1 si échec

Options :
  --mu X       paramètre de masse (défaut 1)
  --lambda X   auto-couplage > 0 (défaut 0.5)
  --n N        nombre de composantes réelles du champ : 1, 2, 4... (défaut 1)
  --points K   nombre de points des courbes (défaut 61)
  --range R    plage du champ en unités de v (défaut 2)
  --sm         préréglage Modèle Standard : v = 246.22 GeV, m_h = 125.25 GeV, N = 4
  --csv FILE   exporte les données de la commande dans FILE
  -h, --help   cette aide
)";
}

int main(int argc, char** argv) {
    Params p;
    std::string cmd = "all";
    bool sm = false;
    try {
        auto need = [&](int& i) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string("valeur manquante après ") + argv[i]);
            return argv[++i];
        };
        auto num = [&](const std::string& s) {
            size_t pos;
            double d = std::stod(s, &pos);
            if (pos != s.size()) throw std::runtime_error("nombre invalide : " + s);
            return d;
        };
        bool first = true;
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a == "-h" || a == "--help") { usage(); return 0; }
            else if (a == "--mu") p.mu = num(need(i));
            else if (a == "--lambda") p.lambda = num(need(i));
            else if (a == "--n") p.N = (int)num(need(i));
            else if (a == "--points") p.points = (int)num(need(i));
            else if (a == "--range") p.range = num(need(i));
            else if (a == "--csv") p.csv = need(i);
            else if (a == "--sm") sm = true;
            else if (first && a[0] != '-') cmd = a;
            else throw std::runtime_error("option inconnue : " + a);
            first = false;
        }
        if (sm) {
            const double v = 246.22, mh = 125.25;
            p.N = 4;
            p.mu = mh / std::sqrt(2.0);
            p.lambda = mh * mh / (2 * v * v);
        }
        if (p.lambda <= 0 || p.mu <= 0) throw std::runtime_error("mu et lambda doivent être > 0");
        if (p.N < 1) throw std::runtime_error("--n doit être >= 1");
        if (p.points < 5) throw std::runtime_error("--points doit être >= 5");
        if (p.range <= 0) throw std::runtime_error("--range doit être > 0");

        Check chk;
        static const std::vector<std::string> known = {"all", "info", "potential", "hat", "vacuum", "oscillation", "kink", "thermal", "test"};
        if (std::find(known.begin(), known.end(), cmd) == known.end())
            throw std::runtime_error("commande inconnue : " + cmd + " (voir --help)");

        std::cout << "Boson de Higgs et brisure spontanée de symétrie\n";
        bool all = cmd == "all" || cmd == "test";
        if (all || cmd == "info") show_parameters(p);
        if (sm && (all || cmd == "info")) gauge_masses(p, 0.652, 0.357, 0.99);
        if (all || cmd == "potential") potential_table(p);
        if (all || cmd == "hat") mexican_hat(p);
        if (all || cmd == "vacuum") vacuum_and_masses(p, chk);
        if (all || cmd == "oscillation") oscillation(p, chk);
        if (all || cmd == "kink") domain_wall(p, chk);
        if (all || cmd == "thermal") thermal(p, chk);

        if (chk.passed + chk.failed > 0) {
            std::cout << "\nVérifications : " << chk.passed << " réussies, " << chk.failed << " échouées\n";
            if (chk.failed) return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "erreur : " << e.what() << "\n";
        return 2;
    }
    return 0;
}
