#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

const double mu = -1.0; // Paramètre de masse
const double lambda = 0.5; // Paramètre d'auto-interaction
const int points = 100; // Nombre de points pour la simulation

// Fonction pour calculer le potentiel de Higgs
double higgs_potential(double phi) {
    return mu * mu * phi * phi + lambda * phi * phi * phi * phi;
}

// Fonction pour calculer les valeurs du potentiel
std::vector<double> calculate_potential(std::vector<double>& phi_values) {
    std::vector<double> potential_values;
    for (double phi : phi_values) {
        potential_values.push_back(higgs_potential(phi));
    }
    return potential_values;
}

// Fonction pour afficher les résultats
void print_results(const std::vector<double>& phi_values, const std::vector<double>& potential_values) {
    std::cout << std::setw(10) << "phi" << std::setw(20) << "V(phi)" << std::endl;
    for (size_t i = 0; i < phi_values.size(); ++i) {
        std::cout << std::setw(10) << phi_values[i] << std::setw(20) << potential_values[i] << std::endl;
    }
}

int main() {
    // Générer les valeurs du champ phi
    std::vector<double> phi_values(points);
    for (int i = 0; i < points; ++i) {
        phi_values[i] = -2.0 + 4.0 * i / (points - 1); // Valeurs de -2 à 2
    }

    // Calculer les valeurs du potentiel de Higgs
    std::vector<double> potential_values = calculate_potential(phi_values);

    // Afficher les résultats
    print_results(phi_values, potential_values);

    return 0;
}