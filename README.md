# Higgs — brisure spontanée de symétrie en C++

Programme pédagogique et vérifié numériquement illustrant le mécanisme de Higgs.
Champ scalaire à `N` composantes réelles (`N=1` : Z2, `N=2` : chapeau mexicain U(1), `N=4` : doublet du Modèle Standard) :

```
V(φ) = −½ μ² |φ|² + ¼ λ |φ|⁴ + μ⁴/(4λ)        (constante : V = 0 dans le vide)
```

> Correction par rapport à la première version : `V = μ²φ² + λφ⁴` avec `μ² > 0` n'a qu'un seul minimum et ne brise aucune symétrie ; le terme quadratique doit être **négatif**.

## Résultats analytiques utilisés

| Grandeur | Formule |
|---|---|
| vev | `v = μ/√λ` |
| masse du Higgs | `m_h = √2 μ` |
| Goldstones | `N−1` modes de masse nulle |
| température critique | `T_c² = 24 μ² / ((N+2) λ)` |
| mur de domaine (N=1) | `φ(x) = v tanh(μx/√2)`, énergie `(2√2/3) μ³/λ` |
| W, Z, fermion (`--sm`) | `g v/2`, `√(g²+g'²) v/2`, `y v/√2` |

## Compilation et utilisation

```
make            # ou : g++ -O2 -std=c++17 higgs.cpp -o higgs
./higgs         # toutes les sections
./higgs --help
```

| Commande | Contenu |
|---|---|
| `info` | paramètres, vev, masses, `T_c` |
| `potential` | table + courbe ASCII de V, V', V'' (`--csv f.csv` pour exporter) |
| `hat` | chapeau mexicain vu de dessus |
| `vacuum` | minimisation (Newton) et masses par Hessien discret : Higgs massif, Goldstones nuls |
| `oscillation` | dynamique RK4 autour du vide ; la fréquence mesurée redonne `m_h` |
| `kink` | mur de domaine 1+1D sur réseau, kink boosté : énergie `γM`, vitesse, conservation |
| `thermal` | `v(T)`, transition du second ordre, restauration de la symétrie |
| `test` | toutes les vérifications ; code de sortie 1 en cas d'échec |

Options : `--mu`, `--lambda`, `--n`, `--points`, `--range`, `--csv`, `--sm` (v = 246,22 GeV, m_h = 125,25 GeV, N = 4).

Exemples :

```
./higgs test --n 2
./higgs info --sm
./higgs thermal --lambda 0.13 --mu 2 --csv thermal.csv
```

## Limites

Théorie classique (arbre) ; la partie thermique utilise l'approximation haute température à l'ordre `T²`. Les corrections quantiques (potentiel de Coleman–Weinberg, boucles) ne sont pas incluses.
