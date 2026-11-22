# exact-exp32 (en français)

Exponentielles en simple précision exactes (`exp2f`, `expf`, `exp10f`) sans test d'arrondi sur processeur, un `exp2f` exact calculé uniquement en float sur carte graphique NVIDIA, et une recherche par programmation génétique (à la Koza) de code exact plus rapide. Chaque candidat est vérifié sur les 2^32 entrées.

Auteur : Maher El Hmissi, novembre 2026. English: [README.md](README.md).

## Résultats

| Partie | Contenu | Vérification | Vitesse (honnêtement) |
| --- | --- | --- | --- |
| `cpu/exp2f.c`, `cpu/expf.c`, `cpu/exp10f.c` | Fichiers au format CORE-MATH : table de 64, polynôme de degré 5, **sans test d'arrondi ni chemin lent** | 0 erreur sur les 2^32 entrées × 4 modes d'arrondi, contre CORE-MATH, avec gcc et clang | Pas de gain net avec l'outil officiel de CORE-MATH sur un Ryzen 7 7800X3D dédié : −9 % (gcc), +9 à +10 % (clang), latence identique |
| `cpu/logf.c` | Même idée pour `logf`, avec un contrôle léger pour 5 entrées qu'aucun calcul en double ne peut trancher | 0 erreur, 2^32 × 4 modes | Pas de gain net |
| `gpu/fexp2.h` | `exp2f` exact (arrondi au plus près) **calculé uniquement en float**, pour les cartes faibles en double précision | 0 erreur sur les 2^32 entrées, sur une RTX 5090 | ×1,02 le temps de l'`exp2f` de CUDA à 1 appel par nombre lu, ×1,17 à 4, ×2,88 à 16, ×5,7 en calcul pur |
| `koza/` | Programmation génétique linéaire sur GPU ; candidats vérifiés sur toutes les entrées, contre-exemples réinjectés | Gagnant exact sur les 2^32 entrées | 15 opérations au lieu de 21 ; 7,2 % plus rapide que la meilleure version écrite à la main (objectif : 10 %) |

L'`exp2f` de CUDA se trompe sur 189 192 458 entrées sur 4 294 967 296 (4,40 %), avec des écarts jusqu'à 2 ulp.

## Ce qui existait déjà

Ce dépôt est une reproduction indépendante, pas une première :

- **Arm (optimized-routines)** publie des `expf`, `exp2f` et `exp10f` exacts construits avec la même recette, calculés en double.
- **LLVM libc** contient un `expf` exact (arrondi au plus près) calculé en float et conçu pour les GPU. Leur `exp2f` en float n'est pas encore exact en novembre 2026 ; celui de `gpu/fexp2.h` l'est.
- **Intel** (arXiv 2605.15547) obtient des `expf` exacts sur ses GPU, surtout en double ; le calcul en 32 bits y est un travail futur.
- Le `log2f` de CORE-MATH est déjà sans test d'arrondi ; RLIBM produit des fonctions exactes par ajustement sous contraintes.

L'apport du dépôt : les données (entrées difficiles, vérifications exhaustives, mesures officielles sur machine dédiée), un `exp2f` exact en float pour CUDA, et la recherche génétique sur GPU. Elle a trouvé une étape compensée en 3 opérations : `sh = fma(th, ph, th); d = sh - th; err ≈ fma(th, ph, -d)`.

## Une leçon sur les mesures

Un banc de mesure maison sur une machine partagée annonçait d'abord 17 à 26 % de gain sur processeur. L'outil officiel de CORE-MATH, sur une machine dédiée, n'a montré aucun gain net. Il faut toujours mesurer avec l'outil des responsables du projet, sur une machine dédiée, avant d'annoncer un gain.

Les commandes pour tout reproduire sont dans [README.md](README.md).

## À propos du moteur de recherche

Le moteur de programmation génétique sur GPU (à la Koza : programmation génétique linéaire, vérification exhaustive, contre-exemples réinjectés) est l'outil de recherche privé de l'auteur ; il n'est pas publié ici. Les programmes qu'il a trouvés, leur re-vérification exhaustive et leurs chronométrages sont publiés.
