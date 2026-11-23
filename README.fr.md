# exact-exp32 (en français)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24485607.svg)](https://doi.org/10.5281/zenodo.24485607)

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

La partie processeur est une reproduction indépendante, pas une première :

- **Arm (optimized-routines)** publie des `expf`, `exp2f` et `exp10f` exacts construits avec la même recette, calculés en double.
- **LLVM libc** contient un `expf` exact (arrondi au plus près) calculé en float et conçu pour les GPU. Leur `exp2f` en float n'est pas encore exact en novembre 2026 ; celui de `gpu/fexp2.h` l'est.
- **Intel** (arXiv 2605.15547) obtient des `expf` exacts sur ses GPU, surtout en double ; le calcul en 32 bits y est un travail futur.
- Le `log2f` de CORE-MATH est déjà sans test d'arrondi ; RLIBM produit des fonctions exactes par ajustement sous contraintes.
- Un `expf` calculé uniquement en entiers, arrondi de façon statique, a été proposé pour LLVM libc en août 2026 (PR #209406) : une autre façon d'éviter le calcul en double, pour `expf` et non `exp2f`.
- La recherche évolutive de code plus rapide existe déjà, par exemple GEVO (optimisation évolutive de noyaux GPU) et STOKE (super-optimisation stochastique de code en nombres flottants, à précision réglable). Aucun des deux n'exige que chaque candidat soit exact sur toutes les entrées.

À notre connaissance (recherche faite en novembre 2026), `gpu/fexp2.h` est le premier `exp2f` correctement arrondi calculé uniquement en float pour CUDA, et `koza/` la première recherche génétique sur GPU où chaque candidat est vérifié sur les 2^32 entrées.

L'apport du dépôt : les données (entrées difficiles, vérifications exhaustives, mesures officielles sur machine dédiée), un `exp2f` exact en float pour CUDA, et la recherche génétique sur GPU. Elle a trouvé une étape compensée en 3 opérations : `sh = fma(th, ph, th); d = sh - th; err ≈ fma(th, ph, -d)`.

## Une leçon sur les mesures

Un banc de mesure maison sur une machine partagée annonçait d'abord 17 à 26 % de gain sur processeur. L'outil officiel de CORE-MATH, sur une machine dédiée, n'a montré aucun gain net. Il faut toujours mesurer avec l'outil des responsables du projet, sur une machine dédiée, avant d'annoncer un gain.

Les commandes pour tout reproduire sont dans [README.md](README.md).

## Voir aussi

[tt-erf-fp32](https://github.com/Maher75017/tt-erf-fp32) : la fonction `erf` FP32 pour le SFPU de Tenstorrent (Blackhole, Wormhole), à 1 ULP au plus sur toutes les entrées, vérifiée avec la même méthode exhaustive sur un modèle bit à bit de l'arithmétique du matériel. Proposée à tt-metal dans [tenstorrent/tt-metal#66809](https://github.com/tenstorrent/tt-metal/pull/66809).

## À propos du moteur de recherche

Le moteur de programmation génétique sur GPU (à la Koza : programmation génétique linéaire, vérification exhaustive, contre-exemples réinjectés) est l'outil de recherche privé de l'auteur ; il n'est pas publié ici. Les programmes qu'il a trouvés, leur re-vérification exhaustive et leurs chronométrages sont publiés.

## Licence

MIT (voir `LICENSE`). Les fichiers issus du projet CORE-MATH (par exemple `cpu/*.c`, les copies de référence `cm_*_ref.c` et `gpu/exp2f_orig.c`) conservent leur mention MIT d'origine et le copyright (c) 2022-2023 Alexei Sibidanov.
