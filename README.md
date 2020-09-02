# Transposition de matrices creuses

- [Transposition de matrices creuses](#transposition-de-matrices-creuses)
  - [Installation](#installation)
  - [Matrices](#matrices)
  - [Utilisation](#utilisation)
    - [Intel Math Kernel Library (MKL)](#intel-math-kernel-library-mkl)
    - [Intel Threading Building Block (TBB)](#intel-threading-building-block-tbb)
  - [Version](#version)
  - [TODO](#todo)
  - [Remarques](#remarques)

## Installation

Testé avec:

- Python 3.8.3: [PSF licence agreement](https://docs.python.org/3/license.html)
- Matplotlib 3.2.1: [The Matplotlib Development Team](https://matplotlib.org/3.2.1/users/license.html)
- NumPy 1.18.5: [Copyright (c) 2005, NumPy Developers](https://numpy.org/doc/stable/license.html)
- pandas 1.0.3: [BSD 3-Clause License](https://pandas.pydata.org/pandas-docs/stable/getting_started/overview.html)

- spTrans: [Wang, Hao and Liu, Weifeng and Hou, Kaixi and Feng, Wu-chun, LGPL-2.1](https://github.com/vtsynergy/sptrans)
- SpaSM: [The SpaSM group, GPL-3.0](https://github.com/cbouilla/spasm/)
- Intel Math Kernel Library (MKL) 2020.0.166: [Intel Simplified Software License](https://software.intel.com/content/www/us/en/develop/articles/end-user-license-agreement.html#inpage-nav-3)
- Intel Threading Building Block (TBB) 2020.0.166: [Intel Simplified Software License](https://software.intel.com/content/www/us/en/develop/articles/end-user-license-agreement.html#inpage-nav-3)
- ICC
- C 11
- C 11 standard librairies
- C++ 11
- C++ 11 STL
- GCC 6.3.0 20170516 (Debian 6.3.0-18+deb9u1) et 10.1.0

## Matrices

Les tailles des matrices sont disponibles [ici](matrices.md).

## Utilisation

### Intel Math Kernel Library (MKL)

Pour utiliser la MKL sur ppti-gpu-1 :

`source /usr/intel/mkl/bin/mklvars.sh intel64 [ilp64]`

La dernière option (qui est facultative) active l'usage d'entiers 64 bits. On peut essayer pour voir si ça résout nos problèmes...

Il faut l'activer dans le Makefile (`USE_TBB=yes`), puis faire `make clean` et recompiler.

Après, il y a plusieurs options (séquentiel, parallèle avec OpenMP, etc). Choisir la bonne dans le Makefile. J'ai utilisé le [Intel® Math Kernel Library Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html) pour les déterminer...

### Intel Threading Building Block (TBB)

Pour utiliser TBB sur ppti-gpu-1 :

`source /usr/intel/tbb/bin/tbbvars.sh intel64`

Ensuite, `make clean` puis recompiler.

## Version

1. count + prefix-sum dans transpose()
2. finalize dans transpose()
3. lo_bucket -> hi_bucket dans transpose_bucket()
4. WCB et OUT en struct entry_t *
5. 2 lignes de cache + OMP sur Rp avec n = 1
6. 3. et OMP sur Rp avec n = 1
7. 6. et maximisation du radix
8. 4. et écriture de Rj et Rx dans la boucle de Rp
9. 8. et test/print
10. 3. avec variation du radix selon nnz

## TODO

- [x] afficher la version de MKL
- [x] afficher la version de TBB
- [x] afficher le nombre de threads utilisés avec MKL
- [x] afficher le nombre de threads utilisés avec TBB
- [x] trouver le remplaçant de task_scheduler_init
- [x] créer un nouveau makefile
- [x] changer les noms des fichiers
- [x] changer la hierarchie des fichiers
- [x] initialiser time et Total
- [x] créer fonction triplet_transpose
- [x] fixer le nombre de threads de omp for avec num_threads
- [x] donner num_threads à MKL
- [x] gérer les benchmarks
- [x] créer un csv
- [x] fonctionne avec gcc/g++ ? OUI
- [x] fonctionne avec gcc/g++ et tbb ? OUI, mais fuites mémoire non gérables dues à TBB
- [x] fonctionne avec gcc/g++ et MKL ? OUI
- [x] fonctionne avec icc/icpc? OUI
- [x] fonctionne avec icc/icpc et tbb ? OUI
- [x] fonctionne avec icpc et mkl ? OUI
- [x] mettre dans un sous-dossier, makefile qui créer deux executables.
- [x] créer regle make benchmarks
- [x] modifications des entetes fallacieux des matrices RSA.
- [x] utiliser valgrind pour corriger les fuites mémoire et les erreurs de driver_wang : sptrans -> segfault, vec -> segfault?
- [x] mettre un assert dans la fonction load en supposant que l'entête est faux à la place de spasm_max
- [x] afficher dans le csv lire use_mkl,compiler_flags
- [x] mettre Have_mkl dans le name dans le csv
- [x] enlever le chemin d'accès des matrices dans le csv
- [x] prendre les algo parallel, moyenner les algo parallel, parcourir les matrices, parcourir les moyennes d'algo parallel (choisir un algo et parcourir le même sur les threads), tracer transpose en fonction de threads, idem pour transpose_tr, récupérer thread = argmax
- [x] parcourir les matrices, parcourir les algo, tracer les boites à moustaches
- [x] améliorer la ligne de compilation (les sauts de lignes dus à -D)
- [x] traiter le csv et afficher des graphes (boxplot pour la variance, pour parallèle : courbe de l'accélération), dataframe pandas
- [x] créer un sous-dossier dans charts automatiquement
- [x] pour la durée : regarder le min sur la moyenne des deux ?
- [x] sauvegarder les fichiers avec des _ au lieux des ' '
- [x] mettre une ligne avec plus de pointillés et des symboles plus petits
- [x] pour résumer : prendre l'addition des medianes par algo sur toutes les matrices
- [x] graphique accélération, calculer l'accélération avec un algo séquentiel pour chaque matrice : Gustavson -O2
- [x] utiliser typedef u32 pour n et m et u32/u64 pour nz
- [x] finalize,spasm_add_entry, spasm_human_format mis dans .h et sans static
- [x] reformater le code
- [x] à remettre: la vérification avec la MKL, 
- [x] trouver le min de la moyenne des deux transpositions dans find_minima
- [x] décrire le fonctionnement de tbb:parallel_sort
- [x] décrire le fonctionnement de std::sort
- [ ] essayer avec O2 et O3 pour chaque algo, O3 defavorables ? Gustavon OK, std::sort OK, tbb::sort OK, MKL , scan, merge
- [ ] citation/licence SpaSM/CADO-NFS ou refactor ?
- [ ] renommer les versions, garder 3 et 10
- [ ] commenter le code et refaire la doc
- [ ] améliorer les prints
- [ ] checker les todo
- [ ] parallèliser et vectoriser avec AVX128 AVX256 AVX512
- [ ] renommer les fichiers contenant des "::"
- [ ] refactor code
- [ ] calculer l'occupation mémoire de scanTrans avec 1 thread par rapport à celle de Gustavson
- [ ] radix sort LSD parallèle
- [ ] radix sort MSD parallèle
- [ ] radix sort MSD+LSD tester sur toutes, utiliser un autre tri lorsque le nombre d'éléments par bucket est petit, demi flush/purge, réduire les branchements conditionnels

## Remarques

std::sort est plus rapide que classical sur pre-transpose[6,7,8,9,10,11,12]. Son écart de durée est d'environ 5 ms (< 10ms). La variabilité inter-algo est alors plus petite que la variablité intra-algo.

Lorsque classical transpose et transpose_tr sont relativement proches, std transpose est plus lent que std transpose_tr

MKL iomp:
 s'arrete au nombre de coeurs physique (44)

MKL iomp 2PV:
 est idéal entre 8 et 24 threads (voire souvent entre 8 et 15, se confirme sur le graphique grobal). Pour les matrices RSA: vers 8 pour les premières puis vers 12 pour les dernières. Pour les matrices de Wang, la courbe atteint parfois un minimum plus tard vers 20-24.
 grande différence selon la forme du rectangle
 plus la matrice est (petite ou dense (première de pre-transpose)) plus le passage à l'echelle est mauvais

MergeTrans:
	est idéal entre 12 et 32 (souvent vers 16, plutôt entre 20 et 32 sur le graphique global)
	peu de différence selon la forme du rectangle sur les matrices de Wang, très grandes différence sur les 10 premières de pre-transpose (pas le même algo ?)

ScanTrans:
	ressemble à MergeTrans, mais le passage à l'échelle semble meilleur car le minimum est atteint un peu plus loin (la décroissance est plus applati jusqu'à 40 et la croissance est moins semble forte)
	sur les 10 premières de pre-transpose les deux transpositions utilisent le même algo
	encore moins de différences entre les deux transpositions que MergeTrans

dans la code de Wang, les asserts ne fonctionne pas avec scantrans^2 et les pre-transpose car celles-ci sont triées par ligne et non par ligne puis par colonne. Ainsi, la sortie colIdx de scantrans qui est triée par ligne puis par colonne n'est pas dans le même ordre que dans le fichier .mtx

TBB:
	sur les 10 premières de pre-transpose, il se passe quelque chose à partir de 80 et grand écart entre les deux transpositions
	courbe très plate, minimum atteint vers 16 (ou 44 sur les dernières)

ScanTrans est globalement le meilleur algo (overall duration et speed up)

Radix sort 3 meilleur que 10 :
 - circuit
 - flickr
 - transient
 - la tr de wiki-talk
 - les tr de pre_transpose[20..53] (peut-être aussi de 15 à 19)
 - inversion des courbes sur pre_transpose[54..58]

ScanTrans meilleur que 10:
 - rajat29
 - stomach
 - circuit
 - Stanford tr
 - torso
 - venkat
 - wiki-talk tr
