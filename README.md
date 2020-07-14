# Transposition de matrices creuses

## Matrices

Le jeu de matrices utilisé par *Wang et. al* est disponible sur la machine ppti-gpu-1 dans le dossier :

	/Infos/lmd/2019/master/ue/MU4IN903-2020fev

Taille : 5Go.

Un autre jeu de matrices (le mien) est disponible dans le dossier :

	/Infos/lmd/2019/master/ue/MU4IN903-2020fev/RSA

Taille : 14Go. Ces matrices sont plus rectangulaires.


## Intel Math Kernel Library (MKL)

Version : 2020.0.166

Pour utiliser la mkl sur ppti-gpu-1 :

	source /usr/intel/mkl/bin/mklvars.sh intel64 [ilp64]  <--- si on veut des entiers 64 bits

La dernière option (qui est facultative) active l'usage d'entiers 64 bits. On peut essayer pour voir si ça résoud nos problèmes...

Il faut l'activer dans le Makefile (`USE_TBB=yes`), puis faire `make clean` et recompiler. 

Après, il y a plusieurs options (séquentiel, parallèle avec OpenMP, etc). Choisir la bonne dans le Makefile. J'ai utilisé le [Intel® Math Kernel Library Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html) pour les déterminer...


## Intel Threading Building Block (TBB)

Version : 2020.0.166

Pour utiliser tbb sur ppti-gpu-1 :

	source /usr/intel/tbb/bin/tbbvars.sh intel64

Ensuite, `make clean` puis recompiler.

## GCC

Version : 6.3.0 20170516 (Debian 6.3.0-18+deb9u1)

## TODO

-[x] afficher la version de MKL
-[x] afficher la version de TBB
-[x] afficher le nombre de threads utilisés avec MKL
-[] afficher le nombre de threads utilisés avec TBB
-[x] trouver le remplaçant de task_scheduler_init
-[] commenter le code
-[] reformater le code
-[] changer les noms des fichiers
-[] changer la hierarchie des fichiers
-[] gérer les benchmarks
-[] décrire le fonctionnement de tbb:parallel_sort
-[] décrire le fonctionnement de std::sort

## Algorithme de tri

Intel TBB parallel_sort : std::sort parallel + std::sort sequentiel avec seuil à 500
std::sort : introsort (quicksort + heap sort avec seuil à ?) + insertion sort avec seuil à 16.

## Liens utiles

