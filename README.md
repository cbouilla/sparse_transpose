# matrices

Le jeu de matrices utilisé par *Wang et. al* est disponible sur la machine ppti-gpu-1 dans le dossier :

	/Infos/lmd/2019/master/ue/MU4IN903-2020fev

Taille : 5Go.

Un autre jeu de matrices (le mien) est disponible dans le dossier :

	/Infos/lmd/2019/master/ue/MU4IN903-2020fev/RSA

Taille : 14Go. Ces matrices sont plus rectangulaires.


# Intel Math Kernel Library (MKL)

Pour utiliser la mkl sur ppti-gpu-1 :

	source /usr/intel/mkl/bin/mklvars.sh intel64 [ilp64]  <--- si on veut des entiers 64 bits

La dernière option (qui est facultative) active l'usage d'entiers 64 bits. On peut essayer pour voir si ça résoud nos problèmes...

Il faut l'activer dans le Makefile (`USE_TBB=yes`), puis faire `make clean` et recompiler. 

Après, il y a plusieurs options (séquentiel, parallèle avec OpenMP, etc). Choisir la bonne dans le Makefile. J'ai utilisé le [Intel® Math Kernel Library Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html) pour les déterminer...


# Intel Threading Building Block (TBB)

Pour utiliser tbb sur ppti-gpu-1 :

	source /usr/intel/tbb/bin/tbbvars.sh intel64

Ensuite, `make clean` puis recompiler.