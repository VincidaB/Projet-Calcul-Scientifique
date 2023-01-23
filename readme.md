# Projet de Calcul Scientifique Numérique
# Etude de la réponse d'un thermocouple par la méthode des différences finies

Ci-joint se trouve le fichier `simulation.py`, les variables pouvant être modifiées  sont présentent au début du programme et sont :
* `T0` la température initiale de la sphère
* `Tinfty` la température à l'infini
* `a` le rayon du thermocouple (de la sphère)
* `K` le facteur de convection
* `Nr`,`Ntheta` et `Nphi` me nombre de mailles selon $e_r$, $e_\theta$ et $e_\varphi$ respectivement
* `rho` la masse volumique du matériau constituant le thermocouple
* `Cp` la chaleur latente de ce matériau
* `Lambda` la conductivité thermique de ce matériau
* `dt` le pas temporel (exprimé en secondes)
* `NB_iterations` le nombres d'itérations de la simulation 

L'animation présentée durera `NB_iterations/(1000*10)` secondes et sera exporté avec le nom `simulation.gif` dans le dossier depuis lequel aura été exécuté la simulation.

Belpois Vincent
Chomette Hugo