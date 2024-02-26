# Created on Fri Feb 5 15:40:27 2021
# GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen

from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.stats import maxwell

def Simulation(champ, iterations):
    # Déclaration de variables influençant le temps d'exécution de la simulation
    Natoms = 100  # change this to have more or fewer atoms
    Ncoeur_coter = 8
    dt = 1E-5  # pas d'incrémentation temporel

    # Déclaration de variables physiques "Typical values"
    mass = 4E-3/6E23  # helium mass
    Ratom = 0.01  # wildly exaggerated size of an atom
    Rcoeur = 0.03
    k = 1.4E-23  # Boltzmann constant # TODO: changer pour une constante de Boltzmann en eV/K :)
    T = 300  # around room temperature
    q = 1.602e-19

    # CANEVAS DE FOND
    L = 1  # container is a cube L on a side
    gray = color.gray(0.7)  # color of edges of container and spheres below
    animation = canvas(width=750, height=500)  # , align='left')
    animation.range = L

    # ARÊTES DE BOÎTE 2D
    d = L/2+Ratom
    r = 0.005
    cadre = curve(color=gray, radius=r)
    cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0), vector(-d, d, 0), vector(-d, -d, 0)])

    # POSITION ET QUANTITÉ DE MOUVEMENT INITIALE DES SPHÈRES
    Atoms = []  # Objet qui contiendra les sphères pour l'animation
    Coeur = []  # Objet qui contiendra les
    p = []  # quantité de mouvement des sphères
    apos = []  # position des sphères
    cpos = []  # position des coeurs
    pavg = sqrt(2*mass*1.5*k*T)  # average kinetic energy p**2/(2mass) = (3/2)kT

    # Création de coordonée pour les sphère fixes
    for i in range(Ncoeur_coter):  # Nombre de ligne
        for n in range(Ncoeur_coter):  # Nombre de ligne
            x = i*L/(Ncoeur_coter-1)-L/2  # position d'un coeur en tenant compte de l'origine au centre de la boîte
            y = n*L/(Ncoeur_coter-1)-L/2
            z = 0
            Coeur.append(simple_sphere(pos=vector(x, y, z), radius=Rcoeur, color=color.red))
            cpos.append(vec(x, y, z))

    for i in range(Natoms):
        x = L*random()-L/2  # position aléatoire qui tient compte que l'origine est au centre de la boîte
        y = L*random()-L/2
        z = 0
        if i == 0:
            Atoms.append(simple_sphere(pos=vector(x, y, z), radius=0.03, color=color.magenta))
        else:
            Atoms.append(simple_sphere(pos=vector(x, y, z), radius=Ratom, color=gray))
        apos.append(vec(x, y, z))

        phi = 2*pi*random()  # direction aléatoire pour la quantité de mouvement
        px = pavg*cos(phi)  # qte de mvt initiale selon l'équipartition
        py = pavg*sin(phi)
        pz = 0
        p.append(vector(px, py, pz))  # liste de la quantité de mouvement initiale de toutes les sphères

    # FONCTION POUR IDENTIFIER LES COLLISIONS
    def checkCollisions():
        hitlist = []  # initialisation
        r2 = Ratom + Rcoeur  # distance critique où les 2 sphères entre en contact à la limite de leur rayon
        r2 *= r2
        for i in range(Ncoeur_coter**2):
            ci = cpos[i]
            for j in range(Natoms):
                aj = apos[j]
                dr = ci - aj
                if mag2(dr) < r2:
                    hitlist.append([i, j])
        return hitlist

    # BOUCLE PRINCIPALE POUR L'ÉVOLUTION TEMPORELLE DE PAS dt
    it = 0
    collision = 0
    liste_p = []
    liste_p_avg = []
    liste_pos_avg = []

    for i in range(iterations):
        rate(100)  # limite la vitesse de calcul de la simulation pour que l'animation soit visible à l'oeil humain!

        pavg = np.mean(np.array([nb_p.mag for nb_p in p]))  # remodifie la valeur du p average pour chaque itéreation
        posavg = np.mean(np.array([pos.y for pos in apos]))

        # DÉPLACE TOUTES LES SPHÈRES D'UN PAS SPATIAL deltax
        vitesse = []  # vitesse instantanée de chaque sphère
        deltax = []  # pas de position de chaque sphère correspondant à l'incrément de temps dt
        for i in range(Natoms):
            vitesse.append(p[i]/mass)
            deltax.append(vitesse[i] * dt)
            Atoms[i].pos = apos[i] = apos[i] + deltax[i]

        # CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS AVEC LES MURS DE LA BOÎTE
        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L/2:
                if loc.x < 0:
                    p[i].x = abs(p[i].x)
                else:
                    p[i].x = -abs(p[i].x)
            if abs(loc.y) > L/2:
                if loc.y < 0:
                    p[i].y = abs(p[i].y)
                else:
                    p[i].y = -abs(p[i].y)
                p[i].y += q * champ * dt

        # LET'S FIND THESE COLLISIONS!!!
        hitlist = checkCollisions()

        liste_p_avg.append(pavg)
        liste_pos_avg.append(posavg)
        collision += len(hitlist)

        temperature = pavg**2/(3*k*mass)
        # CONSERVE LA QUANTITÉ DE MOUVEMENT AUX COLLISIONS ENTRE SPHÈRES
        for ij in hitlist:
            i = ij[1]
            j = ij[0]
            posi = apos[i]
            posj = cpos[j]
            vi = p[i]/mass
            rrel = hat(posi-posj)
            vrel = vi

            dx = dot(posi-posj, vrel.hat)
            dy = cross(posi-posj, vrel.hat).mag
            alpha = asin(dy/(Ratom+Rcoeur))
            d = (Ratom+Rcoeur)*cos(alpha)-dx
            deltat = d/vrel.mag

            # CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION
            angle = pi*(random()-1/2)
            norm_p = mass*maxwell.rvs(scale=sqrt(k*temperature/mass))
            p[i] = norm_p * hat(rotate(rrel, angle=angle))
            apos[i] = posi+(posi-posj)*deltat

            liste_p.append(p[i].mag)

        it += 1

    return [liste_p_avg, liste_pos_avg]

Simulation(0, 700)
