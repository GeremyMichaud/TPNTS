#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 15:40:27 2021

#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood
# Claudine Allen
"""

from vpython import *
import time
import sys# Initialisation des variables (la plupart sont inchangées)
# Initialisation des variables (la plupart sont inchangées)
start_time = time.time()
Natoms = 200
dt = 1E-5
mass = 4E-3/6E23
Ratom = 0.01
k = 1.4E-23
T = 300
L = 1

# Initialisation des électrons
electrons = []
p = []
pavg = sqrt(2*mass*1.5*k*T) 

# Ajout des électrons dans la simulation
for i in range(Natoms):
    x, y = L*random()-L/2, L*random()-L/2
    phi = 2*pi*random()
    px, py = pavg*cos(phi), pavg*sin(phi)
    electrons.append(sphere(pos=vec(x,y,0), radius=Ratom, color=color.blue))
    p.append(vec(px,py,0))

# Ajout des cœurs ioniques fixes
Nions = 49  # Nombre de cœurs ioniqCorrues à définir
ions = []
for i in range(Nions):
    x = ((i % 7)-3) * 15 * Ratom  # Répartition périodique
    y = ((i // 7)-3) * 15 * Ratom
    ions.append(sphere(pos=vec(x, y, 0), radius=Ratom*5, color=color.red))

# Modification de checkCollisions pour ignorer les collisions entre électrons
def checkCollisions():
    hitlist = []
    for electron in electrons:
        for ion in ions:
            dr = electron.pos - ion.pos
            if mag2(dr) < (10*Ratom)**2:
                hitlist.append(electron)
    return hitlist


def g_vi(v, m, T):
    """
    Calcule la distribution de Maxwell-Boltzmann pour la vitesse v.
    """
    kB = 1.38e-23  # Constante de Boltzmann (m^2 kg s^-2 K^-1)
    coef = (m / (2 * pi * kB * T)) ** 0.5
    exponent = exp(-m * v**2 / (2 * kB * T))
    return coef * exponent

# Simulation des collisions inélastiques
def handle_inelastic_collisions(hitlist):
    for electron in hitlist:
        speed = sqrt(3*k*T/mass)  # v = sqrt(3kT/m)
        phi = 2*pi*random()
        electron.p = g_vi(speed, mass, T)*cos(phi), g_vi(speed, mass, T)*sin(phi)

# Boucle principale de la simulation modifiée
for _ in range(500):
    hitlist = checkCollisions()
    handle_inelastic_collisions(hitlist)
    for i, electron in enumerate(electrons):
        electron.pos += p[i]/mass * dt
        # Gérer les collisions avec les parois
        if abs(electron.pos.x) > L/2: p[i].x *= -1
        if abs(electron.pos.y) > L/2: p[i].y *= -1
        

print("--- %.3f seconds ---" % (time.time() - start_time))