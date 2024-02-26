#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from vpython import *
import time
import numpy as np
from statistics import mean
from scipy.stats import maxwell

start_time = time.time()

# Physical constants
mass_electron = 9E-31
mass_core = 4E-3 / 6E23  # Helium mass
R_electron = 0.01
R_core = 0.025
K = 1.4E-23
T = 300

# Simulation parameters
N_electrons = 200
N_cores = 36
dt = 1E-7  # Increase the time increment if you dare! Brace your for a spectacular display of fireworks.
iterations = 1000
electron2follow = 0

# Canvas initialization
L = 1
gray = color.gray(0.7)
blue = color.blue
red = color.red
pink = color.magenta
animation = canvas(width=750, height=500, align='center')
animation.range = L
animation.title = 'Simulation du modèle de Drude'
s = """
Simulation de particules modélisées en sphères pour représenter le modèle de Drude.
Les sphères représentent les électrons mobiles (gris) et les coeurs ioniques fixes (rouges).
"""
animation.caption = s
d = L/2 + max(R_core, R_electron)
cadre = curve(color=blue, radius=0.005)
cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0),
            vector(-d, d, 0), vector(-d, -d, 0)])

def initElectrons():
    """Initialize electrons.

    Returns:
        list: Electrons
        list: Momenta (p)
        list: Positions (apos)
    """
    Electrons = []
    p = []
    apos = []
    pavg = sqrt(2 * mass_electron * 1.5 * K * T)

    for i in range(N_electrons):
        x = L * random() - L / 2
        y = L * random() - L / 2
        z = 0
        if i == electron2follow:
            Electrons.append(simple_sphere(pos=vector(x, y, z), radius=0.02, color=pink))
        Electrons.append(simple_sphere(pos=vector(x, y, z), radius=R_electron, color=gray))
        apos.append(vec(x, y, z))
        phi = 2 * pi * random()
        px = pavg * cos(phi)
        py = pavg * sin(phi)
        pz = 0
        p.append(vector(px, py, pz))

    return Electrons, p, apos

def initCores():
    """Initialize ion cores.

    Returns:
        list: Core positions (cpos)
    """
    Cores = []
    cpos = []
    ion_core_spacing = L / (sqrt(N_cores) + 1)

    for i in range(int(sqrt(N_cores))):
        for j in range(int(sqrt(N_cores))):
            # Create a 2D grid of ion cores with periodic arrangement within the square
            x = (i + 1) * ion_core_spacing - L / 2
            y = (j + 1) * ion_core_spacing - L / 2
            z = 0
            Cores.append(simple_sphere(pos=vector(x, y, z), radius=R_core, color=red))
            cpos.append(vec(x,y,z))

    return cpos

def checkCollisions(cpos, apos):
    """Check for collisions.

    Args:
        cpos (list): Core positions
        apos (list): Electron positions

    Returns:
        list: List of collisions
    """
    hitlist = []
    r2 = (R_electron + R_core) ** 2

    for i in range(len(cpos)):
        ci = cpos[i]
        for j in range(len(apos)):
            aj = apos[j]
            # Does not care for electron-electron collisions
            dr = ci - aj
            if mag2(dr) < r2:
                hitlist.append([i, j])

    return hitlist

def simulate():
    """Main simulation function.

    Returns:
        list: Average momenta over time
        list: Average y positions over time
    """
    Electrons, p, apos = initElectrons()
    cpos = initCores()

    t = 0
    time = []
    pavg_list = []
    posavg_list = []
    ln_list = []
    tau_list = []

    for _ in range(iterations):
        rate(300)

        moveElectrons(Electrons, p, apos)
        conserveMomentum(p, apos)
        hitlist = checkCollisions(cpos, apos)
        collideElectronCore(hitlist, p, apos, cpos)

        p_average = calculateAverageMomentum(p)
        pos_average = calculateAverageposition(apos)
        ln = np.log(p_average / (sqrt(2 * mass_electron * 1.5 * K * T)))
        tau = -t / np.log(p_average / (sqrt(2 * mass_electron * 1.5 * K * T)))

        time.append(t)
        pavg_list.append(p_average)
        posavg_list.append(pos_average)
        ln_list.append(ln)
        tau_list.append(tau)
        t += dt

    return pavg_list, posavg_list

def moveElectrons(Electrons, p, apos):
    """Move electrons.

    Args:
        Electrons (list): Electron objects
        p (list): Momenta
        apos (list): Positions
    """
    for i in range(N_electrons):
        vitesse = p[i] / mass_electron
        deltax = vitesse * dt
        Electrons[i].pos = apos[i] = apos[i] + deltax

def conserveMomentum(p, apos):
    """Conserve momentum with wall collisions.

    Args:
        p (list): Momenta
        apos (list): Positions
    """
    for i in range(N_electrons):
        loc = apos[i]
        if abs(loc.x) > L / 2:
            if loc.x < 0: p[i].x =  abs(p[i].x)
            else: p[i].x =  -abs(p[i].x)
        if abs(loc.y) > L / 2:
            if loc.y < 0: p[i].y = abs(p[i].y)
            else: p[i].y = -abs(p[i].y)

def collideElectronCore(hitlist, p, apos, cpos):
    """Handle electron-core collisions.

    Args:
        hitlist (list): List of collisions
        p (list): Momenta
        apos (list): Positions
        cpos (list): Core positions
    """
    for ion, electron in hitlist:
        mtot = mass_electron + mass_core
        posi = cpos[ion]
        posj = apos[electron]
        Vcom = p[electron] / mtot
        vj = p[electron] / mass_electron
        rrel = posi - posj

        if vj.mag2 == 0 or rrel.mag > R_electron:
            continue

        theta = 2 * pi * random()
        mb_mean_velocity = sqrt((8 * K * T) / (pi * mass_electron))
        norm = maxwell.rvs(scale=mb_mean_velocity)
        collision_velocity_x = norm * cos(theta)
        collision_velocity_y = norm * sin(theta)
        collision_velocity_z = 0
        collision_velocity = vector(collision_velocity_x, collision_velocity_y, collision_velocity_z)
        p[electron] = mass_electron * collision_velocity

        dx = dot(rrel, vj.hat)
        dy = cross(rrel, vj.hat).mag
        alpha = asin(dy / (R_core + R_electron))
        d = (R_core + R_electron) * cos(alpha) - dx
        deltat = d / vj.mag

        posj = posj - vj * deltat
        pcomj = p[electron] - mass_electron * Vcom
        rrel = hat(rrel)
        pcomj = pcomj - 2 * dot(pcomj, rrel) * rrel
        p[electron] = pcomj + mass_electron * Vcom
        apos[electron] = posj + (p[electron] / mass_electron) * deltat

def calculateAverageMomentum(p):
    """
    Calculate the average momentum magnitude of electrons.

    Args:
        p (list): List of vectors representing the momenta of electrons.

    Returns:
        float: The average momentum magnitude.
    """
    p_list = [p[ii].mag for ii in range(N_electrons)]
    p_average = mean(p_list)
    return p_average

def calculateAverageposition(epos):
    """
    Calculate the average y-position of electrons.

    Args:
        epos (list): List of vectors representing the positions of electrons.

    Returns:
        float: The average y-position.
    """
    pos_list = [pos.y for pos in epos]
    pos_average = mean(pos_list)
    return pos_average

if __name__ == "__main__":
    result = simulate()
    print("--- %.3f seconds ---" % (time.time() - start_time))