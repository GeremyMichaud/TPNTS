#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from vpython import *
import time
from statistics import mean
from scipy.stats import maxwell

def DrudeModelSimulator(iterations=1000, champ=0):
    # Simulation parameters
    N_electrons = 200
    N_cores = 64  # It is easier when you choose a squared number for the number of cores
    dt = 1E-7  # Increase the time increment if you dare! Brace your for a spectacular display of fireworks.
    electron2follow = 0
    L = 1

    # Physical constants
    #mass_electron = 9E-31
    #mass_core = 4E-3 / 6E23  # Helium mass
    mass = 4E-3 / 6E23  # Helium mass
    R_electron = 0.01
    R_core = 0.025
    K = 1.4E-23
    T = 300
    q = 1.602e-19

    start_time = time.time()

    # Canvas initialization
    animation = canvas(width=750, height=500, align='center')
    animation.range = L
    animation.title = 'Simulation du modèle de Drude'
    s = """
    Simulation de particules modélisées en sphères pour représenter le modèle de Drude.
    Les sphères représentent les électrons mobiles (gris) et les coeurs ioniques fixes (rouges).
    """
    animation.caption = s
    d = L/2 + max(R_core, R_electron)
    cadre = curve(color=color.blue, radius=0.005)
    cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0),
                    vector(-d, d, 0), vector(-d, -d, 0)])

    # Initialize electrons
    Electrons = []
    p = []
    apos = []
    pavg = sqrt(2 * mass * 1.5 * K * T)

    for i in range(N_electrons):
        x = L * random() - L / 2
        y = L * random() - L / 2
        z = 0
        if i == electron2follow:
            Electrons.append(simple_sphere(pos=vector(x, y, z), radius=0.02, color=color.magenta))
        Electrons.append(simple_sphere(pos=vector(x, y, z), radius=R_electron, color=color.gray(0.7)))
        apos.append(vec(x, y, z))
        phi = 2 * pi * random()
        px = pavg * cos(phi)
        py = pavg * sin(phi)
        pz = 0
        p.append(vector(px, py, pz))

    # Initialize ion cores
    Cores = []
    cpos = []
    ion_core_spacing = L / (sqrt(N_cores) - 1)

    for i in range(int(sqrt(N_cores))):
        for j in range(int(sqrt(N_cores))):
            # Create a 2D grid of ion cores with periodic arrangement within the square
            x = i * ion_core_spacing - L / 2
            y = j * ion_core_spacing - L / 2
            z = 0
            Cores.append(simple_sphere(pos=vector(x, y, z), radius=R_core, color=color.red))
            cpos.append(vec(x, y, z))

    def check_collisions():
        """Check for collisions.

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

    pavg_list = []
    posavg_list = []

    # Main simulation
    for _ in range(iterations):
        rate(300)

        pavg = mean([pp.mag for pp in p])
        posavg = mean([pos.y for pos in apos])

        # Move electrons
        for i in range(N_electrons):
            vitesse = p[i] / mass
            deltax = vitesse * dt
            Electrons[i].pos = apos[i] = apos[i] + deltax

        # Conserve momentum with wall collisions
        for i in range(N_electrons):
            loc = apos[i]
            if abs(loc.x) > L / 2:
                p[i].x = abs(p[i].x) if loc.x < 0 else -abs(p[i].x)
            if abs(loc.y) > L / 2:
                p[i].y = abs(p[i].y) if loc.y < 0 else -abs(p[i].y)
                p[i].y += q * champ * dt

        hitlist = check_collisions()

        pavg_list.append(pavg)
        posavg_list.append(posavg)

        temperature = pavg**2 / (3 * K * mass)
        # Handle electron-core collisions
        for ij in hitlist:
            i = ij[1]
            j = ij[0]
            posi = apos[i]
            posj = cpos[j]
            vi = p[i] / mass
            rrel = hat(posi - posj)
            vrel = vi

            dx = dot(posi - posj, vrel.hat)
            dy = cross(posi - posj, vrel.hat).mag
            alpha = asin(dy / (R_electron + R_core))
            d = (R_electron + R_core) * cos(alpha) - dx
            deltat = d / vrel.mag

            # CHANGE L'INTERPÉNÉTRATION DES SPHÈRES PAR LA CINÉTIQUE DE COLLISION
            angle = pi * (random() - 1 / 2)
            norm_p = mass * maxwell.rvs(scale=sqrt(K * temperature / mass))
            p[i] = norm_p * hat(rotate(rrel, angle=angle))
            apos[i] = posi + (posi - posj) * deltat

        """for ion, electron in hitlist:
            mtot = mass_electron + mass_core
            posi = cpos[ion]
            posj = apos[electron]
            Vcom = p[electron] / mtot
            vj = p[electron] / mass_electron
            vi = p[ion] / mass_core
            rrel = hat(posj - posi)

            if vj.mag2 == 0 or rrel.mag > R_electron:
                continue

            dx = dot(posj-posi, vj.hat)
            dy = cross(posj-posi, vj.hat).mag
            alpha = asin(dy / (R_core + R_electron))
            d = (R_core + R_electron) * cos(alpha) - dx
            deltat = d / vj.mag

            theta = 2 * pi * random()
            mb_mean_velocity = sqrt(K * temperature / mass_electron)
            norm = maxwell.rvs(scale=mb_mean_velocity)
            p[electron] = norm * hat(rotate(rrel, angle=theta))
            apos[electron] = posj + (posj - posi) * deltat
            #collision_velocity_x = norm * cos(theta)
            #collision_velocity_y = norm * sin(theta)
            #collision_velocity_z = 0
            #collision_velocity = vector(collision_velocity_x, collision_velocity_y, collision_velocity_z)
            #p[electron] = mass_electron * collision_velocity

            #posj = posj - vj * deltat
            #pcomj = p[electron] - mass_electron * Vcom
            #rrel = hat(rrel)
            #pcomj = pcomj - 2 * dot(pcomj, rrel) * rrel
            #p[electron] = pcomj + mass_electron * Vcom
            #apos[electron] = posj + (p[electron] / mass_electron) * deltat

        # Calculate the average momentum magnitude of electrons
        #p_list = [p[ii].mag for ii in range(N_electrons)]
        #p_average = mean(p_list)

        # Calculate the average y-position of electrons
        #pos_list = [pos.y for pos in apos]
        #pos_average = mean(pos_list)

        #pavg_list.append(p_average)
        #posavg_list.append(pos_average)"""

    print("--- %.3f seconds ---" % (time.time() - start_time))
    return [pavg_list, posavg_list]

DrudeModelSimulator()