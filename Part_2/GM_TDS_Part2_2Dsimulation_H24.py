# Created on Tuesday Feb 26 2024

# DrudeModelSimulator

# Gérémy Michaud

from vpython import *
import numpy as np
import time
from scipy.stats import maxwell


def DrudeModelSimulator(iterations=500, champ=0, dt=1E-7):
    # Increase the time increment dt if you dare! Brace yourself for a spectacular display of fireworks.

    # Simulation parameters
    N_electrons = 200
    N_cores = 64  # It is easier when you choose a squared number for the number of cores
    electron2follow = 0

    # Physical constants
    mass_electron = 9E-31
    R_electron = 0.01  # wildly exaggerated size of an atom
    R_core = 0.025
    k = 1.4E-23
    T = 300
    q = 1.602E-19

    start_time = time.time()

    # CANEVAS DE FOND
    L = 1
    animation = canvas(width=750, height=500, align='center')
    animation.range = L
    animation.title = 'Simulation du modèle de Drude'
    s = """
    Simulation de particules modélisées en sphères pour représenter le modèle de Drude.
    Les sphères représentent les électrons mobiles (gris) et les coeurs ioniques fixes (rouges).
    """
    animation.caption = s
    d = L / 2 + max(R_core, R_electron)
    cadre = curve(color=color.blue, radius=0.005)
    cadre.append([vector(-d, -d, 0), vector(d, -d, 0), vector(d, d, 0),
                vector(-d, d, 0), vector(-d, -d, 0)])

    # Initialize electrons
    Electrons = []
    p = []
    apos = []
    pavg = sqrt(3 * mass_electron * k * T)

    for i in range(N_electrons):
        x = L * random() - L / 2
        y = L * random() - L / 2
        z = 0
        if i == electron2follow:
            Electrons.append(simple_sphere(pos=vector(x, y, z), radius=0.02, color=color.magenta))
        else:
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

    def checkCollisions():
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
    pmag_electron = [mag(p[electron2follow])]

    # Main simulation
    for i in range(iterations):
        rate(300)

        pavg = np.mean([nb_p.mag for nb_p in p])
        #posavg = np.mean([pos.y for pos in apos])

        # Move electrons
        for i in range(N_electrons):
            vitesse = p[i] / mass_electron
            deltax = vitesse * dt
            Electrons[i].pos = apos[i] = apos[i] + deltax

        # Conserve momentum with wall collisions
        for i in range(N_electrons):
            loc = apos[i]
            if abs(loc.x) > L / 2:
                p[i].x = abs(p[i].x) if loc.x < 0 else -abs(p[i].x)
            if abs(loc.y) > L / 2:
                p[i].y = abs(p[i].y) if loc.y < 0 else -abs(p[i].y)
                #p[i].y += q * champ * dt

        hitlist = checkCollisions()

        pavg_list.append(pavg)
        #posavg_list.append(posavg)

        temperature = pavg**2 / (3 * k * mass_electron)  # Redefine new temperature with lost of momentum (inelatic)

        # Handle electron-core collisions
        for ij in hitlist:
            i = ij[1]
            j = ij[0]
            posi = apos[i]
            posj = cpos[j]
            vi = p[i] / mass_electron
            rrel = hat(posi - posj)
            vrel = vi

            deltat = d / vrel.mag

            theta = 2 * pi * random()
            norm_p = mass_electron * maxwell.rvs(scale=sqrt(k * temperature / mass_electron))
            p[i] = norm_p * hat(rotate(rrel, angle=theta))
            apos[i] = posi + (posi - posj) * deltat

            if i == electron2follow:
                pmag_electron.append(mag(p[i]))

    print("--- %.3f seconds ---" % (time.time() - start_time))

    return [pavg_list, pmag_electron]

DrudeModelSimulator()
