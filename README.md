# Gravitational Simulations

Simulating the three-body problem in C++, visualized in the browser with Three.js.

**[Live Demo](https://mugiln.github.io/three-body-problem/)**

---

## The Three-Body Problem

Take three masses in space pulling on each other through gravity. How do they move?

There is no general formula. Two bodies (Earth and Sun) produce a clean ellipse. Add a third and the motion becomes impossible to predict analytically. A tiny change in starting position produces a completely different trajectory a few orbits later. This is one of the earliest known examples of chaos in physics.

Most configurations just fly apart. The scenarios below are the rare exceptions where something structured happens.

---

## Scenarios

**Figure-Eight** -- Three equal masses chasing each other in a stable loop. Discovered by Chenciner and Montgomery in 2000, one of the very few known periodic solutions to the three-body problem.

**Lagrange** -- Three equal masses at the corners of an equilateral triangle, rotating around their shared center of mass. Found by Lagrange in 1772.

**Chaotic** -- The figure-eight with one body nudged slightly out of the orbital plane. Starts looking nearly periodic, then breaks down. A straightforward demonstration of sensitive dependence on initial conditions.

**Euler** -- Three bodies on a straight line, rotating around their center of mass. Euler's 1767 solution. Valid but unstable; any small disturbance causes it to drift, which shows up in the trail pattern.

**Restricted Three-Body** -- The Sun and Jupiter are fixed on a circular orbit. A spacecraft sits in their combined gravitational field. This scenario computes what that field looks like as a surface, with five equilibrium points (L1 through L5) marked on it. L1, L2, L3 lie on the Sun-Jupiter axis. L4 and L5 sit 60 degrees ahead and behind Jupiter. The Trojan asteroids cluster at L4 and L5 because those two points are gravitationally stable.

---

## How It Works

### Physics (C++)

Uses the velocity Verlet integration method, which conserves energy significantly better than basic Euler integration.

Each time step:
1. Compute acceleration on each body from all others: `a = G * m / r^2`
2. Update positions: `x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2`
3. Recompute acceleration at new positions
4. Update velocities: `v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt`

Each scenario outputs a CSV file of trajectory data.

For the restricted three-body scenario, a separate file computes the effective gravitational potential on a 300x300 grid in the rotating reference frame:

```
phi = -G*M_sun/r1 - G*M_jup/r2 - 0.5*omega^2*(x^2 + y^2)
```

The third term accounts for the centrifugal effect of the rotating frame. 

### Visualizer (Three.js)

A single `index.html` handles everything:

- Glowing bodies with fading trail lines
- Starfield background
- Drag to rotate, scroll to zoom
- Adjustable playback speed
- Restricted three-body view renders the potential as a colored 3D surface with L1-L5 labeled directly on it

---

## Project Structure

```
three-body-problem/
├── index.html
├── README.md
├── .gitattributes
├── src/
│   ├── threebodyproblem.cpp
│   └── 3bodyandrestricted.py
└── data/
    ├── figure-eight.csv
    ├── lagrange.csv
    ├── chaotic.csv
    ├── euler.csv
    └── potential.csv
```

---

## Running Locally

Compile and run the C++ simulation:

```bash
g++ -O2 -o sim threebodyproblem.cpp
./sim
```

Serve the visualizer:

```bash
python -m http.server 8000
```

Open `http://localhost:8000` in your browser.

---

## Stack

- C++ for physics simulation and potential computation
- Three.js for 3D WebGL rendering


---

## References

- Chenciner & Montgomery (2000), figure-eight solution
- Lagrange (1772), equilateral triangle solution
- Euler (1767), collinear solution
