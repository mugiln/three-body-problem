# gravitational simulations

simulating the three-body problem in C++ and rendering it in the browser with Three.js.

**[live demo](https://mugiln.github.io/three-body-problem/)**

---

## what is the three-body problem?

take three masses in space. they all pull on each other through gravity. how do they move?

turns out there is no general formula. unlike two bodies (earth and sun, clean ellipse, done), three bodies produce motion that is impossibly sensitive to starting conditions. change one position by a millimeter and you get a completely different trajectory a few orbits later. this is one of the earliest known examples of chaos theory.

most of the time the system just flies apart. the scenarios below are the rare cases where something interesting happens.

---

## scenarios

**figure-eight** -- three equal masses chasing each other in a stable figure-eight loop. discovered by Chenciner and Montgomery in 2000. one of the very few known periodic solutions.

**lagrange** -- three equal masses at the corners of an equilateral triangle, rotating around their shared center of mass. Lagrange found this in 1772.

**chaotic** -- the figure-eight but with body 3 nudged slightly out of the plane. starts looking almost normal, then falls apart. classic butterfly effect demonstration.

**euler** -- three bodies in a straight line, rotating around their center of mass. Euler's 1767 solution. technically valid but unstable in practice, any small push and it drifts.

**restricted three-body** -- instead of simulating three free bodies, this one fixes the Sun and Jupiter on a circular orbit and asks: what does the gravitational landscape look like for a tiny spacecraft caught between them? the answer is a 3D potential surface with five special equilibrium points called Lagrange points (L1 through L5). L1, L2, L3 sit on the axis between them. L4 and L5 sit 60 degrees ahead and behind Jupiter and are actually stable, which is why the Trojan asteroids cluster there.

---

## how it works

### physics (C++)

written in C++ using the velocity Verlet integrator. better energy conservation than basic Euler integration.

each time step:
1. compute acceleration on each body from all others: `a = G * m / r^2`
2. update positions: `x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2`
3. recompute acceleration at new positions
4. update velocities: `v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt`

outputs trajectory data as CSV files, one per scenario.

for the restricted three-body scenario, a separate file computes the effective gravitational potential on a 300x300 grid in the co-rotating frame:

```
phi_eff = -G*M_sun/r1 - G*M_jup/r2 - 0.5*omega^2*(x^2+y^2)
```

the third term is the centrifugal contribution from the rotating frame. output is `potential.csv`.

### visualizer (Three.js)

one `index.html` file handles everything:

- glowing balls with radial gradient sprites
- fading trail lines for trajectory history
- starfield background
- drag to rotate, scroll to zoom
- adjustable playback speed
- restricted three-body view renders the potential as a colored 3D surface with L1-L5 marked directly on it

---

## project structure

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

## running locally

compile and run the C++ simulation:

```bash
g++ -O2 -o sim threebodyproblem.cpp
./sim
```

serve the visualizer:

```bash
python -m http.server 8000
```

then open `http://localhost:8000`.

---

## stack

- C++ for physics and potential computation
- Three.js for 3D WebGL rendering
- GitHub Pages for hosting

---

## references

- Chenciner & Montgomery (2000), figure-eight solution
- Lagrange (1772), equilateral triangle solution
- Euler (1767), collinear solution
