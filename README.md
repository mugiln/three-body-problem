# Gravitational Simulations

A real-time 3D gravitational simulation of the three-body problem, built from scratch in **C++** and visualized in the browser with **Three.js**.

🔴 **[Live Demo](https://mugiln.github.io/gravitational-simulations)** 🟢 🔵

---

## What is the Three-Body Problem?

The three-body problem asks: given three masses in space pulling on each other through gravity, how do they move over time?

Unlike the two-body problem (e.g. Earth orbiting the Sun), the three-body problem has **no general closed-form solution**. The motion is highly sensitive to initial conditions — a tiny change can produce completely different trajectories. This is one of the earliest known examples of **chaos theory**.

---

## Scenarios

### Figure-Eight
Three equal masses chasing each other in a perfect figure-eight. This exact solution was discovered by Chenciner & Montgomery in 2000. It is one of the few known stable periodic orbits of the three-body problem.

### Lagrange
Three equal masses sitting at the corners of an equilateral triangle, rotating around their common center of mass. This is Lagrange's 1772 solution — one of the five famous Lagrange points used in spacecraft positioning today (e.g. the James Webb Space Telescope sits at L2).

### Chaotic
A perturbed figure-eight — body 3 is nudged slightly out of the orbital plane. The system starts almost ordered, then slowly breaks down into unpredictable, chaotic motion. This demonstrates the **butterfly effect**: a tiny change in initial conditions leads to completely different behavior over time.

### Euler
Three bodies in a collinear configuration (all on a straight line), rotating around their common center of mass. This is Euler's 1767 solution. Unlike Lagrange, this orbit is **unstable** — small perturbations cause it to drift, which is visible in the overlapping trail pattern.

---

## How It Works

### Physics Engine (C++)
The simulation is written in C++ using the **Velocity Verlet** integration algorithm, which conserves energy far better than simple Euler integration.

For each time step:
1. Calculate gravitational acceleration on each body from all others using Newton's law: `a = G * m / r²`
2. Update positions: `x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt²`
3. Recalculate acceleration at new positions
4. Update velocities: `v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*dt`

The C++ program outputs trajectory data as CSV files, one per scenario.

### Visualizer (Three.js)
A single `index.html` file reads the CSV data and renders it in 3D using Three.js:
- Glowing balls with radial gradient sprites
- Fading trail lines showing recent trajectory history
- Starfield background
- Mouse drag to rotate, scroll to zoom
- Adjustable playback speed

---

## Project Structure

```
gravitational-simulations/
├── practice.cpp          # C++ simulation (physics engine)
├── index.html            # 3D visualizer (Three.js)
├── figure-eight.csv      # Trajectory data - figure-eight
├── lagrange.csv          # Trajectory data - Lagrange
├── chaotic.csv           # Trajectory data - chaotic
└── euler.csv             # Trajectory data - Euler
```

---

## Running Locally

### 1. Compile and run the C++ simulation
```bash
g++ -O2 -o simulation practice.cpp
./simulation
```
This generates the 4 CSV files.

### 2. Serve the visualizer locally
```bash
python -m http.server 8000
```
Then open `http://localhost:8000` in your browser.

---

## Tech Stack

- **C++** — physics simulation, Velocity Verlet integrator
- **Three.js** — 3D WebGL rendering
- **HTML/CSS/JS** — visualizer interface
- **GitHub Pages** — deployment

---

## References

- Chenciner, A. & Montgomery, R. (2000). 
- Lagrange, J. L. (1772). 
- Euler, L. (1767).
