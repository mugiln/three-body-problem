import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

G = 1.0
dt = 0.004
total_time = 20.0
num_steps = int(total_time / dt)

# ─── General Three-Body: Scenarios ───────────────────────────────────────────

def initialize_bodies(scenario="figure-eight"):
    if scenario == "figure-eight":
        positions = np.array([
            [-0.97000436,  0.24308753, 0],
            [ 0.0,         0.0,        0],
            [ 0.97000436, -0.24308753, 0]
        ], dtype=float)
        velocities = np.array([
            [ 0.466203685,  0.43236573, 0],
            [-0.93240737,  -0.86473146, 0],
            [ 0.466203685,  0.43236573, 0]
        ], dtype=float)
        masses = np.array([1.0, 1.0, 1.0])

    elif scenario == "lagrange":
        R, m = 1.0, 1.0
        positions = np.array([
            [R, 0, 0],
            [R * np.cos(2*np.pi/3), R * np.sin(2*np.pi/3), 0],
            [R * np.cos(4*np.pi/3), R * np.sin(4*np.pi/3), 0]
        ])
        omega = np.sqrt(G * m / R**3)
        velocities = np.array([omega * np.array([-p[1], p[0], 0]) for p in positions])
        masses = np.array([m, m, m])

    elif scenario == "euler":
        m1, m2, m3 = 1.0, 2.0, 3.0
        masses = np.array([m1, m2, m3])
        positions = np.array([[-1.5, 0.0, 0.0], [0.0, 0.0, 0.0], [1.5, 0.0, 0.0]])
        total_mass = np.sum(masses)
        positions[:, 0] -= np.dot(masses, positions[:, 0]) / total_mass
        r = np.linalg.norm(positions, axis=1)
        omega = np.sqrt(G * total_mass / np.sum(masses * r**2))
        velocities = np.array([[0.0, omega * p[0], 0.0] for p in positions])

    elif scenario == "circular":
        positions = np.array([[1.0, 0.0, 0.0], [-0.5, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=float)
        velocities = np.array([[0.0, 0.7, 0.0], [0.0, -0.8, 0.0], [0.6, 0.0, 0.0]])
        masses = np.array([1.0, 1.5, 0.8])

    elif scenario == "chaotic":
        np.random.seed(42)
        positions = np.random.uniform(-1, 0.9, size=(3, 3))
        velocities = np.random.uniform(-0.5, 0.5, size=(3, 3))
        masses = np.array([1.0, 1.0, 1.0])
        positions -= np.sum(positions.T * masses, axis=1) / np.sum(masses)
        velocities -= np.sum(masses[:, np.newaxis] * velocities, axis=0) / np.sum(masses)

    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    return positions, velocities, masses


def calculate_acceleration(positions, masses):
    acc = np.zeros_like(positions)
    for i in range(len(masses)):
        for j in range(len(masses)):
            if i != j:
                r_ij = positions[j] - positions[i]
                r = np.linalg.norm(r_ij)
                if r > 1e-10:
                    acc[i] += G * masses[j] * r_ij / r**3
    return acc


def velocity_verlet_step(positions, velocities, masses, dt):
    acc = calculate_acceleration(positions, masses)
    pos_new = positions + velocities * dt + 0.5 * acc * dt**2
    acc_new = calculate_acceleration(pos_new, masses)
    vel_new = velocities + 0.5 * (acc + acc_new) * dt
    return pos_new, vel_new


def run_simulation(scenario="figure-eight"):
    positions, velocities, masses = initialize_bodies(scenario)
    num_bodies = len(masses)
    trajectory = np.zeros((num_steps, num_bodies, 3))
    trajectory[0] = positions.copy()
    for step in range(1, num_steps):
        positions, velocities = velocity_verlet_step(positions, velocities, masses, dt)
        trajectory[step] = positions.copy()
    return trajectory


def create_animation(trajectory, filename='three_body_animation.gif'):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    max_range = np.max(np.abs(trajectory))
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title('Three-Body Problem Trajectories')

    colors = ['r', 'g', 'b']
    lines = [ax.plot([], [], [], '-', lw=2, color=c)[0] for c in colors]

    def init():
        for line in lines:
            line.set_data([], [])
            line.set_3d_properties([])
        return lines

    def animate(i):
        frame = min(int((i / 100) * len(trajectory)), len(trajectory) - 1)
        for j, line in enumerate(lines):
            line.set_data(trajectory[:frame+1, j, 0], trajectory[:frame+1, j, 1])
            line.set_3d_properties(trajectory[:frame+1, j, 2])
        ax.set_title(f'Three-Body Trajectories ({(frame / len(trajectory)) * 100:.1f}%)')
        return lines

    anim = FuncAnimation(fig, animate, init_func=init, frames=101, interval=50, blit=True)
    print(f"Saving animation to {filename}...")
    anim.save(filename, writer='pillow', fps=20)
    print(f"Saved to {filename}")
    return anim


def compute_energy(positions, velocities, masses):
    KE = 0.5 * np.sum(masses[:, np.newaxis] * velocities**2)
    PE = sum(
        -G * masses[i] * masses[j] / np.linalg.norm(positions[i] - positions[j])
        for i in range(len(masses)) for j in range(i+1, len(masses))
    )
    return KE + PE


# ─── Restricted Three-Body (CR3BP) ───────────────────────────────────────────

def lagrange_points(mu):
    return {
        'L1': ((1 - mu) - (mu / (3 * (1 - mu)))**(1/3), 0),
        'L2': ((1 - mu) + (mu / (3 * (1 - mu)))**(1/3), 0),
        'L3': (-(1 + (5 * mu) / 12), 0),
        'L4': (0.5 - mu,  np.sqrt(3)/2),
        'L5': (0.5 - mu, -np.sqrt(3)/2),
    }


def cr3bp_acceleration(x, y, vx, vy, mu):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    ax = 2*vy + x - (1 - mu)*(x + mu)/r1**3 - mu*(x - 1 + mu)/r2**3
    ay = -2*vx + y - (1 - mu)*y/r1**3 - mu*y/r2**3
    return ax, ay


def effective_potential(x, y, mu):
    r1 = np.sqrt((x + mu)**2 + y**2)
    r2 = np.sqrt((x - 1 + mu)**2 + y**2)
    return 0.5*(x**2 + y**2) + (1 - mu)/r1 + mu/r2


def simulate_cr3bp(mu, point_name, perturbation=0.001, dt=0.001, T=500):
    N = int(T / dt)
    pts = lagrange_points(mu)
    x, y = pts[point_name]
    x += perturbation; y += perturbation
    vx = vy = 0.0
    ax, ay = cr3bp_acceleration(x, y, vx, vy, mu)

    trajectory = np.zeros((N, 2))
    trajectory[0] = [x, y]

    for i in range(1, N):
        x_new = x + vx*dt + 0.5*ax*dt**2
        y_new = y + vy*dt + 0.5*ay*dt**2
        ax_new, ay_new = cr3bp_acceleration(x_new, y_new, vx, vy, mu)
        vx += 0.5*(ax + ax_new)*dt
        vy += 0.5*(ay + ay_new)*dt
        x, y, ax, ay = x_new, y_new, ax_new, ay_new
        trajectory[i] = [x, y]

    return trajectory, pts


def plot_cr3bp(mu, point_name, trajectory, pts):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor('black')
    ax.plot(trajectory[:, 0], trajectory[:, 1], 'white', lw=1.5, label=f'Trajectory near {point_name}')
    ax.plot(-mu, 0, 'bo', label='Primary m1')
    ax.plot(1 - mu, 0, 'ro', label='Primary m2')
    for name, coord in pts.items():
        ax.plot(*coord, 'y*', markersize=10)
        ax.text(coord[0] + 0.02, coord[1] + 0.02, name, color='white', fontsize=9)
    ax.set_title(f'CR3BP Trajectory near {point_name}', color='white')
    ax.set_xlabel('x', color='white'); ax.set_ylabel('y', color='white')
    ax.tick_params(colors='white')
    ax.set_aspect('equal')
    ax.legend(facecolor='black', edgecolor='white', labelcolor='white')
    plt.show()


def plot_effective_potential(mu):
    x_vals = np.linspace(-1.5, 1.5, 500)
    y_vals = np.linspace(-1.5, 1.5, 500)
    X, Y = np.meshgrid(x_vals, y_vals)
    with np.errstate(divide='ignore', invalid='ignore'):
        Z = effective_potential(X, Y, mu)
        Z[np.isnan(Z)] = np.inf
        Z[Z > 10] = 10

    pts = lagrange_points(mu)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_facecolor('black')
    contour = ax.contour(X, Y, Z, levels=30, cmap='viridis')
    ax.plot(-mu, 0, 'bo', label='Primary m1')
    ax.plot(1 - mu, 0, 'ro', label='Primary m2')
    for name, coord in pts.items():
        ax.plot(*coord, 'y*', markersize=10)
        ax.text(coord[0] + 0.02, coord[1] + 0.02, name, color='white', fontsize=9)
    ax.set_title('Effective Potential (Ω) in the Rotating Frame', color='white')
    ax.set_xlabel('x', color='white'); ax.set_ylabel('y', color='white')
    ax.tick_params(colors='white')
    ax.set_aspect('equal')
    fig.colorbar(contour, ax=ax, label='Ω', shrink=0.8)
    plt.show()


# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # General three-body simulation
    scenario = "chaotic"
    time_map = {"figure-eight": 20.0, "lagrange": 200.0, "euler": 50.0}
    total_time = time_map.get(scenario, 10.0)
    num_steps = int(total_time / dt)

    print(f"Running general three-body: {scenario}")
    trajectory = run_simulation(scenario)
    print(f"Simulation done: {len(trajectory)} steps")

    animation = create_animation(trajectory)

    positions, velocities, masses = initialize_bodies(scenario)
    E0 = compute_energy(positions, velocities, masses)

    for _ in range(num_steps):
        positions, velocities = velocity_verlet_step(positions, velocities, masses, dt)
    Ef = compute_energy(positions, velocities, masses)

    print(f"Initial energy: {E0:.6f}")
    print(f"Final energy:   {Ef:.6f}")
    print(f"Energy drift:   {(Ef - E0)/E0 * 100:.6f}%")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    max_range = np.max(np.abs(trajectory))
    ax.set_xlim(-max_range, max_range)
    ax.set_ylim(-max_range, max_range)
    ax.set_zlim(-max_range, max_range)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.set_title(f'Three-Body Problem — {scenario.title()}')
    colors = ['r', 'g', 'b']
    for i in range(3):
        ax.plot(trajectory[:, i, 0], trajectory[:, i, 1], trajectory[:, i, 2],
                '-', color=colors[i], lw=2, label=f'Body {i+1}')
        ax.plot(trajectory[0, i, 0], trajectory[0, i, 1], trajectory[0, i, 2],
                'o', color=colors[i], markersize=8)
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"three_body_{scenario}_trajectories.png", dpi=300)
    plt.show()

    # Restricted three-body (CR3BP)
    mu = 0.01
    selected_point = 'L4'

    print(f"\nRunning CR3BP near {selected_point}")
    traj, pts = simulate_cr3bp(mu, selected_point)
    plot_cr3bp(mu, selected_point, traj, pts)
    plot_effective_potential(mu)
