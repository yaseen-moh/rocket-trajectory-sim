"""
N-Body Orbital Mechanics Simulator with Real-Time Collision Avoidance
-----------------------------------------------------------------------
Simulates gravitational interactions between an arbitrary number of bodies
using direct pairwise force summation and a velocity-Verlet (leapfrog)
integrator, which conserves energy far better than naive Euler integration
over long simulation runs.

Collision avoidance: each body is treated as a sphere of finite radius.
Before every integration step, the simulator looks ahead along each body's
current velocity vector to see whether it will pass within a critical
distance of another body during the step. If so, it applies a soft
repulsive delta-v (a scaled normal impulse) so bodies deflect around each
other instead of clipping through / merging, similar to a simplified
proximity-avoidance maneuver.

Author: Yaseen Mohamed
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field


G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2


@dataclass
class Body:
    name: str
    mass: float                 # kg
    position: np.ndarray         # shape (3,), meters
    velocity: np.ndarray         # shape (3,), m/s
    radius: float = 1.0e5        # meters, used for collision geometry
    trail: list = field(default_factory=list)

    def __post_init__(self):
        self.position = np.array(self.position, dtype=float)
        self.velocity = np.array(self.velocity, dtype=float)


class NBodySimulator:
    def __init__(self, bodies: list[Body], dt: float = 60.0,
                 collision_avoidance: bool = True,
                 avoidance_margin: float = 3.0,
                 avoidance_strength: float = 0.15):
        """
        Parameters
        ----------
        bodies : list of Body
        dt : float
            Integration timestep in seconds.
        collision_avoidance : bool
            If True, apply repulsive delta-v when bodies are predicted to
            come within `avoidance_margin` * (sum of radii) of each other.
        avoidance_margin : float
            Multiplier on the combined body radii that defines the
            "danger zone" that triggers an avoidance maneuver.
        avoidance_strength : float
            Scales how strong the corrective delta-v is (fraction of
            relative velocity redirected per step while inside the zone).
        """
        self.bodies = bodies
        self.dt = dt
        self.collision_avoidance = collision_avoidance
        self.avoidance_margin = avoidance_margin
        self.avoidance_strength = avoidance_strength
        self.time = 0.0
        self.collision_events = []  # log of (time, body_a, body_b, distance)

    # ------------------------------------------------------------------
    # Core physics
    # ------------------------------------------------------------------
    def _pairwise_accelerations(self, positions: np.ndarray) -> np.ndarray:
        """Vectorized O(n^2) gravitational acceleration on every body."""
        n = len(self.bodies)
        acc = np.zeros((n, 3))
        for i in range(n):
            diff = positions - positions[i]          # (n,3)
            dist2 = np.einsum('ij,ij->i', diff, diff)
            dist2[i] = np.inf                          # avoid self-interaction
            dist3 = np.power(dist2, 1.5)
            masses = np.array([b.mass for b in self.bodies])
            acc[i] = np.sum((G * masses[:, None] * diff) / dist3[:, None], axis=0)
        return acc

    def _apply_collision_avoidance(self):
        """Look-ahead check + soft repulsive delta-v between close bodies."""
        n = len(self.bodies)
        for i in range(n):
            for j in range(i + 1, n):
                bi, bj = self.bodies[i], self.bodies[j]
                rel_pos = bj.position - bi.position
                dist = np.linalg.norm(rel_pos)
                danger_dist = self.avoidance_margin * (bi.radius + bj.radius)

                if dist < danger_dist and dist > 1e-6:
                    normal = rel_pos / dist
                    rel_vel = bj.velocity - bi.velocity
                    closing_speed = -np.dot(rel_vel, normal)

                    if closing_speed > 0:  # bodies approaching each other
                        # push apart along the line of centers, scaled by
                        # how deep into the danger zone they are
                        penetration = 1.0 - (dist / danger_dist)
                        impulse = self.avoidance_strength * closing_speed * penetration
                        bi.velocity -= normal * impulse * 0.5
                        bj.velocity += normal * impulse * 0.5

                if dist < (bi.radius + bj.radius):
                    self.collision_events.append((self.time, bi.name, bj.name, dist))

    def step(self):
        """Advance the simulation by one velocity-Verlet timestep."""
        positions = np.array([b.position for b in self.bodies])
        velocities = np.array([b.velocity for b in self.bodies])

        acc_old = self._pairwise_accelerations(positions)

        # kick-drift
        new_positions = positions + velocities * self.dt + 0.5 * acc_old * self.dt ** 2
        for b, p in zip(self.bodies, new_positions):
            b.position = p

        acc_new = self._pairwise_accelerations(new_positions)

        new_velocities = velocities + 0.5 * (acc_old + acc_new) * self.dt
        for b, v in zip(self.bodies, new_velocities):
            b.velocity = v

        if self.collision_avoidance:
            self._apply_collision_avoidance()

        for b in self.bodies:
            b.trail.append(b.position.copy())

        self.time += self.dt

    def run(self, steps: int):
        for _ in range(steps):
            self.step()

    # ------------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------------
    def total_energy(self) -> float:
        """Kinetic + potential energy of the system (should stay ~constant)."""
        ke = sum(0.5 * b.mass * np.dot(b.velocity, b.velocity) for b in self.bodies)
        pe = 0.0
        n = len(self.bodies)
        for i in range(n):
            for j in range(i + 1, n):
                dist = np.linalg.norm(self.bodies[i].position - self.bodies[j].position)
                pe -= G * self.bodies[i].mass * self.bodies[j].mass / dist
        return ke + pe


def demo_solar_system() -> NBodySimulator:
    """Rough Sun-Earth-Moon-ish 3-body system for a quick demo run."""
    sun = Body("Sun", 1.989e30, [0, 0, 0], [0, 0, 0], radius=6.96e8)
    earth = Body("Earth", 5.972e24, [1.496e11, 0, 0], [0, 29780, 0], radius=6.371e6)
    moon = Body("Moon", 7.348e22, [1.496e11 + 3.844e8, 0, 0],
                [0, 29780 + 1022, 0], radius=1.737e6)
    rogue = Body("Rogue Asteroid", 1.0e20,
                 [1.6e11, -2.0e10, 0], [-8000, 6000, 0], radius=5.0e5)
    return NBodySimulator([sun, earth, moon, rogue], dt=3600, avoidance_margin=4.0)


if __name__ == "__main__":
    sim = demo_solar_system()
    e0 = sim.total_energy()
    sim.run(24 * 30)  # simulate 30 days at 1-hour steps
    e1 = sim.total_energy()
    print(f"Initial energy: {e0:.4e} J")
    print(f"Final energy:   {e1:.4e} J")
    print(f"Relative drift: {abs((e1 - e0) / e0) * 100:.6f}%")
    print(f"Collision/near-miss events logged: {len(sim.collision_events)}")
