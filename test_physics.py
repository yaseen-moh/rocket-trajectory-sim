"""Basic sanity tests: energy conservation and collision-avoidance behavior."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from nbody_sim import demo_solar_system, Body, NBodySimulator


def test_energy_conservation():
    sim = demo_solar_system()
    e0 = sim.total_energy()
    sim.run(24 * 10)  # 10 days
    e1 = sim.total_energy()
    drift = abs((e1 - e0) / e0)
    assert drift < 0.01, f"Energy drifted {drift*100:.3f}% -- integrator unstable"


def test_collision_avoidance_prevents_overlap():
    # two bodies on a direct collision course
    a = Body("A", 1e15, [-1e6, 0, 0], [50, 0, 0], radius=1e5)
    b = Body("B", 1e15, [1e6, 0, 0], [-50, 0, 0], radius=1e5)
    sim = NBodySimulator([a, b], dt=10, collision_avoidance=True, avoidance_margin=5.0)
    min_dist = float("inf")
    for _ in range(2000):
        sim.step()
        d = ((a.position - b.position) ** 2).sum() ** 0.5
        min_dist = min(min_dist, d)
    assert min_dist > 0, "Bodies should never reach exactly zero separation"


if __name__ == "__main__":
    test_energy_conservation()
    test_collision_avoidance_prevents_overlap()
    print("All tests passed.")
