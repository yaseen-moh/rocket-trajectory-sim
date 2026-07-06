"""
Rocket Flight Trajectory Simulator
-----------------------------------
A 3-degree-of-freedom (point-mass, planar) rocket trajectory model inspired
by the workflow used with OpenRocket / Ansys STK for the American Rocketry
Challenge: define body-tube diameter, fin geometry, and a motor thrust
curve, then integrate the trajectory under variable wind and atmospheric
conditions to predict apogee, stability, and landing footprint.

This does NOT reimplement OpenRocket's full 6-DOF finite-element model —
it's a lighter-weight point-mass simulator that is fast enough to run
hundreds of Monte Carlo wind-variation trials, which is useful for a
"how sensitive is my apogee to gusty conditions" style analysis.

Author: Yaseen Mohamed
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_ivp

G0 = 9.80665           # m/s^2, standard gravity
R_AIR = 287.05         # J/(kg K)
GAMMA = 1.4


@dataclass
class FinGeometry:
    root_chord: float      # m
    tip_chord: float       # m
    span: float            # m
    sweep_length: float    # m
    count: int = 4

    def planform_area(self) -> float:
        return 0.5 * (self.root_chord + self.tip_chord) * self.span * self.count


@dataclass
class RocketConfig:
    name: str
    dry_mass: float             # kg, mass with no propellant
    propellant_mass: float      # kg
    body_diameter: float        # m
    length: float                # m
    fins: FinGeometry
    cd_body: float = 0.45        # baseline drag coefficient of body alone
    thrust_curve: list = None     # list of (time_s, thrust_N) tuples
    burn_time: float = 2.0

    def reference_area(self) -> float:
        return np.pi * (self.body_diameter / 2) ** 2

    def total_drag_coefficient(self) -> float:
        # crude but directionally correct: extra fin area increases parasitic
        # drag roughly in proportion to its fraction of the reference area
        fin_drag_fraction = self.fins.planform_area() / self.reference_area()
        return self.cd_body + 0.15 * fin_drag_fraction

    def mass_at(self, t: float) -> float:
        if t >= self.burn_time:
            return self.dry_mass
        # linear propellant burn assumption
        return self.dry_mass + self.propellant_mass * (1 - t / self.burn_time)

    def thrust_at(self, t: float) -> float:
        if not self.thrust_curve or t > self.burn_time:
            return 0.0
        times = np.array([p[0] for p in self.thrust_curve])
        thrusts = np.array([p[1] for p in self.thrust_curve])
        return float(np.interp(t, times, thrusts, left=0.0, right=0.0))


def isa_density(altitude_m: float) -> float:
    """International Standard Atmosphere density approximation (troposphere)."""
    T0, P0, L, rho0 = 288.15, 101325.0, 0.0065, 1.225
    if altitude_m < 11000:
        T = T0 - L * altitude_m
        P = P0 * (T / T0) ** (G0 / (R_AIR * L))
    else:
        T = 216.65
        P = 22632.06 * np.exp(-G0 * (altitude_m - 11000) / (R_AIR * T))
    return P / (R_AIR * T)


def wind_vector(altitude_m: float, base_wind_mps: float, gust_std: float,
                 rng: np.random.Generator) -> float:
    """Simple altitude-independent horizontal wind + Gaussian gust noise."""
    return base_wind_mps + rng.normal(0, gust_std)


def simulate_flight(rocket: RocketConfig, launch_angle_deg: float = 5.0,
                     base_wind_mps: float = 3.0, gust_std: float = 1.0,
                     seed: int | None = None, t_max: float = 60.0):
    """
    Integrates planar (x = downrange, y = altitude) trajectory.
    State vector: [x, y, vx, vy]
    Returns a dict with time series and summary stats (apogee, range, etc.)
    """
    rng = np.random.default_rng(seed)
    theta0 = np.radians(90 - launch_angle_deg)  # angle from horizontal

    def rhs(t, state):
        x, y, vx, vy = state
        y = max(y, 0.0)
        m = rocket.mass_at(t)
        thrust = rocket.thrust_at(t)

        speed = np.hypot(vx, vy)
        if speed > 1e-6:
            dir_x, dir_y = vx / speed, vy / speed
        else:
            dir_x, dir_y = np.cos(theta0), np.sin(theta0)

        rho = isa_density(y)
        cd = rocket.total_drag_coefficient()
        area = rocket.reference_area()
        wind = wind_vector(y, base_wind_mps, gust_std, rng)

        rel_vx = vx - wind
        rel_speed = np.hypot(rel_vx, vy)
        drag = 0.5 * rho * cd * area * rel_speed ** 2
        drag_x = -drag * (rel_vx / rel_speed) if rel_speed > 1e-6 else 0.0
        drag_y = -drag * (vy / rel_speed) if rel_speed > 1e-6 else 0.0

        ax = (thrust * dir_x + drag_x) / m
        ay = (thrust * dir_y + drag_y) / m - G0

        # Launch rail / pad constraint: the rocket sits still until thrust
        # overcomes weight, and can't be dragged below ground level once
        # it has landed.
        if y <= 0.0:
            if thrust <= m * G0:
                return [0.0, 0.0, 0.0, 0.0]
            vx, vy = 0.0, max(vy, 0.0)

        return [vx, vy, ax, ay]

    v0x = 0.0
    v0y = 0.0
    state0 = [0.0, 0.0, v0x, v0y]

    # Integrate the full time window, then trim to the first landing point
    # after apogee in post-processing. This avoids brittle event-detection
    # edge cases right at liftoff (altitude starts at ~0 while thrust ramps
    # up from zero), while still giving a clean "touchdown" trajectory.
    sol = solve_ivp(rhs, [0, t_max], state0, max_step=0.02, dense_output=True)

    t = sol.t
    x, y, vx, vy = sol.y
    y = np.clip(y, 0.0, None)
    apogee_idx = int(np.argmax(y))

    landing_idx = len(t) - 1
    for idx in range(apogee_idx + 1, len(t)):
        if y[idx] <= 0.0:
            landing_idx = idx
            break

    t = t[:landing_idx + 1]
    x = x[:landing_idx + 1]
    y = y[:landing_idx + 1]
    vx = vx[:landing_idx + 1]
    vy = vy[:landing_idx + 1]

    return {
        "t": t, "x": x, "y": y, "vx": vx, "vy": vy,
        "apogee_m": y[apogee_idx],
        "apogee_time_s": t[apogee_idx],
        "range_m": x[-1],
        "flight_time_s": t[-1],
        "max_velocity_mps": float(np.max(np.hypot(vx, vy))),
    }


def monte_carlo_apogee(rocket: RocketConfig, trials: int = 200, **kwargs) -> dict:
    """Runs repeated flights with randomized gusts to assess apogee sensitivity."""
    apogees = []
    for i in range(trials):
        result = simulate_flight(rocket, seed=i, **kwargs)
        apogees.append(result["apogee_m"])
    apogees = np.array(apogees)
    return {
        "mean_apogee_m": apogees.mean(),
        "std_apogee_m": apogees.std(),
        "min_apogee_m": apogees.min(),
        "max_apogee_m": apogees.max(),
        "all_apogees": apogees,
    }


DEMO_MOTOR_THRUST_CURVE = [
    (0.0, 0.0), (0.05, 180.0), (0.2, 220.0), (0.6, 200.0),
    (1.2, 160.0), (1.8, 90.0), (2.0, 0.0),
]


def demo_rocket() -> RocketConfig:
    fins = FinGeometry(root_chord=0.12, tip_chord=0.05, span=0.09,
                        sweep_length=0.06, count=4)
    return RocketConfig(
        name="ARC Demo Rocket",
        dry_mass=0.55, propellant_mass=0.06,
        body_diameter=0.041, length=0.9,
        fins=fins, cd_body=0.4,
        thrust_curve=DEMO_MOTOR_THRUST_CURVE, burn_time=2.0,
    )


if __name__ == "__main__":
    rocket = demo_rocket()
    result = simulate_flight(rocket, launch_angle_deg=5, base_wind_mps=3.0, gust_std=1.5, seed=42)
    print(f"Apogee: {result['apogee_m']:.1f} m at t={result['apogee_time_s']:.2f}s")
    print(f"Downrange distance: {result['range_m']:.1f} m")
    print(f"Total flight time: {result['flight_time_s']:.2f} s")
    print(f"Max velocity: {result['max_velocity_mps']:.1f} m/s")

    print("\nRunning Monte Carlo wind sensitivity analysis (200 trials)...")
    mc = monte_carlo_apogee(rocket, trials=200, base_wind_mps=3.0, gust_std=1.5)
    print(f"Mean apogee: {mc['mean_apogee_m']:.1f} m  (std: {mc['std_apogee_m']:.1f} m)")
    print(f"Range: {mc['min_apogee_m']:.1f} m to {mc['max_apogee_m']:.1f} m")
