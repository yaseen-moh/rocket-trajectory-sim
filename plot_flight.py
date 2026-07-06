"""Generates altitude/velocity/range plots for a single simulated flight."""
import matplotlib.pyplot as plt
from rocket_sim import demo_rocket, simulate_flight

rocket = demo_rocket()
result = simulate_flight(rocket, launch_angle_deg=5, base_wind_mps=3.0, gust_std=1.5, seed=42)

fig, axes = plt.subplots(1, 3, figsize=(14, 4))

axes[0].plot(result["t"], result["y"])
axes[0].set_title("Altitude vs Time")
axes[0].set_xlabel("time (s)")
axes[0].set_ylabel("altitude (m)")

axes[1].plot(result["x"], result["y"])
axes[1].set_title("Flight Path (side view)")
axes[1].set_xlabel("downrange (m)")
axes[1].set_ylabel("altitude (m)")
axes[1].set_aspect("equal")

speed = (result["vx"]**2 + result["vy"]**2) ** 0.5
axes[2].plot(result["t"], speed)
axes[2].set_title("Speed vs Time")
axes[2].set_xlabel("time (s)")
axes[2].set_ylabel("speed (m/s)")

plt.tight_layout()
plt.savefig("flight_profile.png", dpi=150)
print("Saved flight_profile.png")
print(f"Apogee: {result['apogee_m']:.1f} m | Range: {result['range_m']:.1f} m")
