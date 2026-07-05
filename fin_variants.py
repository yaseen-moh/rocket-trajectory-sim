"""
Compares apogee across a few fin-geometry variants -- mirrors the kind of
iteration done in OpenRocket when tuning body tube / fin geometry for
stability vs. drag trade-offs.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rocket_sim import RocketConfig, FinGeometry, DEMO_MOTOR_THRUST_CURVE, monte_carlo_apogee

base_kwargs = dict(
    name="variant", dry_mass=0.55, propellant_mass=0.06,
    body_diameter=0.041, length=0.9,
    cd_body=0.4, thrust_curve=DEMO_MOTOR_THRUST_CURVE, burn_time=2.0,
)

variants = {
    "small_fins": FinGeometry(root_chord=0.08, tip_chord=0.03, span=0.06, sweep_length=0.04),
    "baseline":   FinGeometry(root_chord=0.12, tip_chord=0.05, span=0.09, sweep_length=0.06),
    "large_fins": FinGeometry(root_chord=0.16, tip_chord=0.07, span=0.12, sweep_length=0.08),
}

if __name__ == "__main__":
    print(f"{'Variant':<12} {'Fin area (m²)':>14} {'Mean apogee (m)':>16} {'Std (m)':>10}")
    for name, fins in variants.items():
        rocket = RocketConfig(fins=fins, **base_kwargs)
        mc = monte_carlo_apogee(rocket, trials=60, base_wind_mps=3.0, gust_std=1.5)
        print(f"{name:<12} {fins.planform_area():>14.4f} {mc['mean_apogee_m']:>16.1f} {mc['std_apogee_m']:>10.1f}")
