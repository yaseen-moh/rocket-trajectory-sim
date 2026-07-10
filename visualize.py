"""
Real-time (animated) visualization of the N-body simulator.
Run: python visualize.py
Produces an .mp4/.gif if ffmpeg/imagemagick is available, otherwise
shows a live matplotlib animation window.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from nbody_sim import demo_solar_system

sim = demo_solar_system()
n_frames = 400
sub_steps_per_frame = 6  # advance sim faster than the render rate

fig, ax = plt.subplots(figsize=(7, 7))
ax.set_facecolor("black")
fig.patch.set_facecolor("black")
ax.set_aspect("equal")
lim = 2.2e11
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xticks([])
ax.set_yticks([])

colors = ["gold", "royalblue", "silver", "orangered"]
points = [ax.plot([], [], "o", color=colors[i % len(colors)], markersize=8)[0]
          for i in range(len(sim.bodies))]
lines = [ax.plot([], [], "-", color=colors[i % len(colors)], linewidth=0.6, alpha=0.5)[0]
         for i in range(len(sim.bodies))]
title = ax.text(0.02, 0.95, "", color="white", transform=ax.transAxes)


def init():
    for p, l in zip(points, lines):
        p.set_data([], [])
        l.set_data([], [])
    return points + lines + [title]


def animate(frame):
    sim.run(sub_steps_per_frame)
    for b, p, l in zip(sim.bodies, points, lines):
        p.set_data([b.position[0]], [b.position[1]])
        trail = np.array(b.trail[-200:])
        l.set_data(trail[:, 0], trail[:, 1])
    title.set_text(f"t = {sim.time / 86400:.1f} days | collisions logged: "
                    f"{len(sim.collision_events)}")
    return points + lines + [title]


anim = FuncAnimation(fig, animate, init_func=init, frames=n_frames,
                      interval=30, blit=True)

if __name__ == "__main__":
    try:
        anim.save("nbody_demo.gif", writer="pillow", fps=30)
        print("Saved animation to nbody_demo.gif")
    except Exception as e:
        print(f"Could not save gif ({e}); showing live window instead.")
        plt.show()
