import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button
import multiprocessing

def read_positions(file_path):
    positions = np.loadtxt(file_path, delimiter=',')
    dimensions = positions.shape[1] - 2
    positions_dict = {}
    for row in positions:
        cycle = int(row[0])
        particle_id = int(row[1])
        coordinates = row[2:]
        if cycle not in positions_dict:
            positions_dict[cycle] = []
        positions_dict[cycle].append((particle_id, *coordinates))
    return positions_dict, dimensions

def draw_shape(coordinates, dimension, ax):
    if dimension == 3:
        for coord in coordinates:
            particle_id, x, y, z = coord
            ax.scatter(x, y, z, c='r', marker='o')
            ax.text(x, y, z, str(particle_id), color='black', fontsize=8)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        lim = 2*100
        ax.set_xlim(-lim, lim)  # Set the X-axis limits
        ax.set_ylim(-lim, lim)  # Set the Y-axis limits
        ax.set_zlim(-lim, lim)  # Set the Z-axis limits
    elif dimension == 2:
        for coord in coordinates:
            particle_id, x, y = coord
            circle = plt.Circle((x, y), radius=0.1, color='r')
            ax.add_patch(circle)
            ax.text(x, y, str(particle_id), color='black', fontsize=8)
        ax.set_aspect('equal')
        ax.set_xlim(-10, 10)  # Set the X-axis limits
        ax.set_ylim(-10, 10)  # Set the Y-axis limits

def animate_positions(file_path):
    positions, dimension = read_positions(file_path)
    num_frames = len(positions)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d' if dimension == 3 else None)
    
    def update(frame):
        ax.clear()
        draw_shape(positions[frame], dimension, ax)
        ax.text2D(0.05, 0.95, f"Frame: {frame}", transform=ax.transAxes)
    
    def start_stop_animation(event):
        if animation.event_source is None:
            animation.event_source = fig.canvas.new_timer(interval=10)  # Decreased interval to 10 milliseconds
            animation.event_source.add_callback(animation._step)
            animation.event_source.start()
        else:
            animation.event_source.stop()
            animation.event_source = None
    
    start_stop_button_ax = plt.axes([0.8, 0.05, 0.1, 0.075])
    start_stop_button = Button(start_stop_button_ax, 'Start/Stop')
    start_stop_button.on_clicked(start_stop_animation)

    animation = FuncAnimation(fig, update, frames=num_frames, interval=20, repeat=False)
    plt.show()

# Example usage:
if __name__ == '__main__':
    pool = multiprocessing.Pool()
    pool.apply_async(animate_positions, ("../../NBody AMSC/positions.txt",))
    pool.close()
    pool.join()
