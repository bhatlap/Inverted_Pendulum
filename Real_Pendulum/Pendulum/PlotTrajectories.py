import numpy as np
import matplotlib.pyplot as plt


class PlotTrajectories:
    def __init__(self, dt, stopIter, title=" "):
        self.fig = plt.figure(num=None, figsize=(5, 10), dpi=110) # Create a figure for plotting
        self.fig.suptitle(title) # Title
        self.axes = [0, 0, 0, 0, 0] # Initialization of axis
        self.line = [0, 0, 0, 0, 0] # and line
        self.line2 = [0, 0, 0, 0, 0]

        label1 = ['$\Theta$ sim', '$\dot{\Theta}$ sim', '$x$', '$\dot{x} sim$', '$u$']
        label2 = ['$\Theta$ real', '$\dot{\Theta}$ real', '$x$ real', '$\dot{x}$ real', '$u$']

        for i in range(5):
            self.axes[i] = self.fig.add_subplot(5, 1, i + 1)
            self.line[i], = self.axes[i].plot([], [], label=label1[i])
            self.line2[i], = self.axes[i].plot([], [], label=label2[i])
            self.line2[i].set_color("red")

        self.fig.tight_layout(pad=2) # Makes axis fit to data

        self.axes[0].set_ylabel('$\Theta [rad]$')
        self.axes[1].set_ylabel('$\dot{\Theta} [rad/s]$')
        self.axes[2].set_ylabel('$ p[m]$')
        self.axes[3].set_ylabel('$\dot{p} [m/s]$')
        self.axes[4].set_ylabel('$u [m/s^2]$')

        if stopIter is None:
            stopIter = 10000

        # Setting x axis parameters
        for ax in self.axes:
            ax.set_xlabel("$t$ [s]")
            ax.grid()
            ax.legend()
            ax.set_xlim([0, dt*stopIter])

        # Axis limits
        self.axes[0].set_ylim([-6.5, 6.5])
        self.axes[1].set_ylim([-23, 23])
        self.axes[2].set_ylim([-0.5, 0.5])
        self.axes[3].set_ylim([-10, 10])
        self.axes[4].set_ylim([-23, 23])

        plt.show(block=False)
        self.fig.canvas.draw()

        # Lines for each state and input
        for i in range(5):
            self.axes[i].draw_artist(self.line[i])
            self.axes[i].draw_artist(self.line2[i])

        # draw the animated artist, this uses a cached renderer
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)

        # show the result to the screen, this pushes the updated RGBA buffer from the
        # renderer to the GUI framework so you can see it
        self.fig.canvas.blit(self.fig.bbox)

    def update_plot(self, x_trajectory, u_trajectory, u_trajectory2,dt, extra_trajectory=" "):
        x_trajectory_conc = np.concatenate(x_trajectory, axis=1)
        u_trajectory_conc = np.array(u_trajectory)
        u_trajectory_conc2 = np.array(u_trajectory2)
        time = np.arange(0, x_trajectory_conc.shape[1]) * dt

        for i in range(4):
            self.line[i].set_data(time, x_trajectory_conc[i,:])

        self.line[4].set_data(time[:-1], u_trajectory_conc)
        self.line2[4].set_data(time[:-1], u_trajectory_conc2)

        self.fig.canvas.restore_region(self.bg)

        for i in range(5):

            self.axes[i].draw_artist(self.line[i])
            self.axes[i].draw_artist(self.line2[i])

        self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()
