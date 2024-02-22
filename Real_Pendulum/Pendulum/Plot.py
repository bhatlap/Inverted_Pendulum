import numpy as np
import matplotlib.pyplot as plt


class Plot:
    def __init__(self, dt, stopIter, simulation, title="Online MPC Plotting"):

        self.simulation = simulation

        if not self.simulation:
            self.fig = plt.figure(num=None, figsize=(5, 10), dpi=110)  # Create a figure for plotting
            self.fig.suptitle(title)  # Title
            self.axes = [0, 0, 0, 0, 0, 0]  # Initialization of axis
            self.line_calc = [0, 0, 0, 0, 0, 0]  # and line
            self.line_meas = [0, 0, 0, 0, 0, 0]

            label_meas= ['$\Theta$ measured', '$\dot{\Theta}$ measured ', '$x$ measured',
                      '$\dot{x}$ measured', '$u$ sent','']

            label_calc = ['$\Theta$ calculated ', '$\dot{\Theta}$ calculated ', '$x$ calculated ',
                      '$\dot{x}$ calculated', '$u$ received','']

            self.axes[0] = self.fig.add_subplot(6, 1, 0 + 1)
            self.line_calc[0], = self.axes[0].plot([], [], label=label_meas[0])
            self.line_meas[0], = self.axes[0].plot([], [], label=label_calc[0])
            self.line_meas[0].set_color("red")
            self.line_meas[0].set_linestyle('--')

            for i in range(5):
                self.axes[i+1] = self.fig.add_subplot(6, 1, i + 2, sharex=self.axes[0])
                self.line_calc[i+1], = self.axes[i+1].plot([], [], label=label_meas[i+1])
                self.line_meas[i+1], = self.axes[i+1].plot([], [], label=label_calc[i+1])
                self.line_meas[i+1].set_color("red")
                self.line_meas[i+1].set_linestyle('--')

            self.fig.tight_layout(pad=2)  # Makes axis fit to data

            self.axes[0].set_ylabel('$\Theta [rad]$')
            self.axes[1].set_ylabel('$\dot{\Theta} [rad/s]$')
            self.axes[2].set_ylabel('$ p[m]$')
            self.axes[3].set_ylabel('$\dot{p} [m/s]$')
            self.axes[4].set_ylabel('$u [m/s^2]$')
            self.axes[5].set_ylabel('Solver Status')

        else:
            self.fig = plt.figure(num=None, figsize=(5, 10), dpi=110)  # Create a figure for plotting
            self.fig.suptitle(title)  # Title
            self.axes = [0, 0, 0, 0, 0, 0]  # Initialization of axis
            self.line_calc = [0, 0, 0, 0, 0, 0]  # and line

            label_meas = ['$\Theta$ ', '$\dot{\Theta}$ ', '$x$ ',
                          '$\dot{x}$', '$u$', '']

            self.axes[0] = self.fig.add_subplot(6, 1, 0 + 1)
            self.line_calc[0], = self.axes[0].plot([], [], label=label_meas[0])

            for i in range(5):
                self.axes[i + 1] = self.fig.add_subplot(6, 1, i + 2, sharex=self.axes[0])
                self.line_calc[i + 1], = self.axes[i + 1].plot([], [], label=label_meas[i + 1])

            self.fig.tight_layout(pad=2)  # Makes axis fit to data

            self.axes[0].set_ylabel('$\Theta [rad]$')
            self.axes[1].set_ylabel('$\dot{\Theta} [rad/s]$')
            self.axes[2].set_ylabel('$ p[m]$')
            self.axes[3].set_ylabel('$\dot{p} [m/s]$')
            self.axes[4].set_ylabel('$u [m/s^2]$')
            self.axes[5].set_ylabel('Solver Status')

        if stopIter is None:
            stopIter = 5000

        # Setting x axis parameters
        for ax in self.axes:
            ax.set_xlabel("$t$ [s]")
            ax.grid()
            ax.legend(loc='upper right')
            ax.set_xlim([0, dt*stopIter])

        # Axis limits
        self.axes[0].set_ylim([-6.5, 6.5])
        self.axes[1].set_ylim([-30, 30])
        self.axes[2].set_ylim([-0.5, 0.5])
        self.axes[3].set_ylim([-2.5, 2.5])
        self.axes[4].set_ylim([-25, 25])
        self.axes[5].set_ylim([-1, 5])

        plt.show(block=False)
        self.fig.canvas.draw()

        # Lines for each state and input
        if not self.simulation:
            for i in range(6):
                self.axes[i].draw_artist(self.line_calc[i])
                self.axes[i].draw_artist(self.line_meas[i])
        else:
            for i in range(6):
                self.axes[i].draw_artist(self.line_calc[i])

        # draw the animated artist, this uses a cached renderer
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)

        # show the result to the screen, this pushes the updated RGBA buffer from the
        # renderer to the GUI framework so you can see it
        self.fig.canvas.blit(self.fig.bbox)

    def update_plot(self, dt, x_traj_calc, x_traj_meas, u_traj_calc, u_traj_meas, sixth_axe):

        if not self.simulation:
            x_meas_conc = np.concatenate(x_traj_meas, axis=1)
            x_calc_conc = np.concatenate(x_traj_calc, axis=1)
            u_calc_conc = np.array(u_traj_calc)
            u_meas_conc = np.array(u_traj_meas)
            sixth_axe_conc = np.array(sixth_axe)
            time = np.arange(0, x_meas_conc.shape[1]) * dt

            for i in range(4):
                self.line_meas[i].set_data(time, x_meas_conc[i, :])
                self.line_calc[i].set_data(time, x_calc_conc[i, :])

            self.line_meas[4].set_data(time[:-1], u_meas_conc)
            self.line_calc[4].set_data(time[:-1], u_calc_conc)
            self.line_calc[5].set_data(time[:-1], sixth_axe_conc)

            self.fig.canvas.restore_region(self.bg)

            for i in range(6):

                self.axes[i].draw_artist(self.line_meas[i])
                self.axes[i].draw_artist(self.line_calc[i])

        else:
            x_calc_conc = np.concatenate(x_traj_calc, axis=1)
            u_calc_conc = np.array(u_traj_calc)
            sixth_axe_conc = np.array(sixth_axe)
            time = np.arange(0, x_calc_conc.shape[1]) * dt

            for i in range(4):
                self.line_calc[i].set_data(time, x_calc_conc[i, :])

            self.line_calc[4].set_data(time[:-1], u_calc_conc)
            self.line_calc[5].set_data(time[:-1], sixth_axe_conc)

            self.fig.canvas.restore_region(self.bg)

            for i in range(6):
                self.axes[i].draw_artist(self.line_calc[i])

        self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()


