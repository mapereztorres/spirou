import matplotlib.pyplot as plt
import numpy as np

class Plotter:
    def __init__(self):
        self.lw = 3
        self.multiple = 1
        self.x = None
        self.y_min = None
        self.y_max = None
        self.y_min_ZL = None
        self.y_max_ZL = None

    def plot(self):
        self.setup_plot()
        self.plot_data()
        self.set_axes_limits()
        self.draw_vertical_line()
        self.draw_rms_line()
        self.draw_earth()
        self.print_parameters()
        self.save_plot()

    def setup_plot(self):
        if self.multiple:
            plt.figure(figsize=(8, 11))
            ax1 = plt.subplot2grid((3, 1), (0, 0), rowspan=1, colspan=1)
            ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan=2, colspan=1)
        else:
            plt.figure(figsize=(8, 7.5))
            ax2 = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

    def plot_data(self):
        ax1.plot(self.x, np.log10(M_A), color='k', lw=self.lw)
        ax2.fill_between(self.x, self.y_min_ZL, self.y_max_ZL, color="blue", alpha=0.7, label="Zarka/Lanza model")
        ax2.fill_between(self.x, self.y_min, self.y_max, color="orange", alpha=0.7, label="Saur/Turnpenney model")

    def set_axes_limits(self):
        xmin = np.amin(d_orb) / R_star / R_star
        xmax = np.amax(d_orb) / R_star / R_star
        ax11.set_xlim([xmin, xmax])
        ax1.set_xlim([xmin, xmax])
        ax2.set_xlim([xmin, xmax])
        ax1.set_ylim([-3, 0])
        ax2.set_ylim([-3, 3])

    def draw_vertical_line(self):
        ax1.axvline(x=r_orb / R_star, ls='--', color='k', lw=2)
        ax2.axvline(x=r_orb / R_star, ls='--', color='k', lw=2)
       
    def draw_rms_line(self):
        ax2.axhline(y=np.log10(3 * rms), ls='-.', color='grey', lw=2)
        xpos = d_orb_max / 6
        ypos = np.log10(4 * rms)
        ax2.text(x=xpos, y=ypos, s=r"3$\times$RMS", fontsize='small')

    def draw_earth(self):
        if draw_earth:
            paths = ['./pics/earth.png']
            x = [r_orb / R_star]
            y = [np.log10(3 * rms)]
            for x0, y0, path in zip(x, y, paths):
                ab_earth = AnnotationBbox(spi.getImage(path), (x0, y0), frameon=False)
                ax2.add_artist(ab_earth)

    def print_parameters(self):
        d_diff = np.abs((d_orb - r_orb) / R_star)
        loc_pl = np.where(d_diff == d_diff.min())
        B_pl_loc = round(float(Bp[loc_pl] / (bfield_earth * Tesla2Gauss)), 2)
        xpos = xmax * 0.8
        ypos_offset = (ymax - ymin) / 8
        ypos = (ymax - ymin) / 4 + ypos_offset
        d_ypos = (ymax - ymin) / 12
        ax2.text(x=xpos, y=ypos, s=r"$B_\star$    = " + str(B_star) + " G ", fontsize='small')
        ax2.text(x=xpos, y=ypos - d_ypos, s=r"$B_{\rm planet}$ = " + str(B_pl_loc) + r"$B_{\rm Earth}$",
                 fontsize='small')
        ax2.text(x=xpos, y=ypos - 2 * d_ypos, s=r"$\dot{M}$ = " + str(M_star_dot) + "$M_\odot$", fontsize='small')

    def save_plot(self):
        plt.tight_layout()
        outfilePDF = os.path.join(outdir, outfile + ".pdf")
        plt.savefig(outfilePDF)
        outfilePNG = os.path.join(outdir, outfile + ".png")
        plt.savefig(outfilePNG)
        plt.close()


if __name__ == '__main__':
    plotter = Plotter()
    plotter.plot()
