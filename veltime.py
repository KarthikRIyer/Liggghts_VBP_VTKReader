from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import math
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import sys
import os

n = len(sys.argv)
if n < 2:
    print("Enter path of post folder!!")
    exit(0)


def get_avg_vel(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()

    for line in lines:
        if line.find('AVERAGE VELOCITY') != -1:
            return float(line[line.rfind('=') + 2:]), float(filepath[filepath.rindex(
                'dump') + 4:filepath.rindex('.vtk.postprocessed')])

    return 0.0, 0.0

    # colors = [math.sqrt(a * a + b * b) for a, b in zip(vx, vz)]
    # norm = Normalize()
    # norm.autoscale(colors)
    # colormap = cm.inferno
    #
    # sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    # sm.set_array([])
    # fig_dpi = 80
    # # plt.figure(figsize=(800 / fig_dpi, 6 / fig_dpi), dpi=fig_dpi)
    # fig, ax = plt.subplots(figsize=(1920 / fig_dpi, 1080 / fig_dpi))
    # ax.quiver(px, pz, vx, vz, scale=15, color=colormap(norm(colors)))
    # cbar = plt.colorbar(sm)
    # cbar.ax.get_yaxis().labelpad = 15
    # cbar.ax.set_ylabel('Velocity projection in X-Z plane', rotation=270)
    # ax.set_ylabel('Z')
    # ax.set_xlabel('X')
    # plt.savefig(filepath + '.png', dpi=fig_dpi, bbox_inches='tight')
    # print('Saved file' + filepath + '.png')
    # plt.close('all')


def is_post_processed_file(pppath):
    if not os.path.isfile(pppath):
        return False

    filename, file_extension = os.path.splitext(pppath)
    return file_extension == '.postprocessed'


for i in range(1, n):
    path = sys.argv[i]
    filename, file_extension = os.path.splitext(path)

    if os.path.isfile(path):
        pass
    else:
        avg_vel_list = []
        timesteps_list = []
        path = os.path.join(path, '_post_processed')
        ppfiles = [os.path.join(path, f) for f in os.listdir(path) if
                   is_post_processed_file(os.path.join(path, f))]
        # ppfiles = sorted(ppfiles)
        for ppfile in ppfiles:
            avg_vel, timestep = get_avg_vel(ppfile)
            avg_vel_list.append(avg_vel)
            timesteps_list.append(timestep)

        timesteps_list, avg_vel_list = (list(t) for t in
                                        zip(*sorted(zip(timesteps_list, avg_vel_list))))
        # colors = [math.sqrt(a * a + b * b) for a, b in zip(vx, vz)]
        # norm = Normalize()
        # norm.autoscale(colors)
        # colormap = cm.inferno
        #
        # sm = cm.ScalarMappable(cmap=colormap, norm=norm)
        # sm.set_array([])
        fig_dpi = 80
        plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
        plt.plot(timesteps_list, avg_vel_list)
        # plt.show()
        # fig, ax = plt.subplots(figsize=(1920 / fig_dpi, 1080 / fig_dpi))
        # ax.quiver(px, pz, vx, vz, scale=15, color=colormap(norm(colors)))
        # cbar = plt.colorbar(sm)
        # cbar.ax.get_yaxis().labelpad = 15
        # cbar.ax.set_ylabel('Velocity projection in X-Z plane', rotation=270)
        # ax.set_ylabel('Z')
        # ax.set_xlabel('X')
        plt.xlabel('Timesteps')
        plt.ylabel('Average Velocity (m/s)')
        plt.title('Avg Velocity vs timesteps')
        plt.savefig(os.path.join(path, 'vel_time.png'), dpi=fig_dpi,
                    bbox_inches='tight')
        # print('Saved file' + filepath + '.png')
        # plt.close('all')
