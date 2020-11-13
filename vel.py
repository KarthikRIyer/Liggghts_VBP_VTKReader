from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import math
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import sys
import os

n = len(sys.argv)
if n < 2:
    print("Enter .vel filepath or path of post folder!!")
    exit(0)


def process_file(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()

    px = []
    py = []
    pz = []
    vx = []
    vy = []
    vz = []

    for line in lines:
        x = line.split()
        x = [float(c) for c in x]
        px.append(x[0])
        py.append(x[1])
        pz.append(x[2])
        vx.append(x[3])
        vy.append(x[4])
        vz.append(x[5])

    colors = [math.sqrt(a * a + b * b) for a, b in zip(vx, vz)]
    norm = Normalize()
    norm.autoscale(colors)
    colormap = cm.inferno

    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])
    fig_dpi = 80
    # plt.figure(figsize=(800 / fig_dpi, 6 / fig_dpi), dpi=fig_dpi)
    fig, ax = plt.subplots(figsize=(1920 / fig_dpi, 1080 / fig_dpi))
    ax.quiver(px, pz, vx, vz, scale=15, color=colormap(norm(colors)))
    cbar = plt.colorbar(sm)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Velocity projection in X-Z plane', rotation=270)
    ax.set_ylabel('Z')
    ax.set_xlabel('X')
    plt.savefig(filepath + '.png', dpi=fig_dpi, bbox_inches='tight')
    print('Saved file' + filepath + '.png')
    plt.close('all')


def is_vel_file(velpath):
    if not os.path.isfile(velpath):
        return False

    filename, file_extension = os.path.splitext(velpath)
    return file_extension == '.vel'


for i in range(1, n):
    path = sys.argv[i]
    filename, file_extension = os.path.splitext(path)
    if os.path.isfile(path) and file_extension == '.vel':
        process_file(path)
    elif os.path.isfile(path) and file_extension != '.vel':
        pass
    else:
        path = os.path.join(path, '_post_processed')
        velfiles = [os.path.join(path, f) for f in os.listdir(path) if
                    is_vel_file(os.path.join(path, f))]
        velfiles = sorted(velfiles)
        for velfile in velfiles:
            process_file(velfile)
