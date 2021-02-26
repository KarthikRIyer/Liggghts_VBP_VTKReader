from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import math
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import numpy as np
# from sklearn.metrics import r2_score
import sys
import os

n = len(sys.argv)
if n < 2:
    print("Enter path of post folder!!")
    exit(0)


def get_avg_vel(filepath):
    f = open(filepath, 'r')
    lines = f.readlines()

    avg_vel_bin_x = []
    avg_vel_bin_y = []
    avg_vel_bin_z = []
    avg_vel_bin_x_type0 = []
    avg_vel_bin_y_type0 = []
    avg_vel_bin_z_type0 = []
    avg_vel_bin_x_type1 = []
    avg_vel_bin_y_type1 = []
    avg_vel_bin_z_type1 = []
    ts = 0
    avg_v = 0
    smi = 0

    for line in lines:
        if line.find('AVERAGE VELOCITY') != -1:
            avg_v = float(line[line.rfind('=') + 2:])
            ts = float(
                filepath[filepath.rindex('dump') + 4:filepath.rindex('.vtk.vel')])
        elif line.find('AVG VEL BIN X') != -1:
            words = line.split()
            avg_vel_bin_x.append(float(words[-1]))
        elif line.find('AVG VEL BIN Y') != -1:
            words = line.split()
            avg_vel_bin_y.append(float(words[-1]))
        elif line.find('AVG VEL BIN Z') != -1:
            words = line.split()
            avg_vel_bin_z.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN X 0') != -1:
            words = line.split()
            avg_vel_bin_x_type0.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN Y 0') != -1:
            words = line.split()
            avg_vel_bin_y_type0.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN Z 0') != -1:
            words = line.split()
            avg_vel_bin_z_type0.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN X 1') != -1:
            words = line.split()
            avg_vel_bin_x_type1.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN Y 1') != -1:
            words = line.split()
            avg_vel_bin_y_type1.append(float(words[-1]))
        elif line.find('AVG VEL TYPE BIN Z 1') != -1:
            words = line.split()
            avg_vel_bin_z_type1.append(float(words[-1]))
        elif line.find('SMI') != -1:
            words = line.split()
            smi = float(words[-1])
    avg_vel_bin_x.reverse()
    avg_vel_bin_y.reverse()
    avg_vel_bin_z.reverse()
    avg_vel_bin_x_type0.reverse()
    avg_vel_bin_x_type1.reverse()
    avg_vel_bin_y_type0.reverse()
    avg_vel_bin_y_type1.reverse()
    avg_vel_bin_z_type0.reverse()
    avg_vel_bin_z_type1.reverse()

    return avg_v, ts, avg_vel_bin_x, avg_vel_bin_y, avg_vel_bin_z, avg_vel_bin_x_type0, avg_vel_bin_x_type1, avg_vel_bin_y_type0, avg_vel_bin_y_type1, avg_vel_bin_z_type0, avg_vel_bin_z_type1, smi

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


def is_vel_file(velpath):
    if not os.path.isfile(velpath):
        return False

    filename, file_extension = os.path.splitext(velpath)
    return file_extension == '.vel'


for i in range(1, n):
    path = sys.argv[i]
    filename, file_extension = os.path.splitext(path)
    print('Extracting data')
    if os.path.isfile(path):
        pass
    else:
        avg_vel_list = []
        timesteps_list = []
        avg_vel_bin_x_list = []
        avg_vel_bin_y_list = []
        avg_vel_bin_z_list = []
        avg_vel_bin_x_type0_list = []
        avg_vel_bin_y_type0_list = []
        avg_vel_bin_z_type0_list = []
        avg_vel_bin_x_type1_list = []
        avg_vel_bin_y_type1_list = []
        avg_vel_bin_z_type1_list = []
        smi_list = []
        path = os.path.join(path, '_post_processed')
        ppfiles = [os.path.join(path, f) for f in os.listdir(path) if
                   is_vel_file(os.path.join(path, f))]
        # ppfiles = sorted(ppfiles)
        for ppfile in ppfiles:
            avg_vel, timestep, avg_vel_bin_x, avg_vel_bin_y, avg_vel_bin_z, \
            avg_vel_bin_x_type0, avg_vel_bin_x_type1, avg_vel_bin_y_type0, \
            avg_vel_bin_y_type1, avg_vel_bin_z_type0, \
            avg_vel_bin_z_type1, smi = get_avg_vel(ppfile)
            avg_vel_list.append(avg_vel)
            timesteps_list.append(timestep)
            avg_vel_bin_x_list.append(avg_vel_bin_x)
            avg_vel_bin_y_list.append(avg_vel_bin_y)
            avg_vel_bin_z_list.append(avg_vel_bin_z)
            avg_vel_bin_x_type0_list.append(avg_vel_bin_x_type0)
            avg_vel_bin_x_type1_list.append(avg_vel_bin_x_type1)
            avg_vel_bin_y_type0_list.append(avg_vel_bin_y_type0)
            avg_vel_bin_y_type1_list.append(avg_vel_bin_y_type1)
            avg_vel_bin_z_type0_list.append(avg_vel_bin_z_type0)
            avg_vel_bin_z_type1_list.append(avg_vel_bin_z_type1)
            smi_list.append(smi)

        timesteps_list, avg_vel_list, \
        avg_vel_bin_x_list, avg_vel_bin_y_list, \
        avg_vel_bin_z_list, avg_vel_bin_x_type0_list, \
        avg_vel_bin_x_type1_list, avg_vel_bin_y_type0_list, \
        avg_vel_bin_y_type1_list, avg_vel_bin_z_type0_list, avg_vel_bin_z_type1_list, \
        smi_list = (
            list(t) for t in
            zip(*sorted(
                zip(timesteps_list, avg_vel_list, avg_vel_bin_x_list,
                    avg_vel_bin_y_list, avg_vel_bin_z_list,
                    avg_vel_bin_x_type0_list, avg_vel_bin_x_type1_list,
                    avg_vel_bin_y_type0_list, avg_vel_bin_y_type1_list,
                    avg_vel_bin_z_type0_list, avg_vel_bin_z_type1_list, smi_list))))
        print('Extracted data, now plotting...')
        fig_dpi = 80
        # PLOT AVG_VEL vs TIME
        plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
        # if len(timesteps_list) != 0:
        #     fitmodel = np.poly1d(
        #         np.polyfit(timesteps_list, avg_vel_list, 4))
        #     fitline = np.linspace(timesteps_list[0], timesteps_list[-1],
        #                           len(timesteps_list) * 10)
        #     plt.plot(fitline, fitmodel(fitline))
        plt.plot(timesteps_list, avg_vel_list)
        plt.xlabel('Timesteps')
        plt.ylabel('Average Velocity (m/s)')
        plt.title('Avg Velocity vs timesteps')
        plt.savefig(os.path.join(path, 'vel_time.png'), dpi=fig_dpi,
                    bbox_inches='tight')
        # print('Saved file' + filepath + '.png')
        plt.close('all')

        # PLOT SMI vs TIME
        plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
        # if len(timesteps_list) != 0:
        #     fitmodel = np.poly1d(
        #         np.polyfit(timesteps_list, avg_vel_list, 4))
        #     fitline = np.linspace(timesteps_list[0], timesteps_list[-1],
        #                           len(timesteps_list) * 10)
        #     plt.plot(fitline, fitmodel(fitline))
        plt.plot(timesteps_list, smi_list)
        plt.xlabel('Timesteps')
        plt.ylabel('SMI')
        plt.title('SMI vs timesteps')
        plt.savefig(os.path.join(path, 'smi_time.png'), dpi=fig_dpi,
                    bbox_inches='tight')
        # print('Saved file' + filepath + '.png')
        plt.close('all')

        for (ts, avg_vel_bin_x, avg_vel_bin_y, avg_vel_bin_z, avg_vel_bin_x_type0,
             avg_vel_bin_x_type1, avg_vel_bin_y_type0, avg_vel_bin_y_type1,
             avg_vel_bin_z_type0, avg_vel_bin_z_type1) in zip(timesteps_list,
                                                              avg_vel_bin_x_list,
                                                              avg_vel_bin_y_list,
                                                              avg_vel_bin_z_list,
                                                              avg_vel_bin_x_type0_list,
                                                              avg_vel_bin_x_type1_list,
                                                              avg_vel_bin_y_type0_list,
                                                              avg_vel_bin_y_type1_list,
                                                              avg_vel_bin_z_type0_list,
                                                              avg_vel_bin_z_type1_list):
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_x) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_x)), avg_vel_bin_x, 4))
                fitline = np.linspace(0, len(avg_vel_bin_x) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_x)), avg_vel_bin_x)
            plt.xlabel('X bins')
            plt.ylabel('Average Velocity (m/s) Both particles')
            plt.title('Avg Velocity vs X-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_x_bin.png'), dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_y) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_y)), avg_vel_bin_y, 4))
                fitline = np.linspace(0, len(avg_vel_bin_y) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_y)), avg_vel_bin_y)
            plt.xlabel('Y bins')
            plt.ylabel('Average Velocity (m/s) Both particles')
            plt.title('Avg Velocity vs Y-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_y_bin.png'), dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_z) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_z)), avg_vel_bin_z, 4))
                fitline = np.linspace(0, len(avg_vel_bin_z) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_z)), avg_vel_bin_z)
            plt.xlabel('Z bins')
            plt.ylabel('Average Velocity (m/s) Both particles')
            plt.title('Avg Velocity vs Z-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_z_bin.png'), dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            # fine avg vel vs x
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_x_type0) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_x_type0)), avg_vel_bin_x_type0, 4))
                fitline = np.linspace(0, len(avg_vel_bin_x_type0) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_x_type0)), avg_vel_bin_x_type0)
            plt.xlabel('X bins')
            plt.ylabel('Average Velocity Fines (m/s)')
            plt.title('Avg Velocity Fines vs X-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_x_type0_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')
            # coarse avg vel vs x
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_x_type1) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_x_type1)), avg_vel_bin_x_type1, 4))
                fitline = np.linspace(0, len(avg_vel_bin_x_type1) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_x_type1)), avg_vel_bin_x_type1)
            plt.xlabel('X bins')
            plt.ylabel('Average Velocity Coarse (m/s)')
            plt.title('Avg Velocity Coarse vs X-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_x_type1_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            # fine avg vel vs y
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_y_type0) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_y_type0)), avg_vel_bin_y_type0, 4))
                fitline = np.linspace(0, len(avg_vel_bin_y_type0) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_y_type0)), avg_vel_bin_y_type0)
            plt.xlabel('Y bins')
            plt.ylabel('Average Velocity Fines (m/s)')
            plt.title('Avg Velocity Fines vs Y-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_y_type0_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')
            # coarse avg vel vs y
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_y_type1) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_y_type1)), avg_vel_bin_y_type1, 4))
                fitline = np.linspace(0, len(avg_vel_bin_y_type1) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_y_type1)), avg_vel_bin_y_type1)
            plt.xlabel('Y bins')
            plt.ylabel('Average Velocity Coarse (m/s)')
            plt.title('Avg Velocity Coarse vs Y-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_y_type1_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            # fine avg vel vs z
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_z_type0) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_z_type0)), avg_vel_bin_z_type0, 4))
                fitline = np.linspace(0, len(avg_vel_bin_z_type0) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_z_type0)), avg_vel_bin_z_type0)
            plt.xlabel('Z bins')
            plt.ylabel('Average Velocity Fines (m/s)')
            plt.title('Avg Velocity Fines vs Z-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_z_type0_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')
            # coarse avg vel vs z
            plt.figure(figsize=(1920 / fig_dpi, 1080 / fig_dpi), dpi=fig_dpi)
            if len(avg_vel_bin_z_type1) != 0:
                fitmodel = np.poly1d(
                    np.polyfit(range(len(avg_vel_bin_z_type1)), avg_vel_bin_z_type1, 4))
                fitline = np.linspace(0, len(avg_vel_bin_z_type1) - 1, 100)
                plt.plot(fitline, fitmodel(fitline))
            plt.plot(range(len(avg_vel_bin_z_type1)), avg_vel_bin_z_type1)
            plt.xlabel('Z bins')
            plt.ylabel('Average Velocity Coarse (m/s)')
            plt.title('Avg Velocity Coarse vs Z-Bin')
            plt.savefig(os.path.join(path, str(ts) + '.avg_vel_z_type1_bin.png'),
                        dpi=fig_dpi,
                        bbox_inches='tight')
            plt.close('all')

            print(str(ts) + ' processed')
