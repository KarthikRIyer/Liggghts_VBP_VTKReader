import re
import sys
import os
import matplotlib.pyplot as plt

binwise_xf_exp = [0.842069701280228, 0.710978160416363, 0.631664146851759,
                  0.585001669214415, 0.511177503472495, 0.527220980650584,
                  0.485465581561419, 0.28474028886696, 0.135134146898853, 0]

n = len(sys.argv)
if n < 2:
    print("Enter path of post folder!!")
    exit(0)


def is_pp_file(pp_path):
    if not os.path.isfile(pp_path):
        return False

    filename, file_extension = os.path.splitext(pp_path)
    return file_extension == '.postprocessed'


path = sys.argv[1]
filename, file_extension = os.path.splitext(path)
if os.path.isfile(path):
    pass
else:
    path = os.path.join(path, '_post_processed')
    ppfiles = [os.path.join(path, f) for f in os.listdir(path) if
               is_pp_file(os.path.join(path, f))]
    ts_list = [float(
        filepath[filepath.rindex('dump') + 4:filepath.rindex('.vtk.postprocessed')]) for
        filepath in ppfiles]
    ts_list, ppfiles = (
        list(t) for t in
        zip(*sorted(
            zip(ts_list, ppfiles))))
    final_pp_file_path = ppfiles[-1]
    final_pp_file = open(final_pp_file_path, 'r')
    final_pp_file_lines = final_pp_file.readlines()
    binwise_xf_sim = []
    for i in range(len(final_pp_file_lines)):
        if final_pp_file_lines[i].startswith('Bin'):
            binwise_xf_sim.append(float(
                re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
                           final_pp_file_lines[i + 1])[-1]))
    binwise_xf_sim.reverse()

    fig_dpi = 80

    plt.figure(figsize=(800 / fig_dpi, 600 / fig_dpi), dpi=fig_dpi)
    # if len(timesteps_list) != 0:
    #     fitmodel = np.poly1d(
    #         np.polyfit(timesteps_list, avg_vel_list, 4))
    #     fitline = np.linspace(timesteps_list[0], timesteps_list[-1],
    #                           len(timesteps_list) * 10)
    #     plt.plot(fitline, fitmodel(fitline))
    y_pos = [0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95]
    plt.scatter(binwise_xf_exp, y_pos, label='Experimental')
    plt.scatter(binwise_xf_sim, y_pos, label='Simulation')
    plt.legend(loc="upper right")
    plt.xlabel('Fine mass fraction (xi/xo)')
    plt.ylabel('Normalized axial position (y/h)')
    plt.savefig(os.path.join(path, 'bin_xf.png'), dpi=fig_dpi,
                bbox_inches='tight')
    # print('Saved file' + filepath + '.png')
    plt.close('all')

    abs_error = 0

    for xf_exp, xf_sim in zip(binwise_xf_exp, binwise_xf_sim):
        if xf_exp == 0:
            continue
        abs_error = abs(xf_sim-xf_exp)/xf_exp
        # abs_error += 2 * ((xf_sim - xf_exp) / (abs(xf_exp) + abs(xf_sim)))
    abs_error = abs_error / len(binwise_xf_exp)
    print(abs_error * 100)
