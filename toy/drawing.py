import os
import shutil

from matplotlib import pyplot

pyplot.ioff()

draw_and_save_counter = 0
draw_flag = True
path = os.getcwd()


def erase_img():
    erase_folder(path + "/img/master/")
    erase_folder(path + "/img/dual/")


def erase_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def draw(master, dual, save_name=None):
    fig = pyplot.figure(figsize=(8, 10))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    ax1.title.set_text('Master')
    ax2.title.set_text('Dual')
    pyplot.subplot(ax1)
    master.draw(node_size=600)
    pyplot.subplot(ax2)
    dual.draw(node_size=600)
    if save_name is None:
        pyplot.show()
    else:
        pyplot.savefig(save_name)
