# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import

import os, os.path
import numpy as np
from numpy.linalg import inv
import ielu.geometry

# import csvkit as csv # csvkit handles unicode better
import csv


class ElectrodeSet(object):
    def __init__(self):
        self.elect = dict()

    def read_ielu_csv(self, csvfn):
        elect = self.elect
        with open(csvfn, 'rb') as origfile:
            R = csv.reader(origfile)
            for row in R:
                print(row)
                name = row[0]  # name of electrode
                elect[name] = np.array([float(xx) for xx in row[1:]])

        self.coordsys = 'surfaceRAS'


def electrode_csv2fspointset(csvfn, mrifn):
    """convert an entire ielu electrode file into an unlabeled freesurfer pointset for viewing in freeview
    note """
    if csvfn[-4:] == '.csv':
        # then its ok
        outfn = csvfn[:-4] + '.dat'
    else:
        print("unexpected input: %s should be ielu csv electrode file" % csvfn)
        return
    electrodes = ElectrodeSet()
    electrodes.read_ielu_csv(csvfn)

    # get affine transformations
    vox2ras_tkr = ielu.geometry.get_vox2rasxfm(mrifn, stem='vox2ras-tkr')
    vox2ras = ielu.geometry.get_vox2rasxfm(mrifn, stem='vox2ras')
    # the coordinates system is called "surfaceRAS" in freeview
    rastk2ras = np.dot(vox2ras, inv(vox2ras_tkr))
    # this affine transformation is a simple translation



    locs = electrodes.elect.values()  # get coordinates
    raselect = {}

    # elect is in tkRAS space a translation of scanner (?) RAS

    trans_locs = ielu.geometry.apply_affine(locs, rastk2ras)

    with open(outfn, 'wb+') as newfile:
        W = csv.writer(newfile, delimiter=' ')  # separted by spaces
        for vv in trans_locs:
            row = vv[0], vv[1], vv[2]
            W.writerow(row)
        newfile.write('info\n')
        newfile.write('numpoints %d\n' % len(trans_locs))
        newfile.write('useRealRAS 1\n')


def electrode_csv2mult_fspointset(csvfn, mrifn):
    """convert an entire ielu electrode file into multple unlabeled freesurfer pointset for viewing in freeview
    note. Each electrode base name is given a separate filename

    assume that the electrodes are named by convenstion
    stem_1, stem_2, stem_3 etc
    for example A_1, A_2, A_3 and so on"""
    if csvfn[-4:] == '.csv':
        # then its ok
        outfn = csvfn[:-4] + '.dat'
    else:
        print("unexpected input: %s should be ielu csv electrode file" % csvfn)
        return
    electrodes = ElectrodeSet()
    electrodes.read_ielu_csv(csvfn)

    # get affine transformations
    vox2ras_tkr = ielu.geometry.get_vox2rasxfm(mrifn, stem='vox2ras-tkr')
    vox2ras = ielu.geometry.get_vox2rasxfm(mrifn, stem='vox2ras')
    # the coordinates system is called "surfaceRAS" in freeview
    rastk2ras = np.dot(vox2ras, inv(vox2ras_tkr))
    # this affine transformation is a simple translation

    keys = electrodes.elect.keys()
    keys.sort()
    stems = dict()
    for kk in keys:
        aa = kk.split('_')
        stems[kk] = aa[0]

    for kk in keys:
        print(stems[kk], kk, electrodes.elect[kk])


def main_test():
    Mct = np.array([
        [-0.9375, 0.0000, -0.0000, 127.0940],
        [0.0000, -0.9375, -0.0000, 176.0210],
        [0.0000, 0.0000, 1.0000, -40.3760],
        [0.0000, 0.0000, 0.0000, 1.0000]])

    orig = [
        [-35.2751, 56.7037, -19.4089],
        [-33.6036, 49.9731, -15.8463],
        [-33.3972, 43.8767, -10.8981],
        [-34.6550, 38.8735, -4.4505],
        [-37.8249, 34.6606, 2.7577],
        [-45.0778, 28.6904, 17.0992],
        [-49.3864, 26.5400, 22.3633],
        [-52.0721, 25.0418, 28.4364],
        [-62.0432, 17.5541, 48.4583]]

    # subjects = os.environ['SUBJECTS_DIR']

    # mri= os.path.join(subjects,'../tests/data/brain.mgz')
    mri = '../tests/data/brain.mgz'
    electrode_csv2fspointset('../tests/data/electrodes_jje2016.csv', mri)
    electrode_csv2mult_fspointset('../tests/data/electrodes_jje2016.csv', mri)


if __name__ == '__main__':
    main_test()
