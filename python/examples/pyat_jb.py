import os
import at
import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["figure.figsize"] = (9.0, 6.0)


def prt_ps(str, ps):
    print(str, end='')
    for x in ps:
        #print(x)
        print('{:11.3e}'.format(x), end='')
    print()


def prt_mat(str, a):
    print(str, end='')
    for j in range(0, len(M)):
        for k in range(0, len(M[0])):
            print('{:14.6e}'.format(a[j][k]), end='')
        print()
    

def get_optics(ring):
    [elemdata0, beamdata, elemdata] = \
        at.get_optics(ring, range(len(ring)+1), get_chrom=True)
    print('\n  nu  = [{:5.3f}, {:5.3f}]'
          .format(beamdata.tune[0], beamdata.tune[1]))
    print('  ksi = [{:5.3f}, {:5.3f}]'
          .format(beamdata.chromaticity[0], beamdata.chromaticity[1]))
    [M, _] = at.find_m66(ring, 0)
    prt_mat('\nPoincaré map:\n', M)
    print('\n', at.radiation_parameters(ring))


def example(ring):
    if False:
        print('\n', dir(ring))

    if False:
        for elem in ring:
            print(elem)

    if False:
        refq1 = at.get_cells(ring, at.checktype(at.Quadrupole))
        print('\n', refq1)
        print('\n', list(ring[refq1]))

    if False:
        print('\nelemdata.shape:', elemdata.shape)
        print('elemdata.fields:')
        for fld in elemdata.dtype.fields.keys():
            print(fld)

    if not False:
        ps = np.array([1e-3, -1e-3, 0e0, 0e0, 0e0, 0e0])
        prt_ps('\n', ps)
        at.lattice_pass(sup_per, ps, 1)
        prt_ps('', ps)

    if False:
        fig = plt.figure()
        plt.plot(elemdata.s_pos, elemdata.beta)
        plt.xlabel('s [m]')
        plt.ylabel(r'$\beta$ [m]')
        plt.show()

    if False:
        fig = plt.figure()
        ring.plot_beta()
        plt.show()

    if False:
        print('\n', ring[0])
        ring[0].PassMethod = "CavityPass"
        print('\n', ring[0])

    if False:
        ring.radiation_on(ring)
        [_, beamdata, _] = at.ohmi_envelope(ring)
        print('\nbeamdata.fields:')
        for fld in beamdata.dtype.fields.keys():
            print(fld)


lat_dir = os.environ['LAT']

lat_name = {
    'hmba'      : 'test_matlab/hmba.mat',
    'esrf_ebs'  : 'test_matlab/err.mat',
    'bessy-iii' : lat_dir+'/BESSY-III/NoTG-TGRB-B60-6bend-6sx_JB_tracy.lat'
}

lat = 'bessy-iii'
if lat != 'bessy-iii':
    sup_per = at.load_mat(lat_name[lat])
else:
    sup_per = \
        at.load_tracy(lat_name[lat], harmonic_number=538, lattice_key='cell')

if False:
    for elem in sup_per:
        print(elem)

print('\n', sup_per)

sup_per.radiation_off()
[M, _] = at.find_m66(sup_per, 0e0)
prt_mat('\nPoincaré map:\n', M)

get_optics(sup_per)

#example(sup_per)
