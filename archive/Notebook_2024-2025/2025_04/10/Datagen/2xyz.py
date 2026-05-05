from ase.io import read, write

atoms = read('thiophene_gaussview_default.gjf')
write('thiophene.extxyz', atoms)