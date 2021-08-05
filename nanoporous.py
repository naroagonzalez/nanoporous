import sisl
from hubbard import HubbardHamiltonian, sp2, density, plot
import numpy as np
import os

# Build sisl.Geometry object
geom = sisl.geom.zgnr(8 )
geom = geom.repeat(16  ,axis=0)

geom = geom.remove([240,241,242,253,254,255])
geom = geom.remove([224,225,238,239])
geom = geom.remove([x+214 for x in range(4)])
geom = geom.remove([x+197 for x in range(6)])
geom = geom.remove([x+181 for x in range(6)])
geom = geom.remove([x+166 for x in range(4)])
geom = geom.remove([144,145,158,159])
geom = geom.remove([128,129,130,141,142,143])
geom = geom.remove([112,113,114,115,125,126,127])
geom = geom.remove([96,97,111])
geom = geom.remove([x+87 for x in range(3)])
geom = geom.remove([x+69 for x in range(7)])
geom = geom.remove([x+53 for x in range(7)])
geom = geom.remove([x+39 for x in range(3)])
geom = geom.remove([16,17,31])
geom = geom.remove([x+13 for x in range(3)])
geom = geom.remove([x for x in range(4)])

# # Plot geometry of the unit cell
p = plot.GeometryPlot(geom, cmap='Greys',figsize=(19,40))
p.annotate(size=10);

geom.sc.set_nsc([3,3,1])

# Build sisl.Hamiltonian object using the sp2 function
H_elec = sp2(geom, t1=2.7, t2=0.2, t3=0.18)

# Build the HubbardHamiltonian object with U=3. eV, Main Field Hamiltonian
MFH_elec = HubbardHamiltonian(H_elec, U=3.)

MFH_elec.random_density()

# Converge until a tolerance of tol, print info for each 10 completed iterations
dn = MFH_elec.converge(density.calc_n, tol=1e-7, print_info=True, steps=10)

# E_level = MFH_elec.Etot


np.save("MFH-finite20.npy", MFH_elec)

# p = plot.SpinPolarization(MFH_elec, colorbar=True, vmax=0.4, vmin=-0.4, figsize=(40,80))
# p.savefig('SpinPolarization.png')

# p = plot.Bandstructure(MFH_elec)
# p.savefig('BandStructure.png')

# ev, evec = MFH_elec.eigh(eigvals_only=False, spin=0)

# print('Eigenvalues:')
# print(ev)

# # Plot wavefunctions
# for i in range(6):
#     p = plot.Wavefunction(MFH_elec, 500*evec[:, i], figsize=(10, 3))
#     p.set_title('State %i' % i)
#     p.annotate()
#     p.savefig('state%i.pdf'%i)

#Get bond order: (bonded e - antibonded e)/2
bo = MFH_elec.get_bond_order(format='dense')

# print('Bond order:')
np.save("H_matrix.npy", bo)
