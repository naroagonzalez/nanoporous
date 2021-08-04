import sisl
from hubbard import HubbardHamiltonian, sp2, density, plot
import numpy as np
import os

# Build sisl.Geometry object
geom = sisl.geom.zgnr(10 )
geom = geom.repeat(18  ,axis=0)

geom = geom.remove(359)
geom = geom.remove([x+353 for x in range(3)])
geom = geom.remove([x+342 for x in range(8)])
geom = geom.remove(341)
geom = geom.remove([x+324 for x in range(7)])
geom = geom.remove(323)
geom = geom.remove([x+306 for x in range(7)])
geom = geom.remove(305)
geom = geom.remove([x+288 for x in range(7)])
geom = geom.remove(287)
geom = geom.remove([x+270+3 for x in range(3)])
geom = geom.remove(269)
geom = geom.remove([x+252+11 for x in range(3) ])
geom = geom.remove([x+234+9 for x in range(9) ])
geom = geom.remove([x+216+9 for x in range(9) ])
geom = geom.remove([x+198+9 for x in range(9) ])
geom = geom.remove([x+180+9 for x in range(9) ])
geom = geom.remove([x+180+3 for x in range(3) ])

geom = geom.remove([x+162 for x in range(18) if x!=10 and x!=15 ])
geom = geom.remove([x+144 for x in range(7) ])
geom = geom.remove([x+126 for x in range(7) ])
geom = geom.remove([x+108 for x in range(7) ])
geom = geom.remove([x+90 for x in range(7) if x!=2 ])
geom = geom.remove([x+72+11 for x in range(7) if x!=4 ])
geom = geom.remove([x+54+10 for x in range(8) ])
geom = geom.remove([x+36+10 for x in range(8) ])
geom = geom.remove([x+18+10 for x in range(8) ])
geom = geom.remove([x for x in range(18) if x !=2 and x!=7])

# geom.tile([1,2,1])
# Plot geometry of the unit cell
p = plot.GeometryPlot(geom, cmap='Greys',figsize=(19,40))
# p.annotate(size=10);

geom.sc.set_nsc([1,1,0])

# Build sisl.Hamiltonian object using the sp2 function
H_elec = sp2(geom, t1=2.7, t2=0.2, t3=0.18)

# Build the HubbardHamiltonian object with U=3. eV, Main Field Hamiltonian
MFH_elec = HubbardHamiltonian(H_elec, U=3.)

MFH_elec.random_density()

# Converge until a tolerance of tol, print info for each 10 completed iterations
dn = MFH_elec.converge(density.calc_n, tol=1e-7, print_info=True, steps=10)

E_level = MFH_elec.Etot


np.save("MFH-finite20.npy", MFH_elec)

p = plot.SpinPolarization(MFH_elec, colorbar=True, vmax=0.4, vmin=-0.4, figsize=(40,80))
p.savefig('SpinPolarization.png')

p = plot.Bandstructure(MFH_elec, ymax=3)
p.savefig('BandStructure.png')

ev, evec = MFH_elec.eigh(eigvals_only=False, spin=0)

print('Eigenvalues:')
print(ev)

# Plot wavefunctions
for i in range(6):
    p = plot.Wavefunction(MFH_elec, 500*evec[:, i], figsize=(10, 3))
    p.set_title('State %i' % i)
    p.annotate()
    p.savefig('state%i.pdf'%i)

#Get bond order: (bonded e - antibonded e)/2
bo = MFH_elec.get_bond_order(format='csr')

print('Bond order:')
print(bo)
