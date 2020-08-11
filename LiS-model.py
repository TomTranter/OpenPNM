# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 21:09:12 2020

@author: tom
"""

import openpnm as op
from openpnm.phases import mixtures
import openpnm.models.physics as pm
import numpy as np
import matplotlib.pyplot as plt

net = op.network.Cubic([60, 1, 1], spacing=3e-5)
geo = op.geometry.StickAndBall(network=net)
pore_d = op.models.misc.constant
throat_d = op.models.misc.constant
geo.add_model(propname='pore.diameter', model=pore_d, value=1e-5)
geo.add_model(propname='throat.diameter', model=throat_d, value=5e-6)
geo.regenerate_models()
# op.topotools.plot_vpython(net)
elec = mixtures.GenericMixture(network=net, name='electrolyte')
S_8 = mixtures.GenericSpecies(network=net, name='S_8')
S_4 = mixtures.GenericSpecies(network=net, name='S_4')
S_2 = mixtures.GenericSpecies(network=net, name='S_2')
S_1 = mixtures.GenericSpecies(network=net, name='S_1')
Li = mixtures.GenericSpecies(network=net, name='Li')
A = mixtures.GenericSpecies(network=net, name='A')
# anion = mixtures.GenericSpecies(network=net, name='anion')
S_molecular_weight = 32.07e-3 # kg/mol
Li_molecular_weight = 6.94e-3 # kg/mol
# High plateau Sulfur Species - S8
S_8['pore.molecular_weight'] = S_molecular_weight*8 # kg/mol
S_8['pore.diffusivity'] = 1.0e-09  # m2/s
S_8['throat.diffusivity'] = 1.0e-09  # m2/s
S_8['pore.valence'] = -2
S_8['throat.valence'] = -2
# Low plateau Sulfur Species - S4
S_4['pore.molecular_weight'] = S_molecular_weight*4 # kg/mol
S_4['pore.diffusivity'] = 1.0e-10  # m2/s
S_4['throat.diffusivity'] = 1.0e-10  # m2/s
S_4['pore.valence'] = -2
S_4['throat.valence'] = -2
# Final Sulfur Species - S2
S_2['pore.molecular_weight'] = S_molecular_weight*2 # kg/mol
S_2['pore.diffusivity'] = 1.0e-10  # m2/s
S_2['throat.diffusivity'] = 1.0e-10  # m2/s
S_2['pore.valence'] = -2
S_2['throat.valence'] = -2
# Final Sulfur Species - S2
S_1['pore.molecular_weight'] = S_molecular_weight*1 # kg/mol
S_1['pore.diffusivity'] = 1.0e-10  # m2/s
S_1['throat.diffusivity'] = 1.0e-10  # m2/s
S_1['pore.valence'] = -2
S_1['throat.valence'] = -2
# Lithium ion
Li['pore.molecular_weight'] = Li_molecular_weight*1 # kg/mol
Li['pore.diffusivity'] = 1.0e-9  # m2/s
Li['throat.diffusivity'] = 1.0e-9  # m2/s
Li['pore.valence'] = +1
Li['throat.valence'] = +1
# Anion
A['pore.molecular_weight'] = Li_molecular_weight*1 # kg/mol
A['pore.diffusivity'] = 1.0e-9  # m2/s
A['throat.diffusivity'] = 1.0e-9  # m2/s
A['pore.valence'] = -1
A['throat.valence'] = -1

elec.set_component([S_8, S_4, S_2, S_1, Li, A])
elec.set_concentration(component=S_8, values=19.0)
elec.set_concentration(component=S_4, values=0.0)
elec.set_concentration(component=S_2, values=0.0)
elec.set_concentration(component=S_1, values=0.0)
elec.set_concentration(component=Li, values=1000.0)
elec.set_concentration(component=A, values=1000.0)
elec.update_mole_fractions()
# elec.add_model(propname='pore.salt_concentration',
#                model=mods.misc.summation,
#                props=['pore.concentration.Na_'+self.name,
#                       'pore.concentration.Cl_'+self.name])
# elec.add_model(propname='pore.salinity',
#                model=mods.phases.mixtures.salinity,
#                concentration='pore.concentration.Na_'+self.name)
# elec.add_model(propname='pore.mass_density',
#                model=mods.phases.density.water,
#                salinity='pore.salinity')
# elec.add_model(propname='pore.viscosity',
#                model=mods.phases.viscosity.water,
#                salinity='pore.salinity')
# elec.add_model(propname='pore.permittivity',
#                model=mods.phases.permittivity.water)

# physics
phys = op.physics.GenericPhysics(network=net, phase=elec, geometry=geo)

# flow = op.models.physics.hydraulic_conductance.hagen_poiseuille
# phys.add_model(propname='throat.hydraulic_conductance',
#                pore_viscosity='pore.viscosity',
#                throat_viscosity='throat.viscosity',
#                model=flow, regen_mode='normal')

phys['throat.hydraulic_conductance'] = 1e-6

eA_dif = op.models.physics.diffusive_conductance.ordinary_diffusion
ad_dif_mig = op.models.physics.ad_dif_mig_conductance.ad_dif_mig
s_scheme = 'powerlaw'
for ion in elec.components.values():
    phys.add_model(propname='throat.diffusive_conductance.' + ion.name,
                   pore_diffusivity='pore.diffusivity.' + ion.name,
                   throat_diffusivity='throat.diffusivity.' + ion.name,
                   model=eA_dif, regen_mode='normal')
    phys.add_model(propname='throat.ad_dif_mig_conductance.' + ion.name,
                   pore_pressure='pore.pressure', model=ad_dif_mig,
                   ion=ion.name, s_scheme=s_scheme)
current = op.models.physics.ionic_conductance.electroneutrality
phys.add_model(propname='throat.ionic_conductance', ions=[S_8.name, S_4.name, S_2.name, S_1.name],
               model=current, regen_mode='normal')

# settings for algorithms
setts1 = {'solver_max_iter': 5, 'solver_tol': 1e-08, 'solver_rtol': 1e-08,
          'nlin_max_iter': 10, 'cache_A': False, 'cache_b': False}
setts2 = {'g_tol': 1e-4, 'g_max_iter': 4, 't_output': 500, 't_step': 500,
          't_final': 2000, 't_scheme': 'implicit'}

# algorithms
sf = op.algorithms.StokesFlow(network=net, phase=elec, settings=setts1)
sf.set_value_BC(pores=net.pores('back'), values=0.00)
sf.set_value_BC(pores=net.pores('front'), values=0.00)
sf.run()
elec.update(sf.results())

p = op.algorithms.TransientIonicConduction(network=net, phase=elec,
                                           settings=setts1)
p.set_value_BC(pores=net.pores('back'), values=0.1)
p.set_value_BC(pores=net.pores('front'), values=0.00)
p.settings['charge_conservation'] = 'electroneutrality'

ion_TNP_algs = {}

for ion, obj in elec.components.items():
    eA = op.algorithms.TransientNernstPlanck(network=net, phase=elec, ion=obj.name,
                                             settings=setts1, name=ion+'_TNP')
    eA.set_rate_BC(pores=net.pores('back'), values=0)
    eA.set_rate_BC(pores=net.pores('front'), values=0)
    ion_TNP_algs[eA.name] = eA

# Reactions
net["pore.reaction_area"] = np.random.rand(net.Np)
phys['pore.solid_voltage'] = 1.1
phys['pore.open_circuit_voltage'] = 1.2
elec['pore.electrolyte_concentration'] = np.random.rand(net.Np)
BV_params = {
    "z": 4,
    "j0": 1e-3,
    "c_ref": 1000,
    "alpha_anode": 0.4,
    "alpha_cathode": 0.6
}
phys.add_model(propname='pore.rxn_BV_c',
               model=pm.generic_source_term.butler_volmer_c,
               X="pore.concentration.Li", 
               electrolyte_voltage='pore.potential',
               **BV_params)
phys.add_model(propname='pore.rxn_BV_v',
               model=pm.generic_source_term.butler_volmer_v,
               X="pore.potential", **BV_params)


it = op.algorithms.TransientNernstPlanckMultiphysicsSolver(network=net,
                                                           phase=elec,
                                                           settings=setts2)
it.setup(potential_field=p.name, ions=list(ion_TNP_algs.keys()))
it.run()


elec.update(p.results())
for alg in ion_TNP_algs.values():
    elec.update(alg.results())

plt.figure()
for conc in list(ion_TNP_algs['Li_TNP'].results().values()):
    plt.plot(conc)
