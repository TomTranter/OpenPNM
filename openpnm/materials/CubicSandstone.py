import scipy as sp
from openpnm.utils import logging, Project, GenericSettings
from openpnm.network import Cubic
from openpnm.geometry import GenericGeometry
import openpnm.models as mods
logger = logging.getLogger(__name__)


class SandstoneParams():
    r"""
    The following parameters are taken from the paper by Ioannidis and Chatzis,
    Chemical Engineering Science, 1992.
    """
    berea_108 = {'alpha_t': 1.904e-6,
                 'alpha_b': 6.081e-6,
                 'beta_t': 0.536,
                 'beta_b': 1.18,
                 'lattice_constant': 125.6e-6,
                 'b_min': 24.6e-6,
                 'b_max': 70.4e-6,
                 'b_vo': 34.7e-6,
                 'pore_shape': 1.18,  # shape is beta_b
                 'pore_scale': 6.08e-6,  # scale is alpha_b
                 'pore_loc': 24.6e-6,  # loc is b_min
                 'throat_shape': 0.536,  # shape is beta_t
                 'throat_scale': 1.904e-6,  # scale is alpha_t
                 'throat_loc': 0.7e-6}  # loc is ??
    # f = (beta_b/alpha_b)*((x - b_min)/alpha_b)**(beta_b-1)
    # n1 = exp(-((x-b_min)/alpha_b)**(beta_b))
    # n2 = exp(-((b_max-b_min)/alpha_b)**(beta_b))
    # d = 1 - exp(-((b_max-b_min)/alpha_b)**(beta_b))
    # F = (n1 - n2)/d
    boise = {'alpha_t': 13.050e-6,
             'alpha_b': 8.333e-6,
             'beta_t': 0.536,
             'beta_b': 1.18,
             'lattice_constant': 171.2e-6,
             'b_min': 38.7e-6,
             'b_max': 73.9e-6,
             'b_vo': 42.3e-6,
             'pore_shape': 1.75,
             'pore_scale': 10e-6,
             'pore_loc': 40.81e-6,
             'throat_shape': 0.8,
             'throat_scale': 7.5e-6,
             'throat_loc': 1.1e-6}


class CubicSandstoneSettings(GenericSettings):
    r"""
    The following parameters are used to generate the network and geometry

    Parameters
    ----------
    sandstone : str (default = 'berea')
        The type of sandstone to model
    pore_shape : float
        The shape factor used in the Weibull distribution for pore sizes
    pore_scale : float
        The scale factor used in the Weibull distribution for pore sizes
    pore_loc : float
        The location used in the Weibull distribution for pore sizes
    throat_shape : float
        The shape factor used in the Weibull distribution for throat sizes
    throat_scale : float
        The scale factor used in the Weibull distribution for throat sizes
    throat_loc : float
        The location used in the Weibull distribution for throat sizes
    lattice_constant : float
        The lattice spacing used when creating the network
    """
    sandstone = 'berea'
    pore_shape = 1.18
    pore_scale = 6.08e-6
    pore_loc = 24.6e-6
    throat_shape = 0.536
    throat_scale = 1.904e-6
    throat_loc = 0.7e-6
    lattice_constant = 0.0001712


class CubicSandstone(Project):
    r"""
    A network representing sandstone on a Cubic lattice

    Sandstone is one of the standard materials used on geoscience studies due
    to it's importance in oil reservoir engineering as well as having well
    defined pore structure.  This class creates a Cubic network with the
    appropriate lattice spacing and connectivity, then adds a Geometry object
    with the necessary pore-scale models and prescribed parameters.

    Parameters
    ----------
    shape : array_like
        The number of pores along each direction of the domain.  All other
        aspects of this model are prescribed by the code.
    sandstone : str
        Options are 'berea' (default), 'boise', etc.
    settings : dict
        If ``sandstone`` is given as ``None``, then the user can specify their
        own parameters for the pore and throat size distribution via these
        settings.

    Notes
    -----
    The source code for this Material is relatively straight-forward, so is a
    good example starting point for creating custom materials.

    References
    ----------
    [1] ???

    Examples
    --------

    """

    def __init__(self, shape=[10, 10, 10], sandstone='berea', settings={},
                 **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(CubicSandstoneSettings())
        self.settings.update(settings)
        if sandstone is not None:
            standstone_settings = getattr(SandstoneParams, sandstone)
            self.settings.update(standstone_settings)
        pn = Cubic(shape=shape, spacing=self.settings['lattice_constant'],
                   connectivity=6, project=self)
        geom = GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
        geom['pore.seed'] = sp.rand(pn.Np)
        geom.add_model(propname='throat.seed',
                       model=mods.misc.neighbor_lookups.from_neighbor_pores,
                       pore_prop='pore.seed')

        # Pore throat and pore body characteristic dimensions follow
        # respective Weibull distribution, taking location parameters for
        # Berea 108 sample from table 5.
        geom.add_model(propname='pore.size_z',
                       model=mods.geometry.pore_size.weibull,
                       shape=self.settings['pore_shape'],
                       loc=self.settings['pore_loc'],
                       scale=self.settings['pore_scale'],
                       seeds='pore.seed')
        geom.add_model(propname='throat.size',
                       model=mods.geometry.throat_size.weibull,
                       shape=self.settings['throat_shape'],
                       loc=self.settings['throat_loc'],
                       scale=self.settings['throat_scale'],
                       seeds='throat.seed')

        # All pores in this model are of square x-section
        # All throats are of slit shape x-section
        geom['pore.size_x'] = sp.copy(geom['pore.size_z'])
        geom['pore.size_y'] = geom['pore.size_z']*1.5

        # Fetch copies of conns and coords for subsequent size calcs
        conns = pn['throat.conns']
        coords = pn['pore.coords']
        # Create Nt by 2 array of pore coords
        temp = coords[conns]
        temp = sp.absolute(temp[:, 0] - temp[:, 1])
        # Find orientation of each throat and create a label
        pn['throat.dir_x'] = temp[:, 0] > 0
        pn['throat.dir_y'] = temp[:, 1] > 0
        pn['throat.dir_z'] = temp[:, 2] > 0

        # Find width and length of each throat based on it's orientation
        # Start by initializing arrays with 0's
        geom['throat.size_x'] = 0.0
        geom['throat.size_y'] = 0.0
        geom['throat.size_z'] = 0.0
        geom['throat.length'] = 0.0
        geom['throat.width'] = 0.0
        geom['throat.height'] = 0.0

        # Start with x-directional throats
        Lc = self.settings['lattice_constant']
        Ts = pn.throats('dir_x')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_x'][Ts] = Lc - geom['pore.size_x'][conns[Ts, 0]] \
            - geom['pore.size_x'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with y-directional throats
        Ts = pn.throats('dir_y')
        geom['throat.size_z'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_y'][Ts] = Lc - geom['pore.size_y'][conns[Ts, 0]]/2 \
            - geom['pore.size_y'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_x'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_z'][Ts]

        # Start with z-directional throats
        Ts = pn.throats('dir_z')
        geom['throat.size_x'][Ts] = geom['throat.size'][Ts]
        geom['throat.size_y'][Ts] = geom['throat.size'][Ts]*6
        geom['throat.size_z'][Ts] = Lc - geom['pore.size_z'][conns[Ts, 0]]/2 \
            - geom['pore.size_z'][conns[Ts, 1]]/2
        geom['throat.length'][Ts] = geom['throat.size_z'][Ts]
        geom['throat.width'][Ts] = geom['throat.size_y'][Ts]
        geom['throat.height'][Ts] = geom['throat.size_x'][Ts]

        geom.add_model(propname='throat.area',
                       model=mods.misc.basic_math.product,
                       prop1='throat.height', prop2='throat.width')
        geom.add_model(propname='throat.volume',
                       model=mods.misc.basic_math.product,
                       prop1='throat.area',
                       prop2='throat.length')
        geom.add_model(propname='pore.volume',
                       model=mods.misc.basic_math.product,
                       prop1='pore.size_x', prop2='pore.size_y',
                       prop3='pore.size_z')