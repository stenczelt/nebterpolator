import ase.io
import networkx as nx
import numpy as np

from nebterpolator.core import bond_connectivity
# import numpy as np
from nebterpolator.path_operations import smooth_internal, smooth_cartesian

# import nebterpolator
# from nebterpolator import path_operations
# from matplotlib import  pyplot as plt
# from matplotlib import gridspec


species_colors = {"H": "k", "C": "tab:gray", "O": "tab:red"}


# node_colros = [species_colors[s] for s in chemical_symbols]

def atoms_to_graph(atoms: ase.Atoms):
    at_graph = nx.Graph()

    for i, sym in enumerate(atoms.get_chemical_symbols()):
        at_graph.add_node(i, symbol=sym)

    bonds = bond_connectivity(atoms.get_positions() * 0.1, atoms.get_chemical_symbols())
    at_graph.add_edges_from(bonds)

    return at_graph


def get_molecules_in_graph(at_graph: nx.Graph):
    subsets = list(nx.connected_components(at_graph))
    mol_graphs = []
    for node_set in subsets:
        mol_graphs.append(at_graph.subgraph(node_set))

    return mol_graphs


def get_connectivity_change(frames):
    """Finds where changes happened in the connectivity of a trajectory.

    Connectivity change is understood as the molecular graphs not being isomorphic.
    nb. distances in the connectivity calculator need to be in nanometers, not A.

    Parameters
    ----------
    frames : list(Ase.Atoms)
        trajectory, atom indices need to match

    Returns
    -------
    interest : list of list, shape(n_changes, 2)
        The pairs of consecutive frames where there was a change in connectivity


    """
    graphs_traj = [atoms_to_graph(at) for at in frames]

    interest = []
    for i in range(len(frames) - 1):
        this_change = not nx.algorithms.is_isomorphic(graphs_traj[i], graphs_traj[i + 1])
        if this_change:
            interest.append([i, i + 1])

    return interest


def smooth_traj(filename_in, filename_out, smoothing_width=None, xyz_smoothing_strength=2.0, w_morse=0.0):
    """

    Parameters
    ----------
    filename_in
    filename_out
    smoothing_width: float
        these two parameters are adjustable, and depend on the length of the traj
        cutoff period for the internal coordinate smoother. motions with a shorter
        period than this (higher frequency) will get filtered out
    xyz_smoothing_strength:  float, default=2.0
        the spline smoothing factor used for the cartesian smoothing step, that
        runs after the internal coordinates smoother. The point of this is ONLY
        to correct for "jitters" in the xyz coordinates that are introduced by
        imperfections in the redundant internal coordinate -> xyz coordinate
        step, which runs after smoothing in internal coordinates
    morse

    Returns
    -------


    Notes
    -----
    Adapted from the Nebterpolate.py script. The comments for the arguments and here are
    from the authors of that script.

    ## How smoothing works:
    transform into redundant internal coordinates, apply a fourier based
    smoothing, and then transform back to cartesian.
    the internal -> cartesian bit is the hard step, since there's no
    guarantee that a set of cartesian coordinates even exist that satisfy
    the redundant internal coordinates, after smoothing.

    ## inversion:
    we're using a levenberg-marquardt optimizer to find the "most consistent"
    cartesian coordinates

    currently, the choice of what internal coordinates to use is buried
    a little in the code, in the function path_operations.union_connectivity
    basically, we're using ALL pairwise distances, all of the angles between
    sets of three atoms, a-b-c, that actually get "bonded" during the
    trajectory, and all of the dihedral angles between sets of 4 atoms,
    a-b-c-d, that actually get "bonded" during the trajectory.

    """

    nm_in_angstrom = 0.1

    # read the input file
    frames = ase.io.read(filename_in, ":")
    xyzlist = np.array([fr.get_positions() for fr in frames])
    xyzlist *= nm_in_angstrom
    atom_names = frames[0].get_chemical_symbols()

    if xyzlist.shape[1] < 4:
        ValueError("Interpolator cannot handle less than four atoms.")

    # If smoothing width not provided, default to using the trajectory length
    if smoothing_width is None:
        smoothing_width = len(xyzlist)
        smoothing_width += smoothing_width % 2
        smoothing_width -= 1
    if smoothing_width % 2 != 1:
        raise RuntimeError("Smoothing width must be an odd number")

    # smoothing in internal coordinates
    smoothed, errors = smooth_internal(xyzlist, atom_names, width=smoothing_width, bond_width=smoothing_width,
                                       angle_width=smoothing_width, dihedral_width=smoothing_width, w_morse=w_morse)

    # apply a bit of spline smoothing in cartesian coordinates to
    # correct for jitters
    jitter_free = smooth_cartesian(smoothed,
                                   strength=xyz_smoothing_strength,
                                   weights=1.0 / errors)

    # reset the positions in the file
    final_pos = jitter_free / nm_in_angstrom
    for i, at in enumerate(frames):
        at.set_positions(final_pos[i])

    # write the output file
    print('Saving output to', filename_out)
    ase.io.write(filename_out, frames)
