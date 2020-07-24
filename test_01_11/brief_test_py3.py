#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ase.io
import numpy as np

from nebterpolator import path_operations

# In[2]:


# In[3]:


traj = ase.io.read("traj.xyz", "::4")
ase.io.write("traj_4.xyz", traj)
len(traj)

# In[4]:


pos_nm = np.array([at.get_positions() for at in traj]) * 0.1
pos_save = pos_nm.copy()

# pos_nm.shape  # so this is (n_frames, n_atoms, 3)


# In[5]:


chemical_symbols = traj[0].get_chemical_symbols()


# chemical_symbols


# In[6]:


# result, errors = path_operations.smooth_internal(pos_nm, chemical_symbols, 7)


# In[7]:


# (result - pos_nm).max()


# In[8]:


# new_pos_a = result * 10
# traj_window7 = []
# for pos in result:
#     traj_window7.append(ase.Atoms(positions=pos*10, symbols=chemical_symbols, cell=traj[0].cell))
# ase.io.write("traj_window7.xyz", traj_window7)


# In[9]:


def calc_and_write(i):
    result, errors = path_operations.smooth_internal(pos_nm, chemical_symbols, i)

    traj_window = []
    for pos in result:
        traj_window.append(ase.Atoms(positions=pos * 10, symbols=chemical_symbols, cell=traj[0].cell))
    ase.io.write("traj_window" + str(i) + ".xyz", traj_window)


calc_and_write(41)

# In[ ]:


# raise NotImplementedError
#
#
# # In[ ]:
#
#
# # calc_and_write(15)
#
#
# # In[ ]:
#
#
# # calc_and_write(11)
#
#
# # In[ ]:
#
#
# # calc_and_write(21)
#
#
# # In[ ]:
#
#
# # calc_and_write(31)
#
#
# # # bond connectivity now
#
# # In[ ]:
#
#
# from nebterpolator.core import bond_connectivity
#
#
# # In[ ]:
#
#
# conn = bond_connectivity(pos_nm[0], chemical_symbols)
# conn
#
#
# # In[ ]:
#
#
# minima = ase.io.read("opttimised_every5th.xyz", ":")
#
#
# # In[ ]:
#
#
# pos_minima = np.array([at.get_positions() for at in minima]) * 0.1
#
#
# # In[ ]:
#
#
# conn_minima = [bond_connectivity(pos_minima[i], chemical_symbols) for i in range(len(pos_minima))]
#
#
# # In[ ]:
#
#
# conn_minima
#
#
# # In[ ]:
#
#
# # >>> G = nx.petersen_graph()
# # >>> plt.subplot(121)
# # <matplotlib.axes._subplots.AxesSubplot object at ...>
# # >>> nx.draw(G, with_labels=True, font_weight='bold')
# # >>> plt.subplot(122)
# # <matplotlib.axes._subplots.AxesSubplot object at ...>
# # >>> nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')
#
#
# # In[ ]:
#
#
# import networkx as nx
# from matplotlib import  pyplot as plt
#
#
# # In[ ]:
#
#
# graph = nx.from_edgelist(conn)
#
#
# # In[ ]:
#
#
# # nx.draw(graph)
# graph
#
#
# # In[ ]:
#
#
# colors = {"H":"tab:blue", "C":"tab:gray", "O":"tab:red"}
# node_colros = [colors[s] for s in chemical_symbols]
#
#
# # In[ ]:
#
#
# plt.subplot(121)
# nx.draw(graph, with_labels=True, font_weight='bold', node_color=node_colros)
#
#
# # In[ ]:
#
#
# minima_graphs = [nx.from_edgelist(ed) for ed in conn_minima]
#
#
# # In[ ]:
#
#
# interest = []
#
# for i in range(len(minima) - 1):
#     if not nx.algorithms.is_isomorphic(minima_graphs[i], minima_graphs[i+1]):
#         interest.append([i, i+1])
#
# print(interest)
#
#
# # In[ ]:
#
#
# nx.algorithms.is_isomorphic(minima_graphs[0], minima_graphs[-1])
#
#
# # In[ ]:
#
#
# conn_minima[-1]
#
#
# # In[ ]:
#
#
# plt.subplot(111)
# nx.draw(minima_graphs[-1], with_labels=True, font_weight='bold', node_color=node_colros)
#
#
# # In[ ]:
#
#
#
#
