"""
At the current version, which in this sense is the same as the time of fork, the used covalent radii were in nm with a
local definition. This is checked against the ASE values, which are in A.
"""
from ase import data

from nebterpolator.core import connectivity

# print(data.covalent_radii)   # list[num]=radius (A)
# print(connectivity.COVALENT_RADII) # dict[symbol]=radius (nm)
# print(data.atomic_numbers)  # dict[symbol]=num

print("{:4}{:4}{:8.4}{:8.4}{:8.4}".format("sym", "num", "this", "ase_value", "diff"))
for sym, num in data.atomic_numbers.items():

    try:
        this_package = connectivity.COVALENT_RADII[sym] * 10
        ase_value = data.covalent_radii[num]
        print("{:4}{:4}{:8.4f}{:8.4f}{:8.4f}".format(sym, num, this_package, ase_value, this_package - ase_value))
    except KeyError:
        print("{:4}{:4} gives a KeyError".format(sym, num))
