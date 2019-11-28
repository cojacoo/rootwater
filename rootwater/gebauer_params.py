"""
Gebauer Weibull Paraemeters
===========================

These are the tree-specific parameters of Gebauer et al. (2008) for the 
4-parameter Weibull function for sap velocity distribution in the sapwood. 

References
----------
Gebauer, T., Horna, V., and Leuschner, C.: Variability in radial sap flux
density patterns and sapwood area among seven co-occurring temperate 
broad-leaved tree species, Tree Physiol., 28, 1821–1830, 2008.
"""

gp = {
      'beech': {'name':'Fagus sylvatica', 'a': 2.69, 'b': 3.42, 'c': 1.00, 'd': 2.44},
      'hornbeam': {'name':'Carpinus betulus', 'a': 1.37, 'b': 5.88, 'c': 2.43, 'd': 2.79},
      'limeA': {'name':'Tilia sp. (A)', 'a': 1.62, 'b': 6.35, 'c': 2.71, 'd': 3.28},
      'limeB': {'name':'Tilia sp. (B)', 'a': 1.11, 'b': 4.52, 'c': 1.67, 'd': 1.88},
      'sycamoremaple': {'name':'Acer pseudoplatanus', 'a': 1.44, 'b': 8.98, 'c': 3.47, 'd': 3.42},
      'maple': {'name':'Acer campestre', 'a': 1.74, 'b': 4.86, 'c': 1.94, 'd': 2.50},
      'ash': {'name':'Fraxinus excelsior', 'a': 1.00, 'b': 1.44, 'c': 1.54, 'd': 0.42}
     }