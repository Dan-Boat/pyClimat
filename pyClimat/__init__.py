# Initialise all packages for just importing Package 
try:
    
    from .plots import *
    from .plot_utils import *
    from .analysis import *
    from .data import *
except:
    from plots import *
    from plot_utils import *
    from analysis import *
    from data import *
 