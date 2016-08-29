import matplotlib.path
if hasattr(matplotlib.path.Path, 'contains_points'):
    from matplotlib.path import Path
else:
    from .path import Path
