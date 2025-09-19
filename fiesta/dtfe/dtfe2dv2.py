import numpy as np
from scipy.spatial import Delaunay as scDelaunay

import shift

from .. import boundary
from .. import coords
from .. import src


class Delaunay2Dv2:
    
    
    def __init__(self):
        self.x = None
        self.y = None
        self.f = None
        self.mass = None
        self.points = None
        self.point_type = None
        self.point_dens = None
        self.boundary = None
        self.cross_boundary = None
        self.delaunay = None
        self.delaunay_area = None
        self.simplices = None
        self.simptypes = None
        self.circumcircles = None
        self.angle_total = None
        self.dtfe_mode = None
        self.x0 = None
        self.y0 = None
        self.f0 = None
        self.delf0 = None
        self.verbose = True


    def setup(self, x, y, boundary, mass=None, f=None):
        """
        Assign particles and field values.
        
        Parameters
        ----------
        x : array
            X-axis positions for points.
        y : array
            Y-axis positions for points.
        boundary : list
            List of boundary values given as [xmin, xmax, ymin, ymax].
        mass : array, optional
            Mass weights.
        f : array, optional
            Field values.
        """
        self.x = x
        self.y = y
        self.boundary = boundary
        if mass is None:
            self.mass = np.ones(len(x))
        else:
            self.mass = mass
        self.f = f
        self.points = coords.coord2points([self.x, self.y])

    
    def _check_delaunay(self):
        mask = np.zeros(len(self.points))
        idx = np.unique(self.simplices.flatten())
        mask[idx] = 1.
        cond = np.where(mask != 1.)[0]
        if len(cond) != 0:
            if self.verbose:
                print('Delaunay computation merged some points due to floating type precision, weights for particles in the DTFE are being adjusted to account for this.')
            for idx in cond:
                idxnearest = np.argmin((self.x-self.x[idx])**2. + (self.y-self.y[idx])**2.)
                self.w[idxnearest] += self.w[idx]
                self.w[idx] = 0.
    
    
    def construct(self):
        """
        Construct Delaunay triangulation.
        """
        self.delaunay = scDelaunay(self.points)
        self.simplices = self.delaunay.simplices
        self._check_delaunay()
        

    def _get_simplex_circumcircle(self, simplex):
        """
        Returns the circumcirle center and radius for a given simplex.
        
        Parameters
        ----------
        simplex : array
            An array of length 3 containing the index for the points in the simplex.
            
        Returns
        -------
        Circumcircle : array
            An array of length 3, where the first two values are the x-axis and y-axis circumcircle center
            and the last term is the circumcircle radius.
        """
        vertices = self.points[simplex]
        x0, y0 = vertices[0,0], vertices[0,1]
        x1, y1 = vertices[1,0], vertices[1,1]
        x2, y2 = vertices[2,0], vertices[2,1]
        A0 = x0**2 + y0**2
        A1 = x1**2 + y1**2
        A2 = x2**2 + y2**2
        D = 2*np.linalg.det(np.array([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]]))
        
        circumcircle = np.zeros(3)
        circumcircle[0] = np.linalg.det(np.array([[A0, y0, 1], [A1, y1, 1], [A2, y2, 1]]))/D
        circumcircle[1] = np.linalg.det(np.array([[x0, A0, 1], [x1, A1, 1], [x2, A2, 1]]))/D
        circumcircle[2] = np.sqrt((x0-circumcircle[0])**2 + (y0-circumcircle[1])**2)
        
        return circumcircle
    
    
    def get_circumcircles(self):
        """
        Computes the circumcirles centers and radius for all simplices.
        """
        self.circumcircles = np.array([self._get_simplex_circumcircle(simplex) for simplex in self.simplices])
        
    
    def _get_simplex_angles(self, simplex):
        """
        Returns the simplex angles.
        
        Parameters
        ----------
        simplex : array
            An array of length 3 containing the index for the points in the simplex.
            
        Returns
        -------
        angles : array
            The angles subtended by the simplex from the corresponding point in the simplex.
        """
        vertices = self.points[simplex]
        angles = np.zeros(3)
        x0, y0 = vertices[0,0], vertices[0,1]
        x1, y1 = vertices[1,0], vertices[1,1]
        x2, y2 = vertices[2,0], vertices[2,1]
        u1x = x1 - x0
        u1y = y1 - y0
        u2x = x2 - x0
        u2y = y2 - y0
        angles[0] = np.arctan2(u1x*u2y - u1y*u2x, u1x*u2x + u1y*u2y)
        u0x = x0 - x1
        u0y = y0 - y1
        u2x = x2 - x1
        u2y = y2 - y1
        angles[1] = np.arctan2(u0x*u2y - u0y*u2x, u0x*u2x + u0y*u2y)
        u0x = x0 - x2
        u0y = y0 - y2
        u1x = x1 - x2
        u1y = y1 - y2
        angles[2] = np.arctan2(u0x*u1y - u0y*u1x, u0x*u1x + u0y*u1y)
        return angles
    
    
    def _add_angles_2_point(self, simplex):
        """
        Adds the angles from a simplex to the total angle subtended by simplices attached to each points.
        
        Parameters
        ----------
        simplex : array
            An array of length 3 containing the index for the points in the simplex.
        """
        angles = self._get_simplex_angles(simplex)
        self.angle_total[simplex] += abs(angles)
    
    
    def get_angle_total(self):
        """
        Determines the total angle subtended by simplices attached to each points.
        """
        self.angle_total = np.zeros(len(self.points))
        [self._add_angles_2_point(simplex) for simplex in self.simplices]
        
    
    def _get_boundary_simplex(self, idx):
        """
        Determines whether a simplex circumcircle crosses boundaries.
        """
        xc, yc, r = self.circumcircles[idx][0], self.circumcircles[idx][1], self.circumcircles[idx][2]
        if xc-r <= self.boundary[0]:
            self.cross_boundary[idx][0] = 1.
        if xc+r >= self.boundary[1]:
            self.cross_boundary[idx][1] = 1.
        if yc-r <= self.boundary[2]:
            self.cross_boundary[idx][2] = 1.
        if yc+r >= self.boundary[3]:
            self.cross_boundary[idx][3] = 1.
    
        
    def get_boundary_simplex(self):
        """
        Determines whether a simplex circumcircle crosses boundaries.
        """
        self.cross_boundary = np.zeros((len(self.simplices), 4))
        [self._get_boundary_simplex(idx) for idx in range(len(self.simplices))]
    
    
    def determine_delaunay_boundary(self):
        """
        Determine Delaunay boundary simplices and points.
        """
        # TODO add more comments
        self.simptypes = np.ones(len(self.simplices))
        cond = np.where(np.sum(self.cross_boundary, axis=1) != 0.)[0]
        self.simptypes[cond] = -1. # Circumcircle center is completely outside the boundary
        self.point_type = np.ones(len(self.points))
        # first column is tied to the simplex type and the second to whether a particle's angles subtend the full circle.
        cond = np.where(np.abs(self.angle_total - 2*np.pi) > 1e-3)[0]
        self.point_type[cond] = -1. 
        cond = np.where(self.simptypes == -1.)[0]
        idxs = np.unique(self.simplices[cond].flatten())
        self.point_type[idxs] = -1.
        cond = np.where(self.simptypes == 1.)[0]
        cond = np.where(self.simptypes == 1.)[0]
        simplices = self.simplices[cond]
        point_type0 = self.point_type[simplices[:,0]]
        point_type1 = self.point_type[simplices[:,1]]
        point_type2 = self.point_type[simplices[:,2]]
        cond2 = np.where((point_type0 == -1.) | (point_type1 == -1.) | (point_type2 == -1.))[0]
        self.simptypes[cond[cond2]] = 0.
        cond = np.where(self.simptypes == 0.)[0]
        idxs = np.unique(self.simplices[cond].flatten())
        cond = np.where(self.point_type[idxs] != -1)[0]
        self.point_type[idxs[cond]] = 0.

    def run(self):
        """
        Run Delaunay tesselation and compute boundary information.
        """
        self.construct()
        self.get_circumcircles()
        self.get_angle_total()
        self.get_boundary_simplex()
        self.determine_delaunay_boundary()

    
    def get_area(self):
        """Calculates the area of the delaunay simplex."""
        x, y = self.points[:, 0], self.points[:, 1]
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        self.delaunay_area = src.delaunay_area_2d(x, y, del_vert0, del_vert1=del_vert1, del_vert2=del_vert2)
    
    
    def get_dens(self):
        """Calculates the density of each point in the delaunay tessellation."""
        if self.delaunay_area is None:
            self.get_area()
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        point_area = src.sum_delaunay_area_4_points_2d(self.delaunay_area, del_vert0, del_vert1, del_vert2, len(self.points))
        self.point_dens = self.mass/point_area
        
    
    def set_field(self, f):
        """Sets the field values of the input points.

        Parameters
        ----------
        f : array
            Field values.
        """
        if f is not None:
            lenf = len(f)
            assert lenf == len(self.points), "f must be equal to input points."
            self.f = f
        x, y, f = self.points[:, 0], self.points[:, 1], self.f
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        self.x0 = x[del_vert0]
        self.y0 = y[del_vert0]
        self.f0 = f[del_vert0]
        self.delf0 = src.get_delf0_2d(x, y, f, del_vert0, del_vert1, del_vert2)
        self.dtfe_mode = 'field'
    
    
    def set_field2dens(self):
        self.set_field(self.point_dens)
        self.dtfe_mode = 'density'
    
        
    def find_simplex(self, x, y):
        """Find the simplex the coordinates lie within."""
        points = coords.xy2points(x, y)
        simplices = self.delaunay.find_simplex(points)
        return simplices
    

    def estimate(self, x, y):
        """Estimates a field from the Delaunay tesselation.

        Parameters
        ----------
        x : array
            X-coordinate for the field estimation.
        y : array
            Y-coordinate for the field estimation.

        Returns
        -------
        f_est : array
            Estimates of the field
        """
        simplices_idx = self.find_simplex(x, y)
        f_est = np.nan*np.ones(len(x))
        # if self.dtfe_mode == 'field':
        #     cond = np.where((simplices_idx != -1) & (self.simptypes[simplices_idx] != -1.))[0]
        # if self.dtfe_mode == 'density':
        cond = np.where((simplices_idx != -1) & (self.simptypes[simplices_idx] == 1.))[0]
        f_est[cond] = src.delaunay_estimate_2d(simplices_idx[cond], x[cond], y[cond], self.x0, self.y0, self.f0, self.delf0)
        return f_est


    def get_border(self):
        """
        Returns border points and simplices.
        """
        border = {}
        keep_idx = -np.ones(len(self.x), dtype='int')
        cond = np.where(self.point_type != 1)[0]
        keep_idx[cond] = np.arange(len(cond))
        border['x'] = self.x[cond]
        border['y'] = self.y[cond]
        border['mass'] = self.mass[cond]
        border['ptype'] = self.point_type[cond]
        if self.f is None:
            border['f'] = None
        else:
            border['f'] = self.f[cond]
        border['dtfe_mode'] = self.dtfe_mode
        cond = np.where(self.simptypes != 1)[0]
        border['simplices'] = keep_idx[self.simplices[cond]]
        border['simptypes'] = self.simptypes[cond]
        return border


    def clean(self):
        self.__init__()


class DelaunayJoiner2D:
    
    
    def __init__(self):
        """
        """
        self.x = []
        self.y = []
        self.f = []
        self.mass = []
        self.Np = 0
        self.Nb = 0
        self.pID = []
        self.ptypes = []
        self.simplices = []
        self.simptypes = []
        self.dtfe_mode = None
        self.dtfe = None
        
    
    def add_border(self, border):
        """
        """
        if self.dtfe_mode is None:
            self.dtfe_mode = border['dtfe_mode']
        else:
            assert self.dtfe_mode == border['dtfe_mode'], 'Border DTFE must match preassigned border DTFEs.'
        self.x.append(border['x'])
        self.y.append(border['y'])
        self.f.append(border['f'])
        self.mass.append(border['mass'])
        self.pID.append(self.Nb*np.ones(len(border['x'])))
        self.Nb += 1
        self.ptypes.append(border['ptype'])
        idx2globalidx = np.arange(len(border['x'])) + self.Np
        self.Np += len(border['x'])
        self.simplices.append(idx2globalidx[border['simplices']])
        self.simptypes.append(border['simptypes'])
        
        
    
    def _concatenate(self):
        """
        """
        self.x = np.concatenate(self.x)
        self.y = np.concatenate(self.y)
        self.f = np.concatenate(self.f)
        self.mass = np.concatenate(self.mass)
        self.pID = np.concatenate(self.pID)
        self.ptypes = np.concatenate(self.ptypes)
        self.simplices = np.concatenate(self.simplices)
        self.simptypes = np.concatenate(self.simptypes)
        
        
    def run(self, boundary):
        """
        """
        self._concatenate()
        self.dtfe = fiesta.dtfe.Delaunay2Dv2()
        self.dtfe.setup(self.x, self.y, boundary, mass=self.mass)
        self.dtfe.run()
        if self.dtfe_mode == 'density':
            self.dtfe.get_dens()
            cond = np.where(self.ptypes == 0.)[0]
            self.dtfe.point_dens[cond] = self.f[cond]
            self.dtfe.set_field2dens()
        else:
            self.dtfe.set_field(self.f)
            
        cond = np.where(
            (self.pID[self.dtfe.simplices[:,0]] == self.pID[self.dtfe.simplices[:,1]]) & 
            (self.pID[self.dtfe.simplices[:,1]] == self.pID[self.dtfe.simplices[:,2]]) &
            (self.ptypes[self.dtfe.simplices[:,0]] == 0.) & 
            (self.ptypes[self.dtfe.simplices[:,1]] == 0.) & 
            (self.ptypes[self.dtfe.simplices[:,2]] == 0.)
        )[0]
        self.dtfe.simptypes[cond] = -1.
    
    
    def estimate(self, x, y):
        """"""
        return self.dtfe.estimate(x, y)
    
    
    def get_border(self):
        """"""
        return self.dtfe.get_border()
    
    
    def clean(self):
        """"""
        return self.__init__()