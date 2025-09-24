import numpy as np
from scipy.spatial import Delaunay as scDelaunay

import shift

from .. import boundary
from .. import coords
from .. import src
from .. import maths

from typing import List, Optional, Union


class Delaunay3D:
    
    
    def __init__(self) -> None:
        """
        Initialises the class function.
        """
        self.x = None
        self.y = None
        self.z = None
        self.f = None
        self.mass = None
        self.points = None
        self.point_type = None
        self.point_dens = None
        self.boundary = None
        self.cross_boundary = None
        self.delaunay = None
        self.delaunay_volume = None
        self.simplices = None
        self.simptypes = None
        self.circumcircles = None
        self.solid_angle_total = None
        self.dtfe_mode = None
        self.x0 = None
        self.y0 = None
        self.z0 = None
        self.f0 = None
        self.delf0 = None
        self.verbose = True


    def setup(
            self, x: np.ndarray, y: np.ndarray, z: np.ndarray, boundary: List[float], mass: Optional[np.ndarray] = None,
        ) -> None:
        """
        Assign particles and field values.
        
        Parameters
        ----------
        x : array
            X-axis positions for points.
        y : array
            Y-axis positions for points.
        z : array
            Z-axis positions for points.
        boundary : list
            List of boundary values given as [xmin, xmax, ymin, ymax, zmin, zmax].
        mass : array, optional
            Mass weights.
        """
        self.x = x
        self.y = y
        self.z = z
        self.boundary = boundary
        if mass is None:
            self.mass = np.ones(len(x))
        else:
            self.mass = mass
        self.points = coords.coord2points([self.x, self.y, self.z])

    
    def _check_delaunay(self) -> None:
        """
        This will test whether the delaunay tesselation is fully connected, and test and correct
        for merging of points due to floating point precision.
        """
        mask = np.zeros(len(self.points))
        idx = np.unique(self.simplices.flatten())
        mask[idx] = 1.
        cond = np.where(mask != 1.)[0]
        if len(cond) != 0:
            if self.verbose:
                print('Delaunay computation merged some points due to floating type precision, weights for particles in the DTFE are being adjusted to account for this.')
            for idx in cond:
                idxnearest = np.argmin((self.x-self.x[idx])**2. + (self.y-self.y[idx])**2. + (self.z-self.z[idx])**2.)
                self.mass[idxnearest] += self.mass[idx]
                self.mass[idx] = 0.
    
    
    def construct(self) -> None:
        """
        Construct Delaunay triangulation.
        """
        self.delaunay = scDelaunay(self.points)
        self.simplices = self.delaunay.simplices
        self._check_delaunay()
    

    def get_circumcircles(self) -> None:
        """
        Computes the circumcirles centers and radius for all simplices.
        """
        x0, y0, z0 = self.x[self.simplices[:,0]], self.y[self.simplices[:,0]], self.z[self.simplices[:,0]]
        x1, y1, z1 = self.x[self.simplices[:,1]], self.y[self.simplices[:,1]], self.z[self.simplices[:,1]]
        x2, y2, z2 = self.x[self.simplices[:,2]], self.y[self.simplices[:,2]], self.z[self.simplices[:,2]]
        x3, y3, z3 = self.x[self.simplices[:,3]], self.y[self.simplices[:,3]], self.z[self.simplices[:,3]]
        one = np.ones(len(x0))
        D = 2*maths.det4by4(x0, y0, z0, one, x1, y1, z1, one, x2, y2, z2, one, x3, y3, z3, one)
        cond = np.where(D != 0.)[0]
        A0 = x0**2 + y0**2 + z0**2
        A1 = x1**2 + y1**2 + z1**2
        A2 = x2**2 + y2**2 + z2**2
        A3 = x3**2 + y3**2 + z3**2
        self.circumcircles = np.nan*np.ones((len(x0), 4))
        self.circumcircles[cond,0] = maths.det4by4(
            A0[cond], y0[cond], z0[cond], one[cond], A1[cond], y1[cond], z1[cond], one[cond], 
            A2[cond], y2[cond], z2[cond], one[cond], A3[cond], y3[cond], z3[cond], one[cond]
        )/D[cond]
        self.circumcircles[cond,1] = maths.det4by4(
            x0[cond], A0[cond], z0[cond], one[cond], x1[cond], A1[cond], z1[cond], one[cond], 
            x2[cond], A2[cond], z2[cond], one[cond], x3[cond], A3[cond], z3[cond], one[cond]
        )/D[cond]
        self.circumcircles[cond,2] = maths.det4by4(
            x0[cond], y0[cond], A0[cond], one[cond], x1[cond], y1[cond], A1[cond], one[cond], 
            x2[cond], y2[cond], A2[cond], one[cond], x3[cond], y3[cond], A3[cond], one[cond]
        )/D[cond]
        self.circumcircles[cond,3] = np.sqrt((x0[cond]-self.circumcircles[cond,0])**2 + (y0[cond]-self.circumcircles[cond,1])**2 + (z0[cond]-self.circumcircles[cond,2])**2)


    def get_solid_angle(self, 
            px0: Union[float, np.ndarray], py0: Union[float, np.ndarray], pz0: Union[float, np.ndarray], 
            x1: Union[float, np.ndarray], y1: Union[float, np.ndarray], z1: Union[float, np.ndarray], 
            x2: Union[float, np.ndarray], y2: Union[float, np.ndarray], z2: Union[float, np.ndarray], 
            x3: Union[float, np.ndarray], y3: Union[float, np.ndarray], z3: Union[float, np.ndarray]
        ) -> np.ndarray:
        """
        Computes the solid angle subtended from a point p and three points 1, 2 and 3.

        Parameters
        ----------
        px0: float or array
            X value for point we want to compute the angle from.
        py0: float or array
            Y value for point we want to compute the angle from.
        pz0: float or array
            Z value for point we want to compute the angle from.
        x1 : float or array
            Point 1's x value.
        y1 : float or array
            Point 1's y value.
        z1 : float or array
            Point 1's z value.
        x2 : float or array
            Point 2's x value.
        y2 : float or array
            Point 2's y value.
        z2 : float or array
            Point 2's z value.
        x3 : float or array
            Point 3's x value.
        y3 : float or array
            Point 3's y value.
        z3 : float or array
            Point 3's z value.
            
        Returns
        -------
        solid_angle : array
            Solid angle subtended from a point p and two points 1, 2 and 3.
        """
        ax = x1 - px0
        ay = y1 - py0
        az = z1 - pz0
        bx = x2 - px0
        by = y2 - py0
        bz = z2 - pz0
        cx = x3 - px0
        cy = y3 - py0
        cz = z3 - pz0

        # numerator: a dot (b cross c)
        # b cross c = (bycz - bzcy, bzcx - bxcz, bxcy - bycx)
        # a dot (b cross c) = ax(bycz - bzcy) + ay(bzcx - bxcz) + az(bxcy - bycx)
        num = ax*(by*cz - bz*cy) + ay*(bz*cx - bx*cz) + az*(bx*cy - by*cx)

        # Norms
        na = np.sqrt(ax**2 + ay**2 + az**2)
        nb = np.sqrt(bx**2 + by**2 + bz**2)
        nc = np.sqrt(cx**2 + cy**2 + cz**2)

        # denominator: na*nb*nc + (a dot b)*nc + (b dot c)*na + (c dot a)*nb
        den = na*nb*nc + (ax*bx+ay*by+az*bz)*nc + (bx*cx+by*cy+bz*cz)*na + (cx*ax+cy*ay+cz*az)*nb

        solid_angle = 2*abs(np.arctan2(num, den))
        return solid_angle


    def get_solid_angle_total(self) -> None:
        """
        Determines the total solid angle subtended by simplices attached to each points.
        """
        id0 = self.simplices[:,0]
        id1 = self.simplices[:,1]
        id2 = self.simplices[:,2]
        id3 = self.simplices[:,3]
        x0 = self.x[id0]
        y0 = self.y[id0]
        z0 = self.z[id0]
        x1 = self.x[id1]
        y1 = self.y[id1]
        z1 = self.z[id1]
        x2 = self.x[id2]
        y2 = self.y[id2]
        z2 = self.z[id2]
        x3 = self.x[id3]
        y3 = self.y[id3]
        z3 = self.z[id3]
        self.solid_angle_total = np.zeros(len(self.points))
        solid_angle = self.get_solid_angle(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3)
        np.add.at(self.solid_angle_total, id0, solid_angle)
        solid_angle = self.get_solid_angle(x1, y1, z1, x2, y2, z2, x3, y3, z3, x0, y0, z0)
        np.add.at(self.solid_angle_total, id1, solid_angle)
        solid_angle = self.get_solid_angle(x2, y2, z2, x3, y3, z3, x0, y0, z0, x1, y1, z1)
        np.add.at(self.solid_angle_total, id2, solid_angle)
        solid_angle = self.get_solid_angle(x3, y3, z3, x0, y0, z0, x1, y1, z1, x2, y2, z2)
        np.add.at(self.solid_angle_total, id3, solid_angle)


    def get_boundary_simplex(self) -> None:
        """
        Determines whether a simplex circumcircle crosses boundaries.
        """
        xc, yc, zc, r = self.circumcircles[:,0], self.circumcircles[:,1], self.circumcircles[:,2], self.circumcircles[:,3]
        self.cross_boundary = np.zeros((len(self.simplices), 6))
        cond = np.where(xc-r <= self.boundary[0])[0]
        self.cross_boundary[cond, 0] = 1.
        cond = np.where(xc+r >= self.boundary[1])[0]
        self.cross_boundary[cond, 1] = 1.
        cond = np.where(yc-r <= self.boundary[2])[0]
        self.cross_boundary[cond, 2] = 1.
        cond = np.where(yc+r >= self.boundary[3])[0]
        self.cross_boundary[cond, 3] = 1.
        cond = np.where(zc-r <= self.boundary[4])[0]
        self.cross_boundary[cond, 4] = 1.
        cond = np.where(zc+r >= self.boundary[5])[0]
        self.cross_boundary[cond, 5] = 1.
        
    
    def determine_delaunay_boundary(self) -> None:
        """
        Determine Delaunay boundary simplices and points.
        """
        # Construct simplice type array
        self.simptypes = np.ones(len(self.simplices))
        # find simplices that cross a boundary and/or have circumcircles that are undefined (i.e. NaNs due to determinant=0.)
        cond = np.where((np.sum(self.cross_boundary, axis=1) != 0.) | (np.isfinite(self.circumcircles[:,0]) == False))[0]
        self.simptypes[cond] = -1. # Circumcircle center is completely outside the boundary
        # Define point type array
        self.point_type = np.ones(len(self.points))
        # Determine points that lie at the border of the DTFE, by looking at the angle subtended by each simplex for each point
        cond = np.where(np.abs(self.solid_angle_total - 4*np.pi) > 1e-3)[0]
        self.point_type[cond] = -1. 
        # Assign points in simplex type -1 to have point type of -1 too.
        cond = np.where(self.simptypes == -1.)[0]
        idxs = np.unique(self.simplices[cond].flatten())
        self.point_type[idxs] = -1.
        # For simplices with any simplex of type 1 change the type to type 0 if one point in the simplex is of type -1.
        cond = np.where(self.simptypes == 1.)[0]
        simplices = self.simplices[cond]
        point_type0 = self.point_type[simplices[:,0]]
        point_type1 = self.point_type[simplices[:,1]]
        point_type2 = self.point_type[simplices[:,2]]
        point_type3 = self.point_type[simplices[:,3]]
        cond2 = np.where((point_type0 == -1.) | (point_type1 == -1.) | (point_type2 == -1.) | (point_type3 == -1.))[0]
        self.simptypes[cond[cond2]] = 0.
        # For points in simplices of type 0, assign all points that are not type -1 to type 0.
        cond = np.where(self.simptypes == 0.)[0]
        idxs = np.unique(self.simplices[cond].flatten())
        cond = np.where(self.point_type[idxs] != -1)[0]
        self.point_type[idxs[cond]] = 0.


    def run(self) -> None:
        """
        Run Delaunay tesselation and compute boundary information.
        """
        self.construct()
        self.get_circumcircles()
        self.get_solid_angle_total()
        self.get_boundary_simplex()
        self.determine_delaunay_boundary()

    
    def get_volume(self) -> None:
        """
        Calculates the volume of the delaunay simplex.
        """
        x, y, z = self.x, self.y, self.z
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        del_vert3 = self.simplices[:, 3]
        self.delaunay_volume = src.delaunay_volume_3d(x, y, z, del_vert0, del_vert1, del_vert2, del_vert3)
    
    
    def get_dens(self) -> None:
        """
        Calculates the density of each point in the delaunay tessellation.
        """
        if self.delaunay_volume is None:
            self.get_volume()
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        del_vert3 = self.simplices[:, 3]
        point_volume = src.sum_delaunay_vol_4_points_3d(self.delaunay_volume, del_vert0, del_vert1, del_vert2, del_vert3, len(self.points))
        self.point_dens = self.mass/point_volume
        
    
    def set_field(self, f: np.ndarray) -> None:
        """
        Sets the field values of the input points.

        Parameters
        ----------
        f : array
            Field values.
        """
        if f is not None:
            lenf = len(f)
            assert lenf == len(self.points), "f must be equal to input points."
            self.f = f
        x, y, z, f = self.points[:, 0], self.points[:, 1], self.points[:,2], self.f
        del_vert0 = self.simplices[:, 0]
        del_vert1 = self.simplices[:, 1]
        del_vert2 = self.simplices[:, 2]
        del_vert3 = self.simplices[:, 3]
        self.x0 = x[del_vert0]
        self.y0 = y[del_vert0]
        self.z0 = z[del_vert0]
        self.f0 = f[del_vert0]
        self.delf0 = src.get_delf0_3d(x, y, z, f, del_vert0, del_vert1, del_vert2, del_vert3)
        self.dtfe_mode = 'field'
    
    
    def set_field2dens(self) -> None:
        """
        Sets the field values to the density.
        """
        self.set_field(self.point_dens)
        self.dtfe_mode = 'density'
    
        
    def find_simplex(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray], z: Union[float, np.ndarray]) -> np.ndarray:
        """
        Find the simplex the coordinates lie within.

        Parameters
        ----------
        x : array
            X-coordinate for the field estimation.
        y : array
            Y-coordinate for the field estimation.
        z : array
            Z-coordinate for the field estimation.
        
        Returns
        -------
        simplices : array
            Simplex ID that the point lies in, -1 if point lies outside of all simplices.
        """
        points = coords.xyz2points(x, y, z)
        simplices = self.delaunay.find_simplex(points)
        return simplices
    

    def estimate(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray], z: Union[float, np.ndarray], allow_border: bool = False) -> Union[float, np.ndarray]:
        """
        Estimates a field from the Delaunay tesselation.

        Parameters
        ----------
        x : array
            X-coordinate for the field estimation.
        y : array
            Y-coordinate for the field estimation.
        z : array
            Z-coordinate for the field estimation.
        allow_border : bool, optional
            Allow border interpolation for non-density related interpolation.
        
        Returns
        -------
        f_est : array
            Estimates of the field
        """
        simplices_idx = self.find_simplex(x, y, z)
        f_est = np.nan*np.ones(len(x))
        if allow_border and self.dtfe_mode == 'field':
            cond = np.where((simplices_idx != -1) & (self.simptypes[simplices_idx] != -1.))[0]
        else:
            cond = np.where((simplices_idx != -1) & (self.simptypes[simplices_idx] == 1.))[0]
        f_est[cond] = src.delaunay_estimate_3d(simplices_idx[cond], x[cond], y[cond], z[cond], self.x0, self.y0, self.z0, self.f0, self.delf0)
        return f_est


    def get_border(self) -> dict:
        """
        Returns border points and simplices.

        Returns
        -------
        border : dict
            Dictionary object contained border points and simplices from previously computed Delaunay tesselation.
        """
        border = {}
        keep_idx = -np.ones(len(self.x), dtype='int')
        cond = np.where(self.point_type != 1)[0]
        keep_idx[cond] = np.arange(len(cond))
        border['x'] = self.x[cond]
        border['y'] = self.y[cond]
        border['z'] = self.z[cond]
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


    def clean(self) -> None:
        """
        Reinitialise the class.
        """
        self.__init__()



class DelaunayMerger3D:
    
    
    def __init__(self) -> None:
        """
        Initiates class.
        """
        self.x = []
        self.y = []
        self.z = []
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
        
    
    def add_border(self, border: dict) -> None:
        """
        Add border points and simplices from an already computed Delaunay tesselation.

        Parameters
        ----------
        border : dict
            Dictionary object contained border points and simplices from previously computed Delaunay tesselation.
        """
        if self.dtfe_mode is None:
            self.dtfe_mode = border['dtfe_mode']
        else:
            assert self.dtfe_mode == border['dtfe_mode'], 'Border DTFE must match preassigned border DTFEs.'
        self.x.append(border['x'])
        self.y.append(border['y'])
        self.z.append(border['z'])
        self.f.append(border['f'])
        self.mass.append(border['mass'])
        self.pID.append(self.Nb*np.ones(len(border['x'])))
        self.Nb += 1
        self.ptypes.append(border['ptype'])
        idx2globalidx = np.arange(len(border['x'])) + self.Np
        self.Np += len(border['x'])
        self.simplices.append(idx2globalidx[border['simplices']])
        self.simptypes.append(border['simptypes'])
    
    
    def _concatenate(self) -> None:
        """
        Concatenates all border points and simplices.
        """
        self.x = np.concatenate(self.x)
        self.y = np.concatenate(self.y)
        self.z = np.concatenate(self.z)
        self.f = np.concatenate(self.f)
        self.mass = np.concatenate(self.mass)
        self.pID = np.concatenate(self.pID)
        self.ptypes = np.concatenate(self.ptypes)
        self.simplices = np.concatenate(self.simplices)
        self.simptypes = np.concatenate(self.simptypes)

    
    def _filter(self, boundary) -> None:
        """
        Filters out points and simplices that go beyond the boundary

        Parameters
        ----------
        Boundary : list
            A list containing the minimum and maximum x values, given in this order: [xmin, xmax, ymin, ymax].
        """
        mask4part = np.zeros(len(self.x))
        keep_idx = -np.ones(len(self.x), dtype='int')
        cond = np.where(
            (self.x >= boundary[0]) & (self.x <= boundary[1]) & 
            (self.y >= boundary[2]) & (self.y <= boundary[3]) & 
            (self.z >= boundary[4]) & (self.z <= boundary[5]))[0]
        mask4part[cond] = 1.
        keep_idx[cond] = np.arange(len(cond))
        self.x = self.x[cond]
        self.y = self.y[cond]
        self.z = self.z[cond]
        self.mass = self.mass[cond]
        self.pID = self.pID[cond]
        self.ptypes = self.ptypes[cond]
        if self.f is None:
            self.f = None
        else:
            self.f = self.f[cond]
        cond = np.where((mask4part[self.simplices[:,0]] == 1.) & 
                        (mask4part[self.simplices[:,1]] == 1.) & 
                        (mask4part[self.simplices[:,2]] == 1.) & 
                        (mask4part[self.simplices[:,3]] == 1.))[0]
        self.simplices = keep_idx[self.simplices[cond]]
        self.simptypes = self.simptypes[cond]
    
        
    def run(self, boundary: List[float], apply_filter: bool = False) -> None:
        """
        Constructes the Delaunay tesselation from the border input points.

        Parameters
        ----------
        Boundary : list
            A list containing the minimum and maximum x values, given in this order: [xmin, xmax, ymin, ymax].
        apply_filter : bool, optional
            Apply filter on boundary particles, this will cut out all points and simplices that are beyond or cross
            the boundary. Default set to False.
        """
        self._concatenate()
        if apply_filter:
            self._filter(boundary)
        self.dtfe = Delaunay3D()
        self.dtfe.setup(self.x, self.y, self.z, boundary, mass=self.mass)
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
            (self.pID[self.dtfe.simplices[:,2]] == self.pID[self.dtfe.simplices[:,3]]) &
            (self.ptypes[self.dtfe.simplices[:,0]] == 0.) & 
            (self.ptypes[self.dtfe.simplices[:,1]] == 0.) & 
            (self.ptypes[self.dtfe.simplices[:,2]] == 0.) & 
            (self.ptypes[self.dtfe.simplices[:,3]] == 0.)
        )[0]
        self.dtfe.simptypes[cond] = -1.
    
    
    def estimate(self, x: Union[float, np.ndarray], y: Union[float, np.ndarray], z: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """
        Estimates the field from input coordinates, if coordinate is not in simplex
        or the simplex is of type != 1 then a Nan is returned.

        Parameters
        ----------
        x : float or array
            X-axis coordinate.
        y : float or array
            Y-axis coordinate.
        z : float or array
            Z-axis coordinate
        """
        return self.dtfe.estimate(x, y, z)
    
    
    def get_border(self) -> dict:
        """
        Returns border points and simplices from a merged delaunay tesselation.
        """
        return self.dtfe.get_border()
    
    
    def clean(self) -> None:
        """
        Reinitialises the class.
        """
        return self.__init__()