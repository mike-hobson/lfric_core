#!/usr/bin/env python
# *****************************COPYRIGHT*******************************
# (C) Crown copyright Met Office. All rights reserved.
# For further details please refer to the file LICENCE
# which you should have received as part of this distribution.
# *****************************COPYRIGHT*******************************
"""Generate a UM grid file in SCRIP format from a namelist file."""
import os
from netCDF4 import Dataset
from um_utils.cutout import CoordRotator
import numpy
import iris
import iris.fileformats
import iris.analysis
import iris.coord_systems
import f90nml


class GRID:
    "Holds information about the UM grid"
    def __init__(self):
        self.vname = None     # Title to use in netCDF file
        self.grid = None      # UM grid type: P/U/V points at grid cell centres
        self.l_area = None    # If not 'no', output grid cell surface areas

        self.dlon = None      # Longitude of grid cell centres
        self.dlat = None      # Lattitude of grid cell centres
        self.dclo = None      # Longitude of grid cell corners
        self.dcla = None      # Lattitude of cell corners
        self.darea = None     # Area of each grid cell

        self.ucoor = None     # Coordinate units
        self.uarea = None     # Area units of each grid cell

        self.shape = None     # Shape of grid

        # Information read from namelist
        self.rotated_grid = None  # Flag to indicated grid is rotated
        self.pole_lon = None  # Lambda pole
        self.pole_lat = None  # Phi pole
        self.delx = None  # Delta lambda
        self.dely = None  # Delta phi
        self.npx = None  # Points lambda
        self.npy = None  # Points phi
        self.xorigin = None  # lambda origin
        self.yorigin = None  # Phi origin

    def read_namelist(self, grid_namelist):
        grid_def = f90nml.read(grid_namelist)

        self.rotated_grid = grid_def['grid']['rotated']
        self.pole_lon = grid_def['grid']['lambda_pole']
        self.pole_lat = grid_def['grid']['phi_pole']
        self.delx = grid_def['grid']['delta_lambda_targ']
        self.dely = grid_def['grid']['delta_phi_targ']
        self.npx = grid_def['grid']['points_lambda_targ']
        self.npy = grid_def['grid']['points_phi_targ']
        self.xorigin = grid_def['grid']['lambda_origin_targ']
        self.yorigin = grid_def['grid']['phi_origin_targ']


def UM_guess_bounds(cube):
    """
    Estimate area of grid cells by guessing bounds of 'latitude'
    and 'longitude' coordinates and return areas as an array.

    Args:

    * cube:
        iris.cube.Cube

    Returns:
        iris.cube.Cube with bounds on coordinates

    """

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()

    return cube


def UM_sort_cube_bounds(cube):
    '''
    Sort out the coordinate bounds on a cube, add bounds for each coordinate
    and make sure that latitude bounds to not go outside [-90,90]

    Arguments:

      * cube:
         iris cube to deal with
    '''

    for coord in cube.dim_coords:
        if not coord.has_bounds():
            coord.guess_bounds()

    if len(cube.coord_dims('latitude')) > 0:
        latbounds = cube.coord('latitude').bounds.copy()
        latbounds = numpy.clip(latbounds, -90., 90.)
        cube.coord('latitude').bounds = latbounds

    if len(cube.coord_dims('longitude')) > 0:
        if cube.coord('longitude').bounds.max() == \
           cube.coord('longitude').bounds.min()+360:
            cube.coord('longitude').circular = True


def UM_get_grid(cube):
    """Calculate UM grid"""
    lon = numpy.resize(cube.coord('longitude').points, (cube.shape))
    lat = numpy.zeros((cube.shape))
    for i in range(cube.shape[1]):
        lat[:, i] = cube.coord('latitude').points
    return lon, lat


def UM_gen_area(cube):
    """Calculate area for each grid"""
    return iris.analysis.cartography.area_weights(cube)


def UM_create_cube(my_grid):
    """Create UM cube with the grid definitions"""
    tiny = 1.0e-10

    xcoord = numpy.arange(my_grid.xorigin, my_grid.xorigin +
                          (my_grid.npx * my_grid.delx - tiny), my_grid.delx)
    ycoord = numpy.arange(my_grid.yorigin, my_grid.yorigin +
                          (my_grid.npy * my_grid.dely - tiny), my_grid.dely)

    if my_grid.grid == 'U':
        xcoord -= .5*my_grid.delx
        print(my_grid.npx, xcoord.shape, 'U-grid dims')

    if my_grid.grid == 'V':
        ycoord = numpy.arange(my_grid.yorigin, my_grid.yorigin +
                              ((my_grid.npy + 1) * my_grid.dely - tiny),
                              my_grid.dely)
        print(my_grid.npy, ycoord.shape, 'V-grid dims')
        ycoord -= .5*my_grid.dely
        ycoord[0] = numpy.float32(ycoord[0])
        ycoord[-1] = numpy.float32(ycoord[-1])

    if min(xcoord) < 0.:
        raise Exception('lowest longitude less than zero')

    if min(ycoord) < -90. or min(ycoord) > 90.:
        print(min(ycoord), min(ycoord))
        print('lowest latitude less than -90 or greater than 90')

    mlong = iris.coords.DimCoord(xcoord, standard_name='longitude',
                                 circular=False, units='degrees',
                                 coord_system=iris.coord_systems.GeogCS(
                                     iris.fileformats.pp.EARTH_RADIUS))

    mlat = iris.coords.DimCoord(ycoord, standard_name='latitude',
                                units='degrees',
                                coord_system=iris.coord_systems.GeogCS(
                                    iris.fileformats.pp.EARTH_RADIUS))

    data = numpy.zeros((len(ycoord), len(xcoord)))
    cube = iris.cube.Cube(data, dim_coords_and_dims=[(mlat, 0), (mlong, 1)],
                          var_name='cube')
    UM_guess_bounds(cube)

    UM_sort_cube_bounds(cube)

    return cube


def UM_lat_corners(cube):
    """Calculate UM latitude corners"""

    lat_bounds = numpy.zeros((cube.shape))
    corners = numpy.zeros((cube.shape[0], cube.shape[1], 4))
    for i in range(cube.shape[1]):
        lat_bounds[:, i] = cube.coord('latitude').bounds[:, 0]
    corners[:, :, 0] = lat_bounds
    corners[:, :, 1] = lat_bounds
    for i in range(cube.shape[1]):
        lat_bounds[:, i] = cube.coord('latitude').bounds[:, 1]
    corners[:, :, 2] = lat_bounds
    corners[:, :, 3] = lat_bounds

    return corners


def UM_lon_corners(cube):
    """Calculate UM longitude corners"""

    lon_bounds = numpy.zeros((cube.shape))
    corners = numpy.zeros((cube.shape[0], cube.shape[1], 4))
    for i in range(cube.shape[0]):
        lon_bounds[i, :] = cube.coord('longitude').bounds[:, 0]
    corners[:, :, 0] = lon_bounds
    corners[:, :, 3] = lon_bounds
    for i in range(cube.shape[0]):
        lon_bounds[i, :] = cube.coord('longitude').bounds[:, 1]
    corners[:, :, 1] = lon_bounds
    corners[:, :, 2] = lon_bounds

    return corners


def corners_transformation(field):
    """Transform corners to correct order"""
    return field.transpose((1, 2, 0))


def transform(fin, name):
    """Remove one dimension from the array"""
    shape = fin.shape
    if len(shape) == 2:
        if(shape[0] == 4):
            fout = numpy.transpose(fin)
        elif(shape[1] == 4):
            fout = fin
        else:
            fout = numpy.zeros((shape[0]*shape[1]), dtype=numpy.float64)
            fout = fin.flatten()
    elif len(shape) == 3:
        if shape[0] == 4:
            fin = corners_transformation(fin)
            shape_o = shape
            shape = fin.shape
            print('transforming corner data for ', name, ' from ', shape_o,
                  ' to', shape)
        fout = numpy.zeros((shape[0]*shape[1], 4), dtype=numpy.float64)
        fout = numpy.reshape(fin, (shape[0]*shape[1], 4))
    elif len(shape) == 1:
        fout = fin
    else:
        raise Exception('problem in transform')
    return fout


def transform_and_write(filename, my_grid):
    """Transform arrays and write netcdf output"""

    rank = len(my_grid.shape)
    mrank = my_grid.shape

    if rank == 2:
        grid_size = my_grid.shape[0]*my_grid.shape[1]
    else:
        raise Exception('transform_and_write: ' +
                        'wrong number of dimensions in input data')

    corners = 4
    fileout = filename
    outfile = Dataset(fileout, 'w')
    outfile.createDimension('grid_rank', numpy.int32(rank))
    outfile.createDimension('grid_size', numpy.int32(grid_size))
    outfile.createDimension('grid_corners', numpy.int32(corners))
    rankout = outfile.createVariable('grid_dims',
                                     numpy.dtype('int32').char,
                                     ('grid_rank', ))
    rankout.long_name = "grid_dims"
    rankout[:] = numpy.asarray(mrank).astype(numpy.int32)
    latout = outfile.createVariable('grid_center_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    latout.long_name = "grid_center_lat"
    latout.units = str(my_grid.ucoor)
    latout[:] = transform(my_grid.dlat, 'lat')
    lonout = outfile.createVariable('grid_center_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', ))
    lonout.long_name = "grid_center_lon"
    lonout.units = str(my_grid.ucoor)
    lonout[:] = transform(my_grid.dlon, 'lon')
    maskout = outfile.createVariable('grid_imask',
                                     numpy.dtype('int32').char,
                                     ('grid_size', ))
    maskout.long_name = "grid_imask"
    maskout.units = "unitless"
    maskout[:] = numpy.int32(1)
    if my_grid.l_area != 'no':
        srfout = outfile.createVariable('grid_area',
                                        numpy.dtype('float64').char,
                                        ('grid_size', ))
        srfout.long_name = "grid surface m^2"
        srfout.units = str(my_grid.uarea)
        srfout[:] = transform(my_grid.darea, 'area')
    claout = outfile.createVariable('grid_corner_lat',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    claout.long_name = "grid_corner_lat"
    claout.units = str(my_grid.ucoor)
    claout[:, :] = transform(my_grid.dcla, 'cla')
    cloout = outfile.createVariable('grid_corner_lon',
                                    numpy.dtype('float64').char,
                                    ('grid_size', 'grid_corners'))
    cloout.long_name = "grid_corner_lon"
    cloout.units = str(my_grid.ucoor)
    cloout[:, :] = transform(my_grid.dclo, 'clo')
    outfile.title = my_grid.vname
    outfile.close()
    return


def get_data(my_grid):

    my_grid.vname = "UM_grid"
    my_grid.grid = os.environ.get('GRID')
    my_grid.l_area = os.environ.get('LAREA')
    my_grid.read_namelist(os.environ.get('NLIST_FILE'))

    cube = UM_create_cube(my_grid)
    
    rotated_coords = CoordRotator(my_grid.pole_lon, my_grid.pole_lat)

    # Get longitude and latitude of cell centres
    x, y = UM_get_grid(cube)
    if my_grid.rotated_grid:
        my_grid.dlon, my_grid.dlat = unrotate_coords(x, y, rotated_coords)
    else:
        my_grid.dlon, my_grid.dlat = x, y

    # Get longitude and latitude of cell corners
    x = UM_lon_corners(cube)
    y = UM_lat_corners(cube)
    if my_grid.rotated_grid:
        my_grid.dclo, my_grid.dcla = unrotate_coords(x, y, rotated_coords)
    else:
        my_grid.dclo, my_grid.dcla = x, y

    # Get area
    my_grid.darea = UM_gen_area(cube)

    # Set units
    my_grid.ucoor = "degrees"
    my_grid.uarea = "m2"

    # Set shape of grid
    my_grid.shape = numpy.asarray(my_grid.dlat.shape[::-1])

    return my_grid


def unrotate_coords(lon_rot, lat_rot, cr):
    """Unrotate coordinates"""

    if lon_rot.shape != lat_rot.shape:
        raise Exception('unrotate_coords: lon/lat arrays have different sizes')
    else:
        lon = numpy.empty(lon_rot.shape)
        lat = numpy.empty(lat_rot.shape)

    for index, _ in numpy.ndenumerate(lon_rot):
        x = lon_rot[index]
        y = lat_rot[index]
        lon[index], lat[index] = cr.unrotate(x, y)
        if lon[index] < 0.0:
            lon[index] = lon[index] + 360.0

    return lon, lat


if __name__ == "__main__":

    um_grid = GRID()
    um_grid = get_data(um_grid)
    transform_and_write(os.environ["GRID_PATH_UM"], um_grid)
