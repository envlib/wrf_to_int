from enum import Enum

import numpy as np

# import fortran_io as f_io
from . import fortran_io as f_io

class Projections( Enum ) :
  LATLON       = 0      # Cylindrical equidistant
  MERC         = 1      # Mercator
  LC           = 3      # Lambert conformal
  GAUSS        = 4      # Gaussian
  PS           = 5      # Polar stereographic
  CASSINI      = 6

class IntermediateFile( object ):
  def __init__( self, prefix, datestr ) :
    self.prefix_ = prefix
    self.datestr_ = datestr
    self.filename_ = self.prefix_.strip() + ":" + self.datestr_.strip()
    self.file_ = open( self.filename_, "wb" )

  def close( self ) :
    self.file_.close()

  def write_next_met_field(
                            self,
                            version,
                            nx, ny,
                            iproj,
                            xfcst,
                            xlvl,
                            startlat, startlon, starti, startj,
                            deltalat, deltalon, dx, dy,
                            xlonc, truelat1, truelat2,
                            earth_radius,
                            is_wind_grid_rel,
                            field,
                            hdate,
                            units,
                            map_source,
                            desc,
                            slab
                            ) :

    f_io.unfmt_ftn_rec_write(
                              [ version ],
                              fmt=f_io.StructFormats.INT32,
                              file=self.file_
                            )

    if field == "GHT" :
      field = "HGT"

    if version == 5 :
      # FP check okay here?
      startloc = None
      if starti == 1.0 and startj == 1.0 :
        startloc = "SWCORNER"
      else :
        startloc = "CENTER  "

      f_io.unfmt_ftn_rec_write(
                                [
                                  hdate.ljust(24)[0:24],
                                  xfcst,
                                  map_source.ljust(32)[0:32],
                                  field.ljust(9)[0:9],
                                  units.ljust(25)[0:25],
                                  desc.ljust(46)[0:46],
                                  xlvl, nx, ny, iproj.value ],
                                fmt=( [ f_io.StructFormats.ARRCHAR ] +
                                      [ f_io.StructFormats.FP32 ] +
                                      [ f_io.StructFormats.ARRCHAR ] * 4 +
                                      [ f_io.StructFormats.FP32 ] +
                                      [ f_io.StructFormats.INT32 ] * 3 ),
                                file=self.file_
                                )
      if iproj == Projections.LATLON or iproj == Projections.GAUSS:
        f_io.unfmt_ftn_rec_write(
                                  [
                                    startloc.ljust(8)[0:8],
                                    startlat,
                                    startlon,
                                    deltalat,
                                    deltalon,
                                    earth_radius ],
                                  fmt=[ f_io.StructFormats.ARRCHAR ] + [ f_io.StructFormats.FP32 ] * 5,
                                  file=self.file_
                                )
      elif iproj == Projections.MERC :
        f_io.unfmt_ftn_rec_write(
                                  [
                                    startloc.ljust(8)[0:8],
                                    startlat,
                                    startlon,
                                    dx,
                                    dy,
                                    truelat1,
                                    earth_radius ],
                                  fmt=[ f_io.StructFormats.ARRCHAR ] + [ f_io.StructFormats.FP32 ] * 6,
                                  file=self.file_
                                )
      elif iproj == Projections.LC :
        f_io.unfmt_ftn_rec_write(
                                  [
                                    startloc.ljust(8)[0:8],
                                    startlat,
                                    startlon,
                                    dx,
                                    dy,
                                    xlonc,
                                    truelat1,
                                    truelat2,
                                    earth_radius ],
                                  fmt=[ f_io.StructFormats.ARRCHAR ] + [ f_io.StructFormats.FP32 ] * 8,
                                  file=self.file_
                                )
      elif iproj == Projections.PS :
        f_io.unfmt_ftn_rec_write(
                                  [
                                    startloc.ljust(8)[0:8],
                                    startlat,
                                    startlon,
                                    dx,
                                    dy,
                                    xlonc,
                                    truelat1,
                                    earth_radius ],
                                  fmt=[ f_io.StructFormats.ARRCHAR ] + [ f_io.StructFormats.FP32 ] * 7,
                                  file=self.file_
                                )

      f_io.unfmt_ftn_rec_write( [ is_wind_grid_rel ], fmt=f_io.StructFormats.INT32, file=self.file_ )
      f_io.unfmt_ftn_rec_write( slab,                 fmt=f_io.StructFormats.FP32, file=self.file_ )
      return 0
    else :
      print( "Didn't recognize format number " + str( version ) )
      return 1


class MapProjection:
    """Stores parameters of map projections as used in the WPS intermediate file format."""

    def __init__(self, projType, startLat, startLon, startI, startJ, deltaLat, deltaLon,
                 dx=0.0, dy=0.0, truelat1=0.0, truelat2=0.0, xlonc=0.0):
        self.projType = projType
        self.startLat = startLat
        self.startLon = startLon
        self.startI = startI
        self.startJ = startJ
        self.deltaLat = deltaLat
        self.deltaLon = deltaLon
        self.dx = dx
        self.dy = dy
        self.truelat1 = truelat1
        self.truelat2 = truelat2
        self.xlonc = xlonc


def _ensure_str(val):
    """Decode bytes to str if needed (h5netcdf attrs may return bytes)."""
    if isinstance(val, bytes):
        return val.decode()
    return str(val)


def write_slab(intfile, slab, xlvl, proj, WPSname, hdate, units, map_source, desc):
    """Write a 2D field slab to an opened WPS intermediate file.

    Handles both regular numpy arrays and masked arrays. NaN values are
    converted to the WPS missing value sentinel (-1.0e30).
    String parameters are decoded from bytes if needed (for h5netcdf compatibility).

    WPS intermediate file convention: earth_radius in km, dx/dy in km.
    """
    missing_value = -1.0e30
    data = np.squeeze(np.asarray(slab, dtype=np.float64))
    if not isinstance(data, np.ma.MaskedArray):
        data = np.ma.array(data, mask=np.isnan(data))
    _ = intfile.write_next_met_field(
        5, data.shape[1], data.shape[0], proj.projType, 0.0, xlvl,
        proj.startLat, proj.startLon, proj.startI, proj.startJ,
        proj.deltaLat, proj.deltaLon, proj.dx, proj.dy, proj.xlonc,
        proj.truelat1, proj.truelat2, 6371.229, 0, _ensure_str(WPSname),
        _ensure_str(hdate), _ensure_str(units), _ensure_str(map_source),
        _ensure_str(desc), data.filled(missing_value))
