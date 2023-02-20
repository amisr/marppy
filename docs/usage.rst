.. _Usage:

Usage
=====

The Marp class contains all functions for converting between geodetic, Apex, and MARP coordinates.  Note that Marp inherits the `Apex <https://apexpy.readthedocs.io/en/latest/autoapi/generated/apexpy/index.html#apexpy.Apex>`_ class, so all Apex functions are also available through Marp.

Initialize the Marp Class
-------------------------

The Marp class requires a date, reference height, and either the location of the MARP north pole or null island ((0,0) point) when initialized.  The date can be provided as a decimal year, date object, or datetime object.  The reference height defaults to 0 km when not specified.  A list of three floating point values should be provided for either `pole` or `null`.  For `pole`, this should be the latitude and longitude of the new MARP north pole, plus and additional angle that the coordinate system should be rotated around that new pole so that the prime meridian falls in the correct location.  For `null`, this should be the latitude and longitude of the new MARP null island, plus the bearing direction of North from that point. The `coords` keyword defines whether the `pole` or `null` positions are specified in Apex or geodetic coordinates (default Apex).  If the position is specified in geodetic coordinates, `alt` should be used to specify the point's altitude above the Earth's surface (defaults to the reference height). Some examples are provided below.

Specify location of north pole in Apex coordinates.

.. code-block:: python

  import marppy
  import datetime as dt
  M = marppy.Marp(dt.datetime(2021,7,31), poll=[54.7,33.8,10.0])

Specify location of null island in geodetic coordinates.

.. code-block:: python

  import marppy
  M = marppy.Marp(2019.4, null=[74.7,-94.9,26.0], alt=300., coords='geo')

Note: For most use cases, specifying null will be most intuitive.

Convert between Apex and Marp Coordinates
-----------------------------------------

Use the `apex2marp()` and `marp2apex()` functions to convert between Apex and Marp coordinates.

.. code-block:: python

  mlat, mlon = M.apex2marp(65.0, 17.0)
  print(mlat, mlon)
  31.567175854041146, 12.261423536020837

.. code-block:: python

  alat, alon = M.marp2apex(31.6, 12.3)
  print(alat, alon)
  64.95365215608916, 17.00718613610059

Convert between Geodetic and Marp Coordinates
---------------------------------------------

Use the `geo2marp()` and `marp2geo()` functions to convert between geodetic and Marp coordinates.  Note that height must be specified when converting to and from geodetic coordinates.

.. code-block:: python

  mlat, mlon = M.geo2marp(84.3, -136.8, height=450.)
  print(mlat, mlon)
  8.116045001238708, -7.773699770282615

.. code-block:: python

  glat, glon, _ = M.marp2geo(8.1, -7.8, height=450.)
  print(glat, glon)
  84.31710052490234, -137.10552978515625
