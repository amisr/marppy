Overview
==========

This is a python implementation of conversions to and from the Modified Apex - Rotated Pole (MARP) coordinate system.  MARP is based on the Apex coordinate system (Richmond, 1995), but performs a  transformation to rotate the pole of the coordinate system to an arbitrary location.  This is advantageous for doing calculations at high latitudes because it removes some of the geometric complications.

This code is heavily based on `apexpy <https://github.com/aburrell/apexpy>`_.


Installation
------------

The easiest way to install marppy is from PyPI::

  pip install marppy


Usage
-----

All functionality is available through the marppy.Marp class. Instantiate the class with the date and the new MARP "origin" (mlat=0, mlon=0) coordinates that dictate how the coordinate system will be rotated, then use the various conversion routines::

  >>> from marppy import Marp
  >>> M = Marp(date=2022.5, lam0=80., phi0=30.)
  >>> # geo to marp
  >>> M.geo2marp(50.6, 27.6, 300.)
  (-4.537915474007271, 39.96811884264296)
  >>> # marp to apex
  >>> M.marp2apex(5.8, 10.2)
  (79.06102883542528, 98.18833972513148)
  >>> # apex to marp - entering the original rotation coordinates should return (0,0)
  >>> M.apex2marp(80., 30.)
  (0.0, 0.0)


Documentation
-------------

Full documentation is available on `readthedocs <https://marppy.readthedocs.io/en/latest/index.html>`_.
