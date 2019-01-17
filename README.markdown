Package Name
============

__Version__: 2018-09-03

__Author(s)__
Kevin Chu `<kevin@velexi.com>`

-------------------------------------------------------------------------------

Preface
-------

* This document is written using GitHub Markdown.  It is best viewed using a
  Markdown previewer that supports GitHub Markdown:

  - [Atom Editor](https://atom.io/)

  - [Chrome Markdown Preview Plus][chrome-markdown-preview-plus] extension
    using Github theme

  - [Firefox Markdown Viewer][firefox-markdown-viewer] add-on

-------------------------------------------------------------------------------

Table of Contents
-----------------

1. [Overview][#1]

2. [Installation][#2]

3. [License][#3]

4. [Known Limitations and Issues][#4]

------------------------------------------------------------------------------

1. Overview
-----------

TODO

------------------------------------------------------------------------------

2. Installation
---------------

For installation instructions, see `INSTALL.markdown`.

# 2.1. Software Dependencies

* SAMRAI (>=3.12)
  - Open MPI (>=3.1.1)
  - HDF5 (>=1.10.2)

# 2.2. Installation Instructions

* Install SAMRAI

  - Download SAMRAI

    * https://computation.llnl.gov/projects/samrai/software#download

  - Install Open MPI

    * TODO: verify that recent versions of Open MPI are compatible

  - Install HDF5

    * TODO: verify that recent versions of HDF5 are compatible

------------------------------------------------------------------------------

3. License
----------

TODO

------------------------------------------------------------------------------

4. Algorithm Description
------------------------

### Reinitialization

* Reinitialization can be disabled by setting `reinitialization_interval`
  equal to 0 in the configuration file.

### Boundary Conditions

* At non-periodic boundaries,

  - homogeneous Neumann boundary conditions are applied ?

  - flow boundary conditions are applied ?

### Regridding

* Regridding can be disabled by setting `regrid_interval` equal to 0 in the
  configuration file.

------------------------------------------------------------------------------

5. Known Limitations and Issues
-------------------------------

* Visualization with VisIt.

  - When generating contour plots, VisIt assumes that data is node-centered.
    PQS-AMR generates data that is cell-centered. To avoid generating
    contour artifacts, use the following steps.

    * Create node-centered variable ... Expressions ... recenter(phi)

    * Create contour plot using node-centered variable (instead of raw
      cell-centered data produced by PQS-AMR).

  - _References_

    * http://visitusers.org/index.php?title=VisIt-tutorial-data-analysis

    * https://elist.ornl.gov/mailman/htdig/visit-users/2010-May/006978.html

* Only one Solver may be instantiated per program because PQS package does
  not support customizable object names (and multiple objects with the same
  name may not be registered with the SAMRAI RestartManager).

* TODO

-------------------------------------------------------------------------------

[-----------------------------INTERNAL LINKS-----------------------------]: #

[#1]: #1-overview

[#2]: #2-installation

[#3]: #3-license

[#4]: #4-known-limitations-and-issues

[-----------------------------EXTERNAL LINKS-----------------------------]: #

[chrome-markdown-preview-plus]: https://chrome.google.com/webstore/detail/markdown-preview-plus/febilkbfcbhebfnokafefeacimjdckgl
[firefox-markdown-viewer]: https://addons.mozilla.org/en-US/firefox/addon/markdown-viewer/

[include-what-you-use]: http://include-what-you-use.org/
