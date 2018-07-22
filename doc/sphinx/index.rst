**************************************************
Pairinteraction - A Rydberg Interaction Calculator
**************************************************

.. image:: https://travis-ci.org/pairinteraction/pairinteraction.svg?branch=master
   :target: https://travis-ci.org/pairinteraction/pairinteraction
   :alt: Travis Build Status
.. image:: https://ci.appveyor.com/api/projects/status/t5l4unwt210gq3al/branch/master?svg=true
   :target: https://ci.appveyor.com/project/pairinteraction/pairinteraction/branch/master
   :alt: AppVeyor Build Status
.. image:: images/arXiv-badge.svg
   :target: https://arxiv.org/abs/1612.08053
   :alt: arXiv:1612.08053
.. image:: images/license-badge.svg
  :target: https://opensource.org/licenses/Apache-2.0
  :alt: License

The *pairinteraction* software calculates properties of Rydberg systems.
The software consists of a C++/Python library and a graphical user interface for pair potential calculations.
For usage examples visit the `tutorials`_ section of the documentation.
Stay tuned by signing up for the newsletter so whenever there are updates to the software or new publications about pairinteraction we can contact you. To subscribe click `here`_.

.. _tutorials: https://pairinteraction.github.io/pairinteraction/sphinx/html/tutorials.html

.. _here: https://goo.gl/forms/4bmz3qeuLjKfRlWJ3

Please cite us
    Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth,
    *Tutorial: Calculation of Rydberg interaction potentials*,
    `J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017) <https://doi.org/10.1088/1361-6455/aa743a>`_, `arXiv:1612.08053 <https://arxiv.org/abs/1612.08053>`_

See works citing pairinteraction on `Google Scholar`_ and the `ADS Digital Library`_.

.. _Google Scholar: https://scholar.google.com/scholar?cites=5795867423675717201
.. _ADS Digital Library: http://adsabs.harvard.edu/cgi-bin/nph-ref_query?bibcode=2017JPhB...50m3001W&amp;refs=CITATIONS&amp;db_key=PHY

Installation
============

Binary builds are available for GNU/Linux, Mac OS X, and Windows through `GitHub Releases`_. For more information, read the `installation docs`_.

.. _GitHub Releases: https://github.com/pairinteraction/pairinteraction/releases
.. _installation docs: https://pairinteraction.github.io/pairinteraction/sphinx/html/installation.html

Main Features
=============

The pairinteraction software simulates systems of one or two Rydberg atoms. It supports a lot of neat features, the major of which are listed below:

* Stark and Zeeman maps for single atom states
* Pair potentials taking into account electromagnetic fields in arbitrary directions
* Multipole interaction up to arbitrary order
* Numerical and analytical methods for radial matrix elements
* Automatic exploitation of symmetries
* Calculation of :math:`C_6` coefficients in degenerate subspaces
* High performance C++ backend, Python interface with NumPy support
* MATLAB compatible data export

Screenshots
===========

.. raw:: html
    
    <style>
    .slideshow-container {
        display: block;
        position: relative
    }
    
    .slide {
        max-width: 100%;
        max-height: 100%;
    }
    
    .slide-text {
        position: absolute;
        font-size: .8em;
        bottom: 10px;
        right: 5px;
        margin-left: 5px;
        padding: .25em;
        color: #fafafa;
        background: rgba(102,102,102,0.6);
        border-radius: 2px;
    }
    
    .slide-text a {
        color: #fafafa;
        text-decoration: underline;
    }
    
    .slideshow-btn {
        display: inline-block;
        position: absolute;
        top: 50%;
        font-weight: bold;
        text-align: center;
        color: #fafafa;
        margin-top: -.75em;
        font-size: 2em;
        padding: 10px 10px 15px 10px;
        background: rgba(102,102,102,0.6);
        opacity: .3;
    }
    
    .slideshow-btn:hover {
        cursor: pointer;
        opacity: 1;
    }
    </style>
    
    <script>
    var slideshow = {
        idx: 0,
    
        start: function () {
            this.show(0);
        },
    
        step: function (n) {
            this.show(this.idx += n);
        },
    
        show: function(n) {
            var slides = document.getElementsByClassName("slide");
            this.idx = Math.abs(n % slides.length);
            for (var i = 0; i < slides.length; i++) {
                slides[i].style.display = "none"; 
            }
            slides[this.idx].style.display = "block";
        }
    }
    </script>
    
    <p>
    <div class="slideshow-container">
      <div class="slide">
        <a href="slides/screen-win64.png"><img src="slides/screen-win64.png"></a>
        <div class="slide-text">Main window with sample configuration on Windows 10 64-bit.</div>
      </div>
      <div class="slide">
        <a href="slides/screen-osx-pairpotential.jpg"><img src="slides/screen-osx-pairpotential.jpg"></a>
        <div class="slide-text">Main window with sample configuration on Mac OS X El Capitan.</div>
      </div>
      <div class="slide">
        <a href="slides/screen-osx-starkmap.jpg"><img src="slides/screen-osx-starkmap.jpg"></a>
        <div class="slide-text"><a href="https://dx.doi.org/10.1002/andp.19143480702">Stark map</a> for a single Rydberg atom.</div>
      </div>
      <div class="slide">
        <a href="slides/screen-osx-zeemanmap.jpg"><img src="slides/screen-osx-zeemanmap.jpg"></a>
        <div class="slide-text"><a href="https://dx.doi.org/10.1103/PhysRev.55.52">Quadratic Zeeman effect</a> for a single Rydberg atom.</div>
      </div>
      <div class="slideshow-btn" style="left:0"  onclick="slideshow.step(-1)">«</div>
      <div class="slideshow-btn" style="right:0" onclick="slideshow.step(+1)">»</div>
    </div>
    <script>slideshow.start();</script>
    </p>

.. toctree::
   :numbered:
   :maxdepth: 2
   :caption: Documentation

   installation.rst
   tutorials.rst

.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: Indices and Tables

    modules.rst

Credits
=======

The pairinteraction software was originally developed at the `5th Institute of Physics`_ and the `Institute for Theoretical Physics III`_ of the University of Stuttgart, Germany.
Currently it is maintained by developers at the `Institute for Theoretical Physics III`_ of the University of Stuttgart in Germany, the `Department of Physics`_ of the University of Otago in New Zealand,
and the `Department of Physics, Chemistry and Pharmacy`_ of the University of Southern Denmark in Denmark.

.. _5th Institute of Physics: http://www.pi5.uni-stuttgart.de/
.. _Institute for Theoretical Physics III: http://www.itp3.uni-stuttgart.de/
.. _Department of Physics: http://www.otago.ac.nz/physics/index.html
.. _Department of Physics, Chemistry and Pharmacy: http://www.sdu.dk/en/fkf
