**************************************************
Pairinteraction - A Rydberg Interaction Calculator
**************************************************

|travis| |appveyor| |codecov| |pypi| |arxiv| |license|

The *pairinteraction* software calculates properties of Rydberg systems.
The software consists of a C++/Python library and a graphical user interface for pair potential calculations.
For usage examples visit the :ref:`tutorials <Tutorials>` section of the documentation.
Stay tuned by `signing up`_ for the newsletter so whenever there are updates to the software or new publications about pairinteraction we can contact you.
If you have a question that is related to problems, bugs, or suggests an improvement, consider raising an :github:`issue <issues>` on :github:`GitHub <>`.

The pairinteraction library is licensed under the `LGPL v3`_. The extension for calculating
radial wave functions using Whittaker functions and the graphical user interface are licensed under the `GPL v3`_.
The GPL v3 also applies to the combined work and all provided binary builds.

.. _signing up: https://goo.gl/forms/4bmz3qeuLjKfRlWJ3
.. _LGPL v3: https://www.gnu.org/licenses/lgpl-3.0.html
.. _GPL v3: https://www.gnu.org/licenses/gpl-3.0.html

Please cite us
    Sebastian Weber, Christoph Tresp, Henri Menke, Alban Urvoy, Ofer Firstenberg, Hans Peter Büchler, Sebastian Hofferberth,
    *Tutorial: Calculation of Rydberg interaction potentials*,
    `J. Phys. B: At. Mol. Opt. Phys. 50, 133001 (2017) <https://doi.org/10.1088/1361-6455/aa743a>`_, `arXiv:1612.08053 <https://arxiv.org/abs/1612.08053>`_

    .. code-block:: bibtex

        @article{Weber2017,
          author = {Weber, Sebastian and Tresp, Christoph and Menke, Henri and Urvoy, Alban and Firstenberg, Ofer and B{\"u}chler, Hans Peter and Hofferberth, Sebastian},
          title = {{Tutorial: Calculation of Rydberg interaction potentials}},
          journal = {J. Phys. B: At. Mol. Opt. Phys.},
          volume = {50},
          number = {13},
          pages = {133001},
          year = {2017},
          doi = {10.1088/1361-6455/aa743a},
          url = {https://doi.org/10.1088/1361-6455/aa743a}
        }

See works citing pairinteraction at the `ADS Digital Library`_ and on `Google Scholar`_.

.. _Google Scholar: https://scholar.google.com/scholar?cites=5795867423675717201
.. _ADS Digital Library: https://ui.adsabs.harvard.edu/#abs/2017JPhB...50m3001W/citations

Installation
============

Binary builds are available for GNU/Linux, Mac OS X, and Windows through :github:`GitHub Releases <releases>`. For using pairinteraction as a Python 3 library,
we recommend the installation via pip by calling ``pip install pairinteraction``. If pairinteraction was installed from the command line, the graphical user
interface can be started by executing ``start_pairinteraction_gui``. For more information, read the :ref:`installation docs <Installation>`.

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
the `Institute of Physics`_ of the University of Rostock in Germany, and the `Department of Physics, Chemistry and Pharmacy`_ of the University of Southern Denmark in Denmark.

.. _5th Institute of Physics: http://www.pi5.uni-stuttgart.de/
.. _Institute for Theoretical Physics III: http://www.itp3.uni-stuttgart.de/
.. _Department of Physics: http://www.otago.ac.nz/physics/index.html
.. _Department of Physics, Chemistry and Pharmacy: http://www.sdu.dk/en/fkf
.. _Institute of Physics: https://www.physik.uni-rostock.de/
