Overview
========

TODO UPDATE THIS

Features
--------

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
-----------

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
        <a href="_static/slides/screen-win64.png"><img src="_static/slides/screen-win64.png"></a>
        <div class="slide-text">Main window with sample configuration on Windows 10 64-bit.</div>
      </div>
      <div class="slide">
        <a href="_static/slides/screen-osx-pairpotential.jpg"><img src="_static/slides/screen-osx-pairpotential.jpg"></a>
        <div class="slide-text">Main window with sample configuration on Mac OS X El Capitan.</div>
      </div>
      <div class="slide">
        <a href="_static/slides/screen-osx-starkmap.jpg"><img src="_static/slides/screen-osx-starkmap.jpg"></a>
        <div class="slide-text"><a href="https://dx.doi.org/10.1002/andp.19143480702">Stark map</a> for a single Rydberg atom.</div>
      </div>
      <div class="slide">
        <a href="_static/slides/screen-osx-zeemanmap.jpg"><img src="_static/slides/screen-osx-zeemanmap.jpg"></a>
        <div class="slide-text"><a href="https://dx.doi.org/10.1103/PhysRev.55.52">Quadratic Zeeman effect</a> for a single Rydberg atom.</div>
      </div>
      <div class="slideshow-btn" style="left:0"  onclick="slideshow.step(-1)">«</div>
      <div class="slideshow-btn" style="right:0" onclick="slideshow.step(+1)">»</div>
    </div>
    <script>slideshow.start();</script>
    </p>
