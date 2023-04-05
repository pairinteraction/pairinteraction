# Makefile for Sphinx documentation

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
    @$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
    @$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Livehtml
SPHINXAUTOBUILD = sphinx-autobuild

ALLSPHINXLIVEOPTS = $(ALLSPHINXOPTS) -q \
   --port 0 \
   --host 0.0.0.0 \
   --open-browser \
   --delay 1 \
   --ignore "*.swp" \
   --ignore "*.pdf" \
   --ignore "*.log" \
   --ignore "*.out" \
   --ignore "*.toc" \
   --ignore "*.aux" \
   --ignore "*.idx" \
   --ignore "*.ind" \
   --ignore "*.ilg" \
   --ignore "*.tex" \
   --watch ../.

.PHONY: livehtml
livehtml:
    $(SPHINXAUTOBUILD) -b html $(ALLSPHINXLIVEOPTS) . $(BUILDDIR)
    @echo
    @echo "Build finished. The HTML pages are in $(BUILDDIR)."
