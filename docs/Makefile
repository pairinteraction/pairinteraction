# Makefile for Linux and MacOS
# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
ALLSPHINXOPTS = "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help $(ALLSPHINXOPTS)

.PHONY: help Makefile clean new

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ $(ALLSPHINXOPTS)


# Additional custom targets
#
clean:
	rm -rf $(SOURCEDIR)/_autosummary
	@$(SPHINXBUILD) -M clean $(ALLSPHINXOPTS)

livehtml:
	sphinx-autobuild -b html $(ALLSPHINXOPTS) --watch ../src/pairinteraction/
