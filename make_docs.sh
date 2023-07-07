set -e

# Set up docs folder
sphinx-quickstart --sep --dot _ --language en --suffix .rst --master index --makefile --batchfile --use-make-mode --author "Luca Fiorito" --project "SANDY API" -v 1 --release 0 --ext-autodoc --ext-doctest --ext-githubpages --ext-mathjax --extensions numpydoc api_docs

# Change configuration file
sed -i '13i\   sandy' api_docs/source/index.rst
sed -i "s/alabaster/sphinx_rtd_theme/" api_docs/source/conf.py
# add folder with sources
sed -i '1 i\import os\nimport sys\nsys.path.insert(0, os.path.abspath(os.path.join("..", "..", "sandy")))\nprint(sys.path)' api_docs/source/conf.py


# Create rst files
sphinx-apidoc --separate --force --module-first -o api_docs/source sandy

# Run make
[[ -n $1 ]] && (cd api_docs && make $1)
