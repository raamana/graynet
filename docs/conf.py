from pathlib import Path
import sys

DOCS_DIR = Path(__file__).resolve().parent
REPO_ROOT = DOCS_DIR.parent
sys.path.insert(0, str(REPO_ROOT))

import graynet
project = 'graynet'
copyright = '2017, Pradeep Reddy Raamana'
author = 'Pradeep Reddy Raamana'
version = graynet.__version__
release = graynet.__version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinxarg.ext',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
html_theme = 'alabaster'
html_static_path = ['_static']
htmlhelp_basename = 'graynetdoc'

latex_elements = {
}

latex_documents = [
    (master_doc, 'graynet.tex', 'graynet Documentation',
     'Pradeep Reddy Raamana', 'manual'),
]

man_pages = [
    (master_doc, 'graynet', 'graynet Documentation',
     [author], 1)
]

texinfo_documents = [
    (master_doc, 'graynet', 'graynet Documentation',
     author, 'graynet', 'One line description of project.',
     'Miscellaneous'),
]

intersphinx_mapping = {'python': ('https://docs.python.org/3', None)}
