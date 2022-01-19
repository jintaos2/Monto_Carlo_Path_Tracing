# -*- coding: utf-8 -*-

import sphinx_rtd_theme

project = u'Monto_Carlo_Path_Tracing'
release = '1.0.1'
author = u'jintaos2'
copyright = author
language = None

extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
    'myst_parser',
    'sphinx.ext.githubpages'
]
myst_enable_extensions = ["dollarmath","html_image"]
source_suffix = ['.rst', '.md']

templates_path = ['_templates']
exclude_patterns = []
gettext_compact = False

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'navigation_depth': 5,
    'prev_next_buttons_location': 'bottom'
}
html_title = 'Monto_Carlo_Path_Tracing'





