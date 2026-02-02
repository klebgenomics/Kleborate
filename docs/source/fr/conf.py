# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'KleborateModular'
copyright = '2026, Sylvain Brisse, Mary Maranga, Kathryn Holt'
author = 'Sylvain Brisse (sylvain.brisse@pasteur.fr)'
release = '3.0.0'
version = '3.0.0'

# Set language to French
language = 'fr'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_theme_options = {
    "show_toc_level": 2,
    "home_page_in_toc": True,
    "navigation_depth": 4,
    "navbar_end": ["langswitch.html"],

}

templates_path = ['../_templates']

html_static_path = html_static_path = ['../_static']
html_logo = '../_static/logo.png'
html_title = "Kleborate (FR)"

html_css_files = [
    'langswitch.css',
]
