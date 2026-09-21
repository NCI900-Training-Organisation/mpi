# sphinx-book-theme Template

[Example site](https://test-template-2026.readthedocs.io/en/latest/) built with `sphinx-book-theme` branch.  

This template uses [sphinx-book-theme](https://sphinx-book-theme.readthedocs.io/en). Some useful pages from the theme documentation are:  

- [Theme Config Options](https://sphinx-book-theme.readthedocs.io/en/stable/reference.html)  
- [Math, Code, Admonition Elements](https://sphinx-book-theme.readthedocs.io/en/stable/content/pydata-content-blocks.html)
- [Example Pages](https://sphinx-book-theme.readthedocs.io/en/stable/reference/kitchen-sink/index.html)


The theme is built upon [PyData Sphinx Theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/index.html). Hence options in `html_theme_options` from the parent theme are also available. 

## Configuration Files
The boilerplate files for the sphinx theme include:   
- `docs/requirements.txt` : Tell Readthedocs site to install the packages needed for the build.  
- `docs/source/conf.py` : Configrations for Sphinx to build HTML files. Theme options are also in this file.
- `docs/source/_static` : Storing logos and files needed for the static site.
- `docs/source/_template` : If needed, any template customisation files to use.  
- `docs/source/index.rst`: Update the homepage per project.


## Use This Template 
1. Use this template to create a new repository and clone to local to make changes.
2. Update `conf.py` where marked as `TODO` in comments.
3. Add content files under `docs/source/`. Content can be `rst`, `md` or `ipynb` files, thanks to the `myst-nb` extension. 

> ℹ️ **Info**  
> By default, this theme renders Markdown as [MyST flavor](https://myst-parser.readthedocs.io/en/latest/live-preview.html), which combines some features from rst and markdown. To change flavor, change `nb_render_markdown_format` value in `conf.py` file.

4. After making changes, use [Sphinx](https://www.sphinx-doc.org/en/master/usage/installation.html) to build the docs locally for preview.  
`pip install sphinx-book-theme`  
`pip install -U sphinx`  
`sphinx-build -M html <sourcedir> <outputdir>`

> ℹ️ **Info**  
> 
> By default, the output HTML will be available in the `_build/html/` directory after running a Sphinx build. You can view it locally by opening `index.html` in a browser.
> 
> For more on Sphinx commands and usage, see the [Sphinx documentation](https://www.sphinx-doc.org/en/master/usage/index.html).

5. When you are happy with the build, push changed files under `docs/source` to the repo.
6. Login [Readthedocs](https://app.readthedocs.org/dashboard/) with team account and build a new project using the repo.

-----
## Useful Links
- [reStructuredText in Sphinx](https://app.readthedocs.org/dashboard/)
- [sphinx-book-theme](https://sphinx-book-theme.readthedocs.io/en/stable/index.html)
- [PyData Theme](https://pydata-sphinx-theme.readthedocs.io/en/stable/user_guide/index.html)
-[MyST](https://myst-parser.readthedocs.io/en/latest/syntax/typography.html)


**Common ARE errors** : https://opus.nci.org.au/spaces/NCIT/pages/399803051/Common+ARE+Errors
