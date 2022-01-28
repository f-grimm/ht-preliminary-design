
Ressources
==========

This site was created using `Sphinx <https://www.sphinx-doc.org/en/master/>`_, a tool that can generate the documentation for Python modules automatically based on the provided doc-strings. The following page describes how to add information or create a similar documentation from scratch.

Setup
-----

::

    pip install sphinx
    pip install sphinx_rtd_theme


In the *docs* directory:

#. Initialize Sphinx. (Use default values by repeatedly pressing ``Enter``)
   
   ::
   
   	   sphinx-quickstart

#. Edit **conf.py** to set up directories, extentions, themes, etc.

#. Create the reST files, one for each module and a combined page **modules.rst**
   
   ::
   
       sphinx-apidoc -o <output_folder> <py_folder>
    
#. Edit **index.rst** to include the desired pages


Maintain
--------

#. Edit index, modules, or add new pages (in the reStructured-Text format)
  
#. Build HTML and PDF files  

   In the *docs* directory:
   
   ::
   
       make html
       make latexpdf


Guides
------

* `Getting started with Sphinx <https://docs.readthedocs.io/en/stable/intro/getting-started-with-sphinx.html>`_
* `reStructured-Text <https://draft-edx-style-guide.readthedocs.io/en/latest/ExampleRSTFile.html>`_