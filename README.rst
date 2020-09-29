NuSTAR General Utilities
========================


|Readthedocs| |Astropy|



Repository for (hopefully) useful NuSTAR python tools
--------------------------------------

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Nuclear Spectroscopic Telescope ARray
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. image:: https://www.nustar.caltech.edu/system/avm_image_sqls/binaries/26/page/nustar_artistconcept_2.jpg?1393022433
    :target: http://www.nustar.caltech.edu
    :alt: NuSTAR



Installation
------------

We require at least python 3.7 for this installation as well as the other dependencies
listed in requirements.txt.

To install these via conda, do:

.. code-block:: bash

    conda install --yes --file requirements.txt

In the current version, a development installation is recommended.
From the shell:

.. code-block:: bash

    pip install -e .

Documentation
-------------
To build the documentation, from the shell:

.. code-block:: bash

    sphinx-build docs docs/_build

And then open the docs/_build/index.html in your browser.

Or, see the readthedocs documentation: |Readthedocs|


Contributing
------------

We love contributions! nustar-gen-utils is open source,
built on open source, and we'd love to have you hang out in our community.

If you'd like to say hello, drop onto the NuSTAR Slack channel from the
`NuSTAR Observer's page <https://www.nustar.caltech.edu/page/observers>`_ .

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


**Imposter syndrome disclaimer**: We want your help. No, really.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There may be a little voice inside your head that is telling you that you're not
ready to be an open source contributor; that your skills aren't nearly good
enough to contribute. What could you possibly offer a project like this one?

We assure you - the little voice in your head is wrong. If you can write code at
all, you can contribute code to open source. Contributing to open source
projects is a fantastic way to advance one's coding skills. Writing perfect code
isn't the measure of a good developer (that would disqualify all of us!); it's
trying to create something, making mistakes, and learning from those
mistakes. That's how we all improve, and we are happy to help others learn.

Being an open source contributor doesn't just mean writing code, either. You can
help out by writing documentation, tests, or even giving feedback about the
project (and yes - that includes giving feedback about the contribution
process). Some of these contributions may be the most valuable to the project as
a whole, because you're coming to the project with fresh eyes, so you can see
the errors and assumptions that seasoned contributors have glossed over.

Note: This disclaimer was originally written by
`Adrienne Lowe <https://github.com/adriennefriend>`_ for a
`PyCon talk <https://www.youtube.com/watch?v=6Uj746j9Heo>`_ , and was adapted by
nustar-gen-util based on its use in the README file for the
`MetPy project <https://github.com/Unidata/MetPy>`_.



.. |Readthedocs| image:: https://img.shields.io/badge/docs-latest-brightgreen.svg?style=flat
    :target: https://nustar-gen-utils.readthedocs.io/en/latest/
    
.. |Astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org
    :alt: Powered by Astropy Badge
