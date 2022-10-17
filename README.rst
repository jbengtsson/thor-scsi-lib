Thor-SCSI – Self-Consistent Symplectic Integrator for Charged Particle Beam Dynamics
------------------------------------------------------------------------------------

Requirements
------------
* (GNU compatible) C/C++ compiler
* GNU autoconf/automake environment and libtool.
* GNU Scientific Library (GSL): https://www.gnu.org/software/gsl.
* Armadillo (for linear algebra): http://arma.sourceforge.net.
* Python https://www.python.org/ for the python interface

The library uses the range checking inmplementation of e.g. `std::vector` as provided by GNU C++; thus its dependency on the GNU compiler collections.

To install
----------

Setup of repository
===================

Dowload the repository and checkout the proper branch. Here it's assumed you will use the directoy `git_repos/tracy-3.5` in your home directory for the tracy code tree.

For this use the following commands to create the directoy `git_repos` and to clone the tree into the tracy-3.5 directory.

.. code:: shell

   mkdir git_repos
   cd git_repos
   git clone git@github.com:jbengtsson/Thor_scsi.git
   cd Thor_scsi

Then select the proper tree by

.. code:: shell

   git checkout Thor_scsi



C++ library
===========

First create environment variable $TRACY_LIB. This will be the prefix where the built library and include files will be installed later on e.g:

.. code:: shell

   export TRACY_LIB=$HOME/git_repos/Thor_scsi


To build the library use:

.. code:: shell

   cd Thor_scsi
   (g)libtoolize
   ./bootstrap
   ./configure --prefix=$Thor_scsi_LIB
   make
   make install

Please note: using the dynamic library in non standard location will require proper set up of the environment later on (e.g. adding the directory where the library is located to `LD_LIBRARY_PATH` environment variable).


Python interface
================

The python interface is based on https://github.com/pybind/pybind11. Building this interface requires to select the proper directory

.. code:: shell

  cd git_repos
  cd Thor_scsi/python

Install proper dependencies

.. code:: shell

    pip3 install -r requirements.txt


And build the extension e.g.

.. code:: shell

    python3 setup.py build
    python3 setup.py install

For further details of the build system see https://pypi.org/project/setuptools/


To run the regression tests
---------------------------

All regression tests can be run using

.. code:: shell

    pip3 install nose
    python3 setup.py nosetests

To run the demo/test program
----------------------------


.. code:: shell

    python3 examples/tst.py


References
==========

* Python interface:

  Guidelines & automated regression testing bootstrapped by Pierre Schnizer:
  
    P\. Schnizer, W\. Khail, J\. Bengtsson *Small Talk on AT* `IPAC 2022 TUPOST029.`_
    
    .. _`IPAC 2022 TUPOST029.`: http://dx.doi.org/10.18429/JACoW-IPAC2022-TUPOST029

  Initial demo/prototype & guidelines by Jan Chrin, PSI, 2017:
  
    J\. Chrin *Channel Access from Cython (and other Cython use cases)* `EPICS Collaboration Meeting 2017.`_
  
    .. _`EPICS Collaboration Meeting 2017.`: https://indico.esss.lu.se/event/889/contributions/7038/attachments/6800/9762/Cython_EpicsTM_Oct2017_Barcelona.pdf#page=32

* Model Server:

    P\. Schnizer, J\. Bengtsson, T\. Birke, L\. Ramirez *Online Model Developments for BESSY II and MLS* `IPAC 2021 WEPAB317.`_
    
    .. _`IPAC 2021 WEPAB317.`: http://dx.doi.org/10.18429/JACoW-IPAC2021-WEPAB317

* Use Cases:

    J\. Bengtsson, T\. Nicholls, W\. Rogers *A CAD Tool for Linear Optics Design: A Use Case Approach* `IPAC 2021 MOPAB047.`_
    
    .. _`IPAC 2021 MOPAB047.`: http://dx.doi.org/10.18429/JACoW-IPAC2021-MOPAB047
    
    J\. Bengtsson, M\. Davidsaver *An Accelerator Physics - Software Engineering Collaboration* `EPICS Collaboration Meeting EPICS 2016.`_
    
    .. _`EPICS Collaboration Meeting EPICS 2016.`: https://indico.esss.lu.se/event/507/contributions/3830

* Beam Dynamics Model:

    J\. Bengtsson, W\. Rogers, T\. Nicholls *A CAD Tool for Linear Optics Design: A Controls Engineer's Geometric Approach to Hill's Equation* `10.48550/arXiv.2109.15066 (2021).`_
    
    .. _`10.48550/arXiv.2109.15066 (2021).`: http://dx.doi.org/10.48550/arXiv.2109.15066
    
    J\. Bengtsson, D\. Briggs, G\. Portmann *A Linear Control Theory Analysis of Transverse Coherent Bunch Instabilities Feedback Systems (The Control Theory Approach to Hill's Equation)* `CBP Tech Note-026, PEP-II AP Note 28-93 (1993).`_
    
    .. _`CBP Tech Note-026, PEP-II AP Note 28-93 (1993).`: https://escholarship.org/uc/item/64s937sf


Thor
====

Author: Johan Bengtsson

Having implemented *DA-Pascal* in the early 1990s, see ref. below, based on a *recursive* approach – i.e., *automatable* by an *universal Turing machine* (any modern *digital computer*) – for *nonlinear beam dynamics analysis* by utilizing *Lie series* on a *beam line object*; *automated* to arbitrary order by  *Truncated Power Serias Algegra* (TPSA). In particular, implemented as a *Pascal module/software library* by extending the *standard procedures & functions* for N. Wirth's *Pascal-S compiler/interpreter*. Hence, in 1994, rather than participating in a "talkshop" for the *CLASSIC collaboration*, we instead prototyped & implemented a C++ *beam line class* based on a *polymorphic number object with reference counting*. For which the latter mechism provides for *garbage collection*; since it is not provided by C++, vs. e.g. *Lisp* & *Smalltalk*.

To quote Forest in:

  E\. Forest, F. Schmidt, E. McIntosh *Introduction to the Polymorphic Tracking Code* `CERN-SL-2002-044 (AP), KEK-Report 2002-3 (2002).`_

  .. _`CERN-SL-2002-044 (AP), KEK-Report 2002-3 (2002).`: https://cds.cern.ch/record/573082/files/CERN-SL-2002-044-AP.pdf

  *Therefore the proper implementation of a fibre bundle was a sine qua non, non-negotiable point, which
  Forest and Bengtsson insisted upon in the early days of pre-CLASSIC C++ collaborations. A quick look
  at CLASSIC (MAD9) shows that the CLASSIC structure does not satisfy this condition. In particular,
  patches are derived from the beam element class and are thus of the same nature as the element. Patches
  are generally a figment of one’s mathematical imagination, useful tools, but they are not physical elements
  which can be ordered from a factory. We believe this is a major flaw in the CLASSIC design. It is perhaps
  the result of placing too much emphasis on implementation using C++ capabilities, rather than the basic
  mathematical framework. We believe this accounts for the excessive number of classes and the complexity
  of CLASSIC (MAD9).*

  ...

  *Here, in PTC, as well as in the original C++ classes that Bengtsson dreamt up in collaboration with
  Forest, the geometrical nature of the fibre reigns supreme. The next step is to introduce a magnet, namely
  EL and/or ELP. The propagator of the full fibre, if well-defined, inherits properties from the chart itself.
  In other words a magnet exists first as a piece of material junk. It can be rotated, translated and drawn.
  The chart provides the connection between this magnet/junk and the external three dimensional space.
  Obviously this exists independently of the existence of single particle propagators associated to EL and/or
  ELP. It is a remarkable mathematical feature that these propagators, under certain conditions, inherit the
  transformational properties of the chart. Of course PTC is set up to take advantage of this.*
 
  ...

  *TRACYII was based on the belief that a dumb user interface should be built on the foundation of a smart
  user interface. In this way complex situations could always be handled. This was so successful that, in the
  2 years of the PEPB design (SLAC), Robin and Bengtsson recompiled TRACYII no more than 2 or 3 times.*

  ...

  *In the case of TRACYII, this was realized by separating the lattice input file (dumb user) from the
  command input file (smart user). This idea, originally from Nishimura, was turned into an uncompromising
  product by Bengtsson. In PTC the same can be achieved by stripping all the core routines from any dumb
  user idiosyncracies. One example common to TRACYII and PTC is the absence of quadrupoles in the core.* 

  ...

  *In addition, as we shall see, if some user’s algorithm uses PTC extended definition
  of the ray to compute the equivalent of the “synchrotron integrals,” then it will be correctly computed under
  any possible mispowering and misaligning of the elements. PTC is a faithful representation of a part of
  nature, just as Seurat’s painting is a faithful representation of some aspect of a scene. In addition, just as
  pointillism adds to the natural setting a seemingly unnatural element, PTC adds properties to the ray being
  tracked which do not exist in nature. In the case of PTC, thanks to a polymorphic type first dreamt up by
  Bengtsson, the electron carries with itself a potential Taylor Series whose variable space is nearly infinite.*

  ...

  *My views have been, at least since the C++ business got underway, that the flow through the magnet
  must be elevated to the status of a mathematical object. And then, it must find its counterpart on the
  silicon canvas, whether painted in C++ or any other language. Polymorphism, Bengtsson’s pointillism, will
  take care of the rest. This is achieved by a local “s” -dependent theory which is shaped around individual
  magnets. The global system is then patched together. The mathematicians gave us the tools to manipulate
  this object: the fibre bundle. PTC simply creates a restricted fibre bundle on the computer, one which is
  relevant to particle accelerators. This structure is incompatible with standard Courant-Snyder theory and
  other similar constructs like Sand’s integrals.*
 
  ...

  *Besides the two individuals whose names appear on this paper and Aimin Xiao who collaborated on the very
  first prototype, I would like to thank Johan Bengtsson (of parts unknown) for convincing me that, at least
  in C++, one could go ahead and make a reasonable job of polymorphism and fibre bundles.*

  .. image:: images/seurat.png

I.e., eventually, he re-implemented the strategy/approach in *Fortran-90*; since by then that *grammar* had been
extended to support *operator overloading*.

Tracy-2
=======

The symplectic integrator for realistic modeling of magnetic lattices for ring-based synchrotrons was initially implemented as a *Pascal module/beam dynamics software library*, by the author 1990, as an *on-line model* to guide the ALS commissioning. In particular, care was taken for the software architecture & resulting records/modules – akin to *objects* although not explicitly supported by the artificial language grammar – to reflect the structure of the mathematical objects describing the underlying *beam dynamics model*.

Hence, the code was also benchmarked & calibrated as part of the ALS commissioning:

  J\. Bengtsson, M. Meddahi *Modeling of Beam Dynamics and Comparison with Measurements for the Advanced Light Source (ALS)* `EPAC 1994.`_

  .. _`EPAC 1994.`: https://accelconf.web.cern.ch/e94/PDF/EPAC1994_1021.PDF

Malika Meddahi was a postdoc at ALS who did her thesis at CERN (for which she received the *Prix Daniel Guinier* in France that year):

  *Effets faisceau-faisceau dans le collisionneur protons-antiprotons du SPS* `CERN-SL-91-30-BI (1991)`_

  .. _`CERN-SL-91-30-BI (1991)`: http://cds.cern.ch/record/223301)

The resulting C code, see below, has now been re-factored by introducing a C++ *beam line class*; i.e., to recover the transparency & simplicity of the original *beam dynamics model*.

Remark: Although the entire *beam dynamics model* had to be replaced & the model/code/"approach" re-architectured & structured – for a reusable approach – as a *Pascal beam dynamics libary* (standard practise in software engineering), the code was named *Tracy-2*, i.e., inspired by the, somewhat archaic demo/prototype/concept *Tracy*:

(E.g. by not having figured out how to pass records as function/procedure variables – vs. scalars only – for the beam dynamics library, the API was rather poor/sloppy; i.e., not scalable/ill suitable to cope with complexity)

  H\. Nishimura "*TRACY, A Tool for Accelerator Design and Analysis*" `EPAC 1988.`_

  .. _`EPAC 1988.`: https://accelconf.web.cern.ch/e88/PDF/EPAC1988_0803.PDF

  for which the *beam dynamics model* was based on the *linearized quadratic Hamiltonian*:

  .. image:: images/H_2.png

for *linear optics design*; i.e., for a *bare lattice* with *mid-plane symmetry*.

E.g. by not having figured out how to pass records as function/procedure variables – vs. scalars only – for the beam dynamics library, the API was rather poor/sloppy; i.e., not scalable/ill suited to cope with the complexity of a realistic model. As expressed by Forest in the title of:

  E\. Forest *A Hamiltonian-Free Description of Single Particle Dynamics for Hopelessly Complex Periodic Systems* `J. Math. Phys. 31 (1990).`_
  
  .. _`J. Math. Phys. 31 (1990).`: http://dx.doi.org/10.1063/1.528795%7D

Hence, the one thing we did find useful for a realistic on-line model – having already implemented one as a sci fellow for LEAR, CERN, in the late 1980s and before that having worked as a teaching assistent at the dept. Software Engineering, Lund Inst. of Tech. while pursuing a MsSci EE – and, hence, adopted for the ALS on-line model. Was the implementation of the beam dynamics model as an extension of the *standard procedures & functions* for the *Pascal-S compiler/interpreter* by N. Wirth; architected as a Pascal software library/module:
  
  M\. Rees, D\. Robson *Practical Compiling with Pascal-S* `(Addison-Wesley, 1988).`_
  
  .. _`(Addison-Wesley, 1988).`: https://books.google.com/books?id=hLomAAAAMAAJ
  
  S\. Pemberton, M\. Daniels *The P4 Compiler and Interpreter* `(1982).`_
  
  .. _`(1982).`: https://homepages.cwi.nl/~steven/pascal/book/pascalimplementation.html

  N\. Wirth *PASCAL-S: A Subset and its Implementation* `Institut für Informatik, ETH, Zürich (1975).`_

  .. _`Institut für Informatik, ETH, Zürich (1975).`: http://pascal.hansotten.com/uploads/pascals/PASCAL-S%20A%20subset%20and%20its%20Implementation%20012.pdf
  
  *Pascal-P6* https://sourceforge.net/projects/pascal-p6.

In other words, since 1994 our *toolkit* – althout it based on one model: the *Hamiltonian for a charged particle in an external electromagnetic field* & a *symplectic intrator* for *magnetic multipoles* & *insertion devices* for ditto – it was implemented as two different codes: Tracy-2 & Thor. Hence, eventually, these were consolidated by using C++ *templates* for the *polymorphich number object* and *beam line class*; aka Tracy-2,3.

Contributions
-------------
* The symplectic integrator for *RADIA kick maps*:

    P\. Elleaume *A New Approach to the Electron Beam Dynamics in Undulators and Wigglers* `EPAC 1992.`_

    .. _`EPAC 1992.`: https://accelconf.web.cern.ch/e92/PDF/EPAC1992_0661.PDF

  was implemented by Laurent Nadolski, SOLEIL, 2002.

* The original •Pascal library/code• was machine translated to C and re-used to implement a *model server* for the SLS commissioning:

    M\. Böge *Update on TRACY-2 Documentation* `SLS Tech Note SLS-TME-TA-1999-0002 (1999).`_

    .. _`SLS Tech Note SLS-TME-TA-1999-0002 (1999).`: http://ados.web.psi.ch/slsnotes/tmeta9902.pdf

    M\. Böge, J. Chrin *A CORBA Based Client-Server Model for Beam Dynamics Applications* `ICALEPCS 1999.`_

    .. _`ICALEPCS 1999.`: https://accelconf.web.cern.ch/ica99/papers/mc1p61.pdf

  with `p2c.`_

    .. _`p2c.`: http://users.fred.net/tds/lab/p2c/historic/daves.index-2012Jul25-20-44-55.html

* Similarly, James Rowland re-used the C version to implement a *Virtual Accelerator* interfaced to EPICS as a *Virtual Input Output Controller* (VIOC):

    M\. Heron, J. Rowland, et al *Progress on the Implementation of the DIAMOND Control System* `ICALEPCS 2005.`_

    .. _`ICALEPCS 2005.`: https://accelconf.web.cern.ch/ica05/proceed-ings/pdf/P1_018.pdf

* Besides, the internal *numerical engine* was manually translated to C and re-used for:

    A\. Terebilo *Accelerator Toolbox for MATLAB* `SLAC-PUB-8732 (2001).`_
  
    .. _`SLAC-PUB-8732 (2001).`: http://www-public.slac.stanford.edu/sciDoc/docMeta.aspx?slacPubNumber=SLAC-PUB-8732
