# cMoP_tmp
Self-Consistent Mori Projector Solver

## List of Contents

1. Introduction
2. Method
   1. Problem
   2. Simplifications
3. Code Layout
   1. dev
      1. The ``rho`` cell array
      2. The ``input`` struct
      3. Setup functions
      4. ``CreateInteractions``
      5. ``TimeIter``
      6. ``CalculateSteadyState``
   2. bin
   3. test
   4. utils
4. References

## 1. Introduction

cMoP_tmp provides a set of Matlab functions which solves many-body driven-dissipative quantum systems 
using the self-consistent Mori projector method [[1]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.245108).  The 
code is designed to solve time-dependent and steady-state problems and can be used for arbitrary
onsite Liouvillian and any bilinear interaction.  This code was written by Edmund Owen at Heriot-Watt University.

## 2. Method

### 2.i Problem

cMoP_tmp is designed to calculate the reduced density matrix
<a href="https://www.codecogs.com/eqnedit.php?latex=\rho_0" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\rho_0" title="\rho_0" /></a>
of a many-body system where the equation of motion is given by the Lindblad equation

<a href="https://www.codecogs.com/eqnedit.php?latex=\dot{R}(t)&space;=&space;(\mathcal{L}_0&space;&plus;&space;\mathcal{L}_I)&space;R(t)" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\dot{R}(t)&space;=&space;(\mathcal{L}_0&space;&plus;&space;\mathcal{L}_I)&space;R(t)" title="\dot{R}(t) = (\mathcal{L}_0 + \mathcal{L}_I) R(t)" /></a>
<br /> <br />
 
where 
<a href="https://www.codecogs.com/eqnedit.php?latex=R(t)" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?R(t)" title="R(t)" /></a>
is the density matrix of the total system, 
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{L}_0" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\mathcal{L}_0" title="\mathcal{L}_0" /></a>
is the tensor product of on-site Liouvillians (i.e. the Liouvillian for the non-interacting system) and
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{L}_I" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\mathcal{L}_I" title="\mathcal{L}_I" /></a>
is the Liouvillian which describes bilinear interactions between the sites:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{L}_I&space;=&space;\sum_{<i,j>}&space;A_i&space;\cdot&space;B_j" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\mathcal{L}_I&space;=&space;\sum_{<i,j>}&space;A_i&space;\cdot&space;B_j" title="\mathcal{L}_I = \sum_{<i,j>} A_i \cdot B_j" /></a>
<br /><br />

where
<a href="https://www.codecogs.com/eqnedit.php?latex=A_i" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?A_i" title="A_i" /></a>
and
<a href="https://www.codecogs.com/eqnedit.php?latex=B_j" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?B_j" title="B_j" /></a>
are single site operators operating on sites 
<a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?i" title="i" /></a>
and
<a href="https://www.codecogs.com/eqnedit.php?latex=j" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?j" title="j" /></a>
respectively.  The cMoP method manipulates the Lindblad equation of motion into a form which describes
the equation of motion of the reduced density matrix
<a href="https://www.codecogs.com/eqnedit.php?latex=\rho_0&space;(t)&space;=&space;\mathrm{Tr}_{\o}&space;\{&space;R&space;(t)&space;\}" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\rho_0&space;(t)&space;=&space;\mathrm{Tr}_{\o}&space;\{&space;R&space;(t)&space;\}" title="\rho_0 (t) = \mathrm{Tr}_{\o} \{ R (t) \}" /></a>
where 
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathrm{Tr}_{\o}&space;\{&space;\cdot&space;\}" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\mathrm{Tr}_{\o}&space;\{&space;\cdot&space;\}" title="\mathrm{Tr}_{\o} \{ \cdot \}" /></a>
is a trace over all sites other than the site of interest 0.  The resulting equations are non-linear
and non-Markovian and require simplifications in order to acquire a numerical solution.

### 2.ii Simplifications

The cMoP_tmp code solves the equation

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;\rho_0&space;(t)}{\partial&space;t}&space;=&space;\sum_i&space;\mathcal{L}_i&space;\,&space;\rho_0&space;(t)" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\frac{\partial&space;\rho_0&space;(t)}{\partial&space;t}&space;=&space;\sum_i&space;\mathcal{L}_i&space;\,&space;\rho_0&space;(t)" title="\frac{\partial \rho_0 (t)}{\partial t} = \sum_i \mathcal{L}_i \, \rho_0 (t)" /></a>
<br /><br />

where 
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{L}_i" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\mathcal{L}_i" title="\mathcal{L}_i" /></a>
can be any non-linear superoperator.  At the moment, only onsite, mean-field and Born
operators are implemented (see Ref. [[1]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.245108)).  Time-independent 
versions of these superoperators also exist for steady state calculations.

The code is designed for finite lattices where translational symmetry can be used to reduce the number of 
degrees of freedom.  In such systems, the reduced density matrices of each site is assumed to be the
same.  However, the code does include a feature for partitioning a system into a finite number of different
density matrices which can interact with each other (e.g. if the system has an A and B sublattice then
the density matrices for each sublattice can be calculated simulateously with interactions coupling the
A sublattice to the B sublattice and vice versa, see Ref [[2]](http://iopscience.iop.org/article/10.1088/1367-2630/aab7d3/meta)).

## 3. Layout

In this section, I will attempt to describe the layout of the code along with design choices and
general explanations of what the various functions do.

### 3.i dev

This folder contains all of the functions which the cMoP code needs in order to run.

For a first-time user; the most important functions are the ``Setup``*name*, ``CreateInteractions``, 
``TimeIter`` and ``CalculateSteadyState`` functions.  In addition, the code uses a struct which is 
labelled ``input`` in the examples which contains the input parameters of the problem, and a struct
called ``solution`` which is an auxilliary struct used to contain temporary copies of the density matrix.

#### 3.i.a The ``rho`` cell array

The density matrix is stored as a cell array.  This allows multiple density matrices to be stored
in the same variable when the system is partitioned into two or more different regions.  Note
that ``rho`` must still be a (1x1) cell array even when there is only one partition

#### 3.i.b The ``input`` struct

As the inputs may be applied to different partitions differently, ``input`` is actually an array
containing a ``subinput`` field which is a cell array of ``input`` structs.  This is a bit of a
mess right now as ``input`` does contain other fields even though it probably shouldn't.  *Be
careful to make sure that interactions are system parameters are in the ``subinput`` field!*

#### 3.i.c Setup functions

The ``Setup``*name* functions (e.g. ``SetupSystem``) take the ``input`` struct and a cell array
of string-value pairs and initialise the system, returning the new version of ``input``.  These
functions will set defaults.  Available input parameters are not documented but should be relatively
clear from context and the parameter parser contained within the functions.

#### 3.i.d ``CreateInteractions``

This function creates interactions between partitions and between copies of the density matrix and
itself (in the case of translationally invariant systems).  The inputs for this funciton are the
``input`` struct and a ``Map`` containing key-value pairs.  The necessary key set is

``keySet = {'interactionStrength', 'A', 'B', 'correlations'}``

although there are additional options.  ``interactionStrength`` is the strength of the interaction,
``A`` and ``B`` are the two parts of the bilinear interaction and ``correlations`` is a correlation
matrix which determines which pairs of interactions are summed over in the Born term.  The bilinear
operators are structs with the fields ``Index``, which is the partition index which the operator
acts on, ``SiteOperator`` which contains the operator for a given site and ``SiteLabel`` which is
the site in the given partition that the operator acts on.

NOTE: *the interactions should be put in the subinput and not in input*! (see bin/AnisotropicPartition.m)

#### 3.i.e ``TimeIter``

Time-integrates the cMoP equations of motion based on the parameters defined within the ``input`` 
struct for the starting density matrix ``rho``.  The outputs are defined by the variable ``probelist``
which is set using the ``SetupTimeIter`` function.

#### 3.i.f ``CalculateSteadyState``

Calculates steady state solutions for the cMoP equations by finding the zero of the non-linear
equation:

<a href="https://www.codecogs.com/eqnedit.php?latex=\sum_i&space;\mathcal{L}_i&space;\,&space;\rho_0&space;=&space;0" target="_blank">
<img src="https://latex.codecogs.com/gif.latex?\sum_i&space;\mathcal{L}_i&space;\,&space;\rho_0&space;=&space;0" title="\sum_i \mathcal{L}_i \, \rho_0 = 0" /></a>
<br />

see Ref. [[1]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.245108) for details.
Note that this function needs the auxilliary struct ``solution`` which is changed within
``solution`` but will not be changed after ``CalculateSteadyState`` has run.

### 3.ii bin

This folder contains complete solutions for some example problems:

Anisotropic_Partition: This file reproduces the results of Ref. [[2]](http://iopscience.iop.org/article/10.1088/1367-2630/aab7d3/meta) 
by calculating cMoP results for the anisotropic Heisenberg lattice for lattices with any number of
dimensions.  This problem demonstrates the cMoP method and includes partitioning onto multiple different
lattices with different reduced density matrices.

OneDOpticalLattice: This file produces mean-field results for a one-dimensional optical lattice.
This problem demonstrates a system which contains both unitary and dissipative interactions for
a mean-field interaction Liouvillian.

### 3.iii test

This folder contains Matlab tests to test various functions of the cMoP_tmp code.  To run these tests,
move to the test folder and enter the runtests command.

### 3.iv utils

Some utility functions written by other parties.  Non-essential.

## 4. References

[[1]](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.245108) Original paper detailing the cMoP method

[[2]](http://iopscience.iop.org/article/10.1088/1367-2630/aab7d3/meta) Limit cycles and cMoP calculations for a partitioned lattice
