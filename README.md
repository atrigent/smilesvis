# SMILES string visualizer

This is a fairly simple utility that takes a SMILES string as input and
attempts to create a reasonable 3D visualization of the compound being
described. To do so, it uses the following libraries:

* pyparsing, for parsing SMILES
* pyopengl, for OpenGL APIs (including GLU and GLUT)
* numpy, for linear algebra tasks

SMILES is a simple way of describing a compound by describing its constituents
and how they are connected. Unfortunately, without implementing all of quantum
mechanics and doing some very expensive computations, it is impossible to
determine the 3D structure of a molecule using only a SMILES string. Therefore,
this program uses VSEPR theory, a simple set of rules for determining 3D
structure that works for many common classes of molecules. VSEPR theory is
often used in the teaching of chemistry and works especially well with organic
compounds.

One particular complication in the determination of 3D structure is the
handling of single bonds. The relative orientations of the two atoms in a
single bond can change freely because very little energy is required to twist
these bonds. Therefore, compounds with single bonds can exist in many different
"conformations" depending on the relative orientations of every two atoms
participating in these bonds. The number of conformations increases very
quickly as the number of single bonds increases. An extreme manifestation of
this is in the protein folding problem, where expensive computations are needed
to compute the geometries of very large molecules with many single bonds. The
same problem, at a smaller scale, is also relevant to this program. For
example, cyclic compounds can only exist in certain conformations.

As noted above, the program uses OpenGL for 3D rendering. However, the specific
methods that it uses are somewhat outdated. Its use of GLU and GLUT in
particular make it a bad representation of current OpenGL development
practices. These older, simpler methods were chosen because the newer, more
flexible methods were not required.

In order for all of the elements of the 3D rendering to be drawn correctly
relative to each other, these elements have to be rotated in the appropriate
ways. One problem was discovering how to rotate the 3D representation of a bond
such that it would be pointing along a bond vector. Another problem involved
rotating the 3D representation of an atom such that it would meet up with an
already drawn atom. This problem was modeled as rotating a bond vector such
that it would be drawn pointing in the opposite direction of an existing bond
vector. The last major problem concerned the correct drawing of double bonds.
The two atoms involved in a double bond will not rotate relative to each other
like atoms in a single bond, so the two double bonded atoms should be drawn in
a specific way relative to each other.

Despite all of the progress made, the program has limitations. The problem of
single bonds, as described above, remains unsolved. In particular, solving this
problem requires finding a way to represent the relative orientations of two
atoms in a single bond, and there isn't a single obvious way to do that. The
lack of support for cyclic compounds is another large hole in the feature-set.
Solving this problem would probably involve adding the ability for geometries
to deviate slightly from their VSEPR ideals, in addition to solving the single
bond problem. Other miscellaneous limitations include incomplete SMILES support
and a somewhat primitive UI.
