#!/usr/bin/python3

from itertools import product
import pyparsing as pp
import argparse
import random
import math
import sys

import OpenGL.GL as gl
import OpenGL.GLU as glu
import OpenGL.GLUT as glut

from linalg import Vector, rotation_vector
import chemistry
import smiles

class MoleculeVisualizer:
    def __init__(self, molecule, hide_hydrogens=False, randomize_single_bonds=False,
                       bond_width=.1, bond_length=2):
        self.molecule = molecule
        self.hide_hydrogens = hide_hydrogens
        self.randomize_single_bonds = randomize_single_bonds
        self.rotate_atoms(self.molecule)

        self.bond_width = bond_width
        self.bond_length = bond_length

        self.quad = glu.gluNewQuadric()
        glu.gluQuadricNormals(self.quad, glu.GLU_SMOOTH)
        glu.gluQuadricDrawStyle(self.quad, glu.GLU_FILL)

        self.cylinder_resolution = 100
        self.sphere_resolution = 100

        self.zoom = 5
        self.rotations = []

        self.scroll_up_button = 3
        self.scroll_down_button = 4
        self.scroll_buttons = [self.scroll_up_button,
                               self.scroll_down_button]
        self.zoom_increment = 0.2

        self.default_element_color = (221/255, 119/255, 1)
        self.bond_color = (128/255,) * 3

        self.drag_start = None

    def _draw_cylinders(self, length, diffs):
        for x_translate, y_translate in diffs:
            gl.glPushMatrix()
            gl.glTranslate(x_translate, y_translate, 0)
            glu.gluCylinder(self.quad, self.bond_width, self.bond_width,
                            length, self.cylinder_resolution, 1)
            gl.glPopMatrix()

    def draw_bond(self, start, end, bond_type, twist_to=None):
        double_diff = self.bond_width + self.bond_width/2
        triple_diff = 2 * self.bond_width + self.bond_width/2
        quadruple_diffs = [double_diff, -double_diff]

        diffs = {
            'single': [(0, 0)],
            'double': [(double_diff, 0), (-double_diff, 0)],
            'triple': [(triple_diff, 0), (0, 0), (-triple_diff, 0)],
            'quadruple': product(quadruple_diffs, quadruple_diffs)
        }[bond_type]

        cylinder_default = Vector(0, 0, 1)
        vec = end - start

        angle = cylinder_default.angle_to(vec)
        rv = rotation_vector(cylinder_default, vec)

        gl.glPushMatrix()
        gl.glTranslate(*start)

        if twist_to is not None:
            rotated_up = Vector(0, 1, 0).rotate(angle, rv)
            twist_angle = rotated_up.right_hand_rule_angle_to(twist_to, vec)

            gl.glRotate(twist_angle, *vec)

        gl.glRotate(angle, *rv)
        gl.glColor(*self.bond_color)
        self._draw_cylinders(vec.length(), diffs)
        gl.glPopMatrix()

    def draw_element(self, element, position):
        gl.glPushMatrix()
        gl.glTranslate(*position)

        if element.color:
            gl.glColor(*element.color)
        else:
            gl.glColor(*self.default_element_color)

        glu.gluSphere(self.quad, 4 * self.bond_width,
                      self.sphere_resolution,
                      self.sphere_resolution)
        gl.glPopMatrix()

    def keyboard(self, key, x, y):
        glut.glutDestroyWindow(self.win)
        sys.exit(0)

    def reshape(self, w, h):
        gl.glViewport(0, 0, w, h)

        gl.glMatrixMode(gl.GL_PROJECTION)
        gl.glLoadIdentity()
        glu.gluPerspective(45, w/h, 1, 100)
        gl.glMatrixMode(gl.GL_MODELVIEW)

    def rotate_atoms(self, atom, previous=None, seen=None):
        if seen is None:
            seen = set()

        seen.add(atom)

        if previous:
            rotate_from = atom.geometry[atom.find_bond_to(previous)]
            rotate_to = -previous.geometry[previous.find_bond_to(atom)]

            angle = rotate_from.angle_to(rotate_to)
            rv = rotation_vector(rotate_from, rotate_to)

            atom.rotate(angle, rv)

            bond_with_previous = atom.bonds[atom.find_bond_to(previous)].order()
            if bond_with_previous == 'single':
                if self.randomize_single_bonds:
                    atom.rotate(random.randrange(0, 360), rotate_to)
                else:
                    cur_plane_normal = atom.plane_normal()
                    prev_plane_normal = previous.plane_normal()

                    if cur_plane_normal is not None and prev_plane_normal is not None:
                        plane_match_angle = cur_plane_normal.right_hand_rule_angle_to(-prev_plane_normal, rotate_to)
                        atom.rotate(plane_match_angle, rotate_to)
            else:
                cur_plane_normal = atom.plane_normal()
                prev_plane_normal = previous.plane_normal()

                if cur_plane_normal is not None and prev_plane_normal is not None:
                    plane_match_angle = cur_plane_normal.right_hand_rule_angle_to(prev_plane_normal, rotate_to)
                    atom.rotate(plane_match_angle, rotate_to)

        for bond in atom.bonds:
            if bond.atom in seen:
                continue

            self.rotate_atoms(bond.atom, atom, seen)

    def draw_atom(self, atom, omit_bonds=None):
        self.draw_element(atom.element, (0, 0, 0))

        for i, (vec, bond) in enumerate(zip(atom.geometry, atom.bonds)):
            if omit_bonds and i in omit_bonds:
                continue

            self.draw_bond(Vector(0, 0, 0), vec * self.bond_length,
                           bond.order(),
                           atom.plane_normal())

    def draw_atoms(self, atom, cur_vector=Vector(0, 0, 0), seen=None):
        if seen is None:
            seen = set()

        seen.add(atom)

        gl.glPushMatrix()
        gl.glTranslate(*cur_vector)

        hydrogen = chemistry.get_element(1)
        is_carbon = atom.element == chemistry.get_element(6)
        omit_bonds = [i for i, bond in enumerate(atom.bonds)
                        if (bond.atom in seen) or
                           (bond.atom.element == hydrogen and
                            is_carbon and self.hide_hydrogens)]
        self.draw_atom(atom, omit_bonds)
        gl.glPopMatrix()

        for i, (bond, vec) in enumerate(zip(atom.bonds, atom.geometry)):
            if i in omit_bonds or bond.atom in seen:
                continue

            next_vector = cur_vector + vec * self.bond_length
            self.draw_atoms(bond.atom, next_vector, seen)

    def display(self):
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)
        gl.glLoadIdentity()

        # do zoom
        glu.gluLookAt(0.0, 0.0, self.zoom, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)

        # do rotation
        for x_rot, y_rot in reversed(self.rotations):
            gl.glRotatef(math.sqrt(x_rot**2 + y_rot**2), y_rot, x_rot, 0)
        #gl.glRotate(x_rot, 0, 1, 0)
        #gl.glRotate(y_rot, 1, 0, 0)

        """
        one, two, three, four = (1, 1, 0), (1, -1, 0), (-1, -1, 0), (-1, 1, 0)
        draw_element(chemistry.elements[0], one)
        draw_single_bond(one, two)
        draw_element(chemistry.elements[0], two)
        draw_double_bond(two, three)
        draw_element(chemistry.elements[0], three)
        draw_triple_bond(three, four)
        draw_element(chemistry.elements[0], four)
        draw_quadruple_bond(four, one)
        """
        #self.draw_atoms(self.molecule)
        gl.glCallList(self.displist)

        #gl.glFlush()
        glut.glutSwapBuffers()

    def wheel(self, direction):
        if direction == self.scroll_up_button:
            if self.zoom > 2:
                self.zoom -= self.zoom_increment
        elif direction == self.scroll_down_button:
            self.zoom += self.zoom_increment

        glut.glutPostRedisplay()

    def mouse(self, button, state, x, y):
        if button in self.scroll_buttons:
            if state == glut.GLUT_UP:
                self.wheel(button)

            return

        if button == glut.GLUT_LEFT_BUTTON:
            if state == glut.GLUT_DOWN and not self.drag_start:
                self.drag_start = (x, y)
                self.rotations.append((0, 0))
            elif state == glut.GLUT_UP and self.drag_start:
                self.drag_start = None

    def rotate(self, cur_x, cur_y):
        if self.drag_start:
            drag_x, drag_y = self.drag_start
            diff_x, diff_y = cur_x - drag_x, cur_y - drag_y
            #before_rot_x, before_rot_y = rotation_before_drag
            self.rotations[-1] = diff_x, diff_y

            glut.glutPostRedisplay()

    def run(self):
        glut.glutInit(sys.argv)
        glut.glutInitWindowSize(500, 500)
        glut.glutInitDisplayMode(glut.GLUT_RGBA | glut.GLUT_DOUBLE | glut.GLUT_DEPTH)
        self.win = glut.glutCreateWindow(b'Molecular visualizer')

        glut.glutDisplayFunc(self.display)
        glut.glutReshapeFunc(self.reshape)
        glut.glutKeyboardFunc(self.keyboard)
        glut.glutMotionFunc(self.rotate)
        glut.glutMouseFunc(self.mouse)

        gl.glClearColor(0, 0, 0, 1)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glEnable(gl.GL_LIGHTING)
        gl.glEnable(gl.GL_LIGHT0)

        gl.glColorMaterial(gl.GL_FRONT, gl.GL_DIFFUSE)
        gl.glEnable(gl.GL_COLOR_MATERIAL)

        # very diffuse and dark specular highlights, to make black visible
        gl.glMaterial(gl.GL_FRONT, gl.GL_SPECULAR, (.1, .1, .1, 1))
        gl.glMaterial(gl.GL_FRONT, gl.GL_SHININESS, 5)

        self.displist = gl.glGenLists(1)
        gl.glNewList(self.displist, gl.GL_COMPILE)
        self.draw_atoms(self.molecule)
        gl.glEndList()

        glut.glutMainLoop()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
                 description='Visualizes molecules specified by a SMILES string.'
             )
    parser.add_argument('--hide-hydrogens', action='store_true', default=False,
                        help='Whether to hide hydrogen atoms bonded to carbon '
                             'atoms in the visualization.')
    parser.add_argument('--randomize-single-bonds', action='store_true', default=False,
                        help='Whether to randomize single bond twists in order to '
                             'give some impression of the behavior of single bonds.')
    parser.add_argument('smiles', help='The SMILES string to visualize.')
    args = parser.parse_args()

    try:
        parsed = smiles.smiles.parseString(args.smiles)
    except pp.ParseException:
        print('Unable to parse the given SMILES string. This could be '
              'because it is invalid or because you used a feature that '
              'the parser cannot currently handle. Such features currently include '
              'aromaticity and non-tetrahedral stereochemistry specifications.')
        sys.exit(1)

    if len(parsed) > 1:
        print('Multiple compounds are not currently supported. Sorry! '
              'Only visualizing the first compound.')

    parsed = parsed[0]
    cyclic_linkages = {}
    graph = smiles.construct_atom_graph(parsed, cyclic_linkages)

    if len(cyclic_linkages) > 0:
        print('Warning: cyclic compounds are not current supported. '
              'This will most likely not render correctly.')

        for num, atoms in cyclic_linkages.items():
            if len(atoms) != 2:
                print('Invalid cyclic compound specification: exactly two '
                      'atoms must have a certain number!')

            left, right = atoms
            left.replace_bond(num, right)
            right.replace_bond(num, left)

    MoleculeVisualizer(graph, args.hide_hydrogens, args.randomize_single_bonds).run()
