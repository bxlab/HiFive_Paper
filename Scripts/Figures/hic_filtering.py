#!/usr/bin/env python

import sys

from pyx import *


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular(labelattrs=[text.size(-3)], titleattrs=[text.size(-2)])
painter2 = graph.axis.painter.regular(labelattrs=None, titleattrs=[text.size(-2)])


def main():
    out_fname = sys.argv[1]
    c = canvas.canvas()
    c.insert(sub_figA(), [trafo.translate(0, 0)])
    c.writePDFfile(out_fname)
    return None


def sub_figA():
    b1 = box(['Raw HiC read pairs'])
    b2 = box(['Reference genome'])
    rb0 = rbox(['Map read ends indendently',
                'to reference genome'])
    b3 = box(['Mapped read pairs'])
    b4 = box(['Fend Data Object'])
    rb1 = rbox(['Assign reads to fends based on',
                'coordinates and orientation'])
    d1 = diamond(['Does fend pair contain read pair',
                  'with duplicate coordinates?'])
    d2 = diamond([r'$\mid$Downstream RE site - read coord1$\mid$',
                  r'+ $\mid$Downstream RE site - read coord2$\mid$',
                  r"$\leq$ insert size"])
    rb2 = rbox(['Count fend pair occurences'])
    d3 = diamond(['Do fends in pair originate',
                  'from same RE fragment?'])
    d4 = diamond(['Do fends in pair originate',
                  'from adjacent fragments?'])
    d5 = diamond(['Are fends in the',
                  'same orientation?'])
    b5 = box(['HiC Data Object'])
    rb3 = rbox(['Sum interactions between valid fends',
               r'$\leq$ minimum distance apart for each fend'])
    d6 = diamond([r'Each fend interaction',
                  r'count $\geq$ minimum count$?$'])
    rb4 = rbox(['Remove fends with interaction',
               r'counts $<$ minimum count'])
    b6 = box(['Final filtered set of fends and interactions'])
    a0 = arrow(height=1.0)
    a1 = arrow()
    a2 = arrow(in_text='no')
    a3 = arrow(in_text='yes')
    a4 = arrow(height=1.0, in_text='no')

    c = canvas.canvas()
    c.insert(b1[0])
    pos = -b1[2] / 2.0
    c.insert(a0[0], [trafo.translate(0, pos)])
    c.insert(b2[0], [trafo.translate(-b2[1] / 2.0 - 1.0, pos - a0[1] / 2.0)])
    c.stroke(path.line(-1, pos - a0[1] / 2.0, 0, pos - a0[1] / 2.0), [style.linewidth.Thick])
    pos -= a0[1]
    c.insert(rb0[0], [trafo.translate(0, pos - rb0[2] / 2.0)])
    pos -= rb0[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(b3[0], [trafo.translate(0, pos - b3[2] / 2.0)])
    pos -= b3[2]
    c.insert(a0[0], [trafo.translate(0, pos)])
    c.insert(b4[0], [trafo.translate(-b4[1] / 2.0 - 1.0, pos - a0[1] / 2.0)])
    c.stroke(path.line(-1, pos - a0[1] / 2.0, 0, pos - a0[1] / 2.0), [style.linewidth.Thick])
    pos -= a0[1]
    c.insert(rb1[0], [trafo.translate(0, pos - rb1[2] / 2.0)])
    pos -= rb1[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(d1[0], [trafo.translate(0, pos - d1[2] / 2.0)])
    pos -= d1[2]
    c.insert(a2[0], [trafo.translate(0, pos)])
    pos -= a2[1]
    c.insert(d2[0], [trafo.translate(0, pos - d2[2] / 2.0)])
    pos -= d2[2]
    c.insert(a3[0], [trafo.translate(0, pos)])
    pos -= a3[1]
    c.insert(rb2[0], [trafo.translate(0, pos - rb2[2] / 2.0)])
    pos -= rb2[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(d3[0], [trafo.translate(0, pos - d3[2] / 2.0)])
    pos -= d3[2]
    c.insert(a2[0], [trafo.translate(0, pos)])
    pos -= a2[1]
    c.insert(d4[0], [trafo.translate(0, pos - d4[2] / 2.0)])
    pos -= d4[2]
    temp_pos = pos + d4[2] / 4.0 - a4[1] / 2.0 - b5[2] / 4.0
    c.insert(d5[0], [trafo.translate(-d5[1] / 2.0 - 1.2, temp_pos)])
    c.insert(a4[0], [trafo.translate(0, pos)])
    c.stroke(path.path(path.moveto(-d4[1] / 2.0, pos + d4[2] / 2.0),
        path.curveto(-d4[1] / 2.0 - 0.4, pos + d4[2] / 2.0,
                     -d5[1] / 2.0 - 1.2, temp_pos + d5[2] / 2.0 + 0.5,   
                     -d5[1] / 2.0 - 1.2, temp_pos + d5[2] / 2.0)),   
                     [style.linewidth.Thick, deco.earrow.normal])
    c.text(-d4[1] / 2.0 - 0.1, pos + d4[2] / 2.0 + 0.1, 'yes',
        [text.halign.right, text.valign.bottom, text.size(-3)])
    c.stroke(path.path(path.moveto(-d5[1] / 2.0 - 1.2, temp_pos - d5[2] / 2.0),
        path.curveto(-d5[1] / 2.0 - 1.2, temp_pos - d5[2] / 2.0 - 0.6,
                     -b5[1] / 2.0 - 0.4, pos - a4[1] - b5[2] / 2.0,
                     -b5[1] / 2.0, pos - a4[1] - b5[2] / 2.0)),
                     [style.linewidth.Thick, deco.earrow.normal])
    c.text(-b5[1] / 2.0 - 1, pos - a4[1] - b5[2] / 2.0 - 0.15, 'no',
        [text.halign.right, text.valign.top, text.size(-3)])
    pos -= a4[1]
    c.insert(b5[0], [trafo.translate(0, pos - b5[2] / 2.0)])
    pos -= b5[2]
    c.stroke(path.line(-6, pos - a1[1] / 2.0, 6, pos - a1[1] / 2.0), [style.linewidth.Thick, style.linestyle.dotted])
    c.text(2.5, pos - a1[1] / 2.0 + 0.2, "Filtering at Data Object creation",
           [text.halign.left, text.valign.bottom, text.size(-3)])
    c.text(2.5, pos - a1[1] / 2.0 - 0.2, "Filtering within HiC Object",
           [text.halign.left, text.valign.top, text.size(-3)])
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(rb3[0], [trafo.translate(0, pos - rb3[2] / 2.0)])
    pos -= rb3[2]
    c.insert(a0[0], [trafo.translate(0, pos)])
    pos -= a0[1]
    temp_pos = pos - d6[2] / 4.0 + a0[1] / 2.0 + rb3[2] / 4.0
    c.insert(rb4[0], [trafo.translate(-rb4[1] / 2.0 - 1.5, temp_pos)])
    c.stroke(path.path(path.moveto(-d6[1] / 2.0, pos - d6[2] / 2.0),
        path.curveto(-d6[1] / 2.0 - 0.4, pos - d6[2] / 2.0,
                     -rb4[1] / 2.0 - 1.5, temp_pos - rb4[2] / 2.0 - 0.4,
                     -rb4[1] / 2.0 - 1.5, temp_pos - rb4[2] / 2.0)),
                    [deco.earrow.normal, style.linewidth.Thick])
    c.text(-rb4[1] / 2.0 - 1.6, temp_pos - rb4[2] / 2.0 - 0.3, 'no',
        [text.halign.right, text.valign.top, text.size(-3)])
    c.stroke(path.path(path.moveto(-rb4[1] / 2.0 - 1.5, temp_pos + rb4[2] / 2.0),
        path.curveto(-rb4[1] / 2.0 - 1.5, temp_pos + rb4[2] / 2.0 + 0.4,
                     -rb3[1] / 2.0 - 0.4, pos + a0[1] + rb3[2] / 2.0,
                     -rb3[1] / 2.0, pos + a0[1] + rb3[2] / 2.0)),
                    [deco.earrow.normal, style.linewidth.Thick])
    c.insert(d6[0], [trafo.translate(0, pos - d6[2] / 2.0)])
    pos -= d6[2]
    c.insert(a3[0], [trafo.translate(0, pos)])
    pos -= a3[1]
    c.insert(b6[0], [trafo.translate(0, pos - b6[2] / 2.0)])
    return c


def box(in_text):
    c = canvas.canvas()
    width = 0
    height = len(in_text) * 0.3
    top = height / 2.0 - 0.15
    for i in range(len(in_text)):
        temp = c.text(0, top - i * 0.3, in_text[i], [text.halign.center, text.valign.middle, text.size(-3)])
        width = max(width, unit.tocm(temp.width))
    width += 0.2
    height += 0.2
    c.stroke(path.rect(-width / 2.0, -height / 2.0, width, height), [style.linewidth.Thick])
    return [c, width, height]


def rbox(in_text):
    c = canvas.canvas()
    width = 0
    height = len(in_text) * 0.3
    top = height / 2.0 - 0.15
    for i in range(len(in_text)):
        temp = c.text(0, top - i * 0.3, in_text[i], [text.halign.center, text.valign.middle, text.size(-3)])
        width = max(width, unit.tocm(temp.width))
    width += 0.2
    height += 0.2
    c.stroke(path.rect(-width / 2.0, -height / 2.0, width, height), [style.linewidth.Thick, deformer.smoothed(0.5)])
    return [c, width, height]


def diamond(in_text):
    c = canvas.canvas()
    width = 0
    height = len(in_text) * 0.3
    top = height / 2.0 - 0.15
    for i in range(len(in_text)):
        temp = c.text(0, top - i * 0.3, in_text[i], [text.halign.center, text.valign.middle, text.size(-3)])
        width = max(width, unit.tocm(temp.width))
    width += height * 4.0
    height = width / 4.0
    c.stroke(path.path(path.moveto(-width / 2.0, 0), path.lineto(0, height / 2.0), path.lineto(width / 2.0, 0),
        path.lineto(0, -height / 2.0), path.closepath()), [style.linewidth.Thick])
    return [c, width, height]


def arrow(height=0.5, in_text=None):
    c = canvas.canvas()
    c.stroke(path.line(0, 0, 0, -height), [style.linewidth.Thick, deco.earrow.normal])
    if not in_text is None:
        c.text(0.2, -height / 2.0, in_text, [text.halign.left, text.valign.middle, text.size(-3)])
    return [c, height]


if __name__ == "__main__":
    main()
