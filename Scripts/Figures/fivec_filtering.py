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
    b1 = box(['Raw 5C read pairs'])
    b2 = box(['Fragment probe sequences'])
    rb1 = rbox(['Map read ends independently',
                'to probe sequences'])
    b3 = box(['Mapped read pairs'])
    rb2 = rbox(['Count fragment pair occurences'])
    d1 = diamond(['Are fragment probes in',
                  'the same orientation?'])
    b4 = box(['5C data object'])
    rb3 = rbox(['Sum interactions between valid fragments',
                'in the same defined genome region and',
               r'$\leq$ minimum distance apart for each fragment'])
    d2 = diamond([r'Each fragment interaction',
                  r'count $\geq$ minimum count$?$'])
    rb4 = rbox(['Remove fragments with interaction',
               r'counts $<$ minimum count'])
    b5 = box(['Final filtered set of fragments and interactions'])

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
    c.insert(rb1[0], [trafo.translate(0, pos - rb1[2] / 2.0)])
    pos -= rb1[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(b3[0], [trafo.translate(0, pos - b3[2] / 2.0)])
    pos -= b3[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(rb2[0], [trafo.translate(0, pos - rb2[2] / 2.0)])
    pos -= rb2[2]
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(d1[0], [trafo.translate(0, pos - d1[2] / 2.0)])
    pos -= d1[2]
    c.insert(a2[0], [trafo.translate(0, pos)])
    pos -= a2[1]
    c.insert(b4[0], [trafo.translate(0, pos - b4[2] / 2.0)])
    pos -= b4[2]
    c.stroke(path.line(-6, pos - a1[1] / 2.0, 6, pos - a1[1] / 2.0), [style.linewidth.Thick, style.linestyle.dotted])
    c.text(2.6, pos - a1[1] / 2.0 + 0.2, "Filtering at Data Object creation",
           [text.halign.left, text.valign.bottom, text.size(-3)])
    c.text(2.6, pos - a1[1] / 2.0 - 0.2, "Filtering within 5C Object",
           [text.halign.left, text.valign.top, text.size(-3)])
    c.insert(a1[0], [trafo.translate(0, pos)])
    pos -= a1[1]
    c.insert(rb3[0], [trafo.translate(0, pos - rb3[2] / 2.0)])
    pos -= rb3[2]
    c.insert(a0[0], [trafo.translate(0, pos)])
    pos -= a0[1]
    temp_pos = pos - d2[2] / 4.0 + a0[1] / 2.0 + rb3[2] / 4.0
    c.insert(rb4[0], [trafo.translate(-rb4[1] / 2.0 - 1.2, temp_pos)])
    c.stroke(path.path(path.moveto(-d2[1] / 2.0, pos - d2[2] / 2.0),
        path.curveto(-d2[1] / 2.0 - 0.4, pos - d2[2] / 2.0,
                     -rb4[1] / 2.0 - 1.2, temp_pos - rb4[2] / 2.0 - 0.4,
                     -rb4[1] / 2.0 - 1.2, temp_pos - rb4[2] / 2.0)),
                    [deco.earrow.normal, style.linewidth.Thick])
    c.text(-rb4[1] / 2.0 - 1.2, temp_pos - rb4[2] / 2.0 - 0.4, 'no',
        [text.halign.right, text.valign.top, text.size(-3)])
    c.stroke(path.path(path.moveto(-rb4[1] / 2.0 - 1.2, temp_pos + rb4[2] / 2.0),
        path.curveto(-rb4[1] / 2.0 - 1.2, temp_pos + rb4[2] / 2.0 + 0.4,
                     -rb3[1] / 2.0 - 0.4, pos + a0[1] + rb3[2] / 2.0,
                     -rb3[1] / 2.0, pos + a0[1] + rb3[2] / 2.0)),
                    [deco.earrow.normal, style.linewidth.Thick])
    c.insert(d2[0], [trafo.translate(0, pos - d2[2] / 2.0)])
    pos -= d2[2]
    c.insert(a3[0], [trafo.translate(0, pos)])
    pos -= a3[1]
    c.insert(b5[0], [trafo.translate(0, pos - b5[2] / 2.0)])
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
