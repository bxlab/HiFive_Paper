#!/usr/bin/env python
#(c) 2013 Emory University. All Rights Reserved
# Code written by: Michael Sauria (mgehrin@emory.edu)

import sys
import os

from pyx import *
from PIL import Image
import numpy
from math import log, sin, cos


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labelattrs=[text.size(-3)], titleattrs=[text.size(-2)] )
pi = numpy.pi

def main():
    out_fname = sys.argv[1]
    c = canvas.canvas()
    c.insert(plot_subfigA())
    c.writePDFfile(out_fname)


def plot_subfigA():
    c = canvas.canvas()
    r = (2.54 * 6) / 20.0
    w = r * 0.5
    x = 0.0
    y = 0.1 + r
    z = 1.5 * w
    a0 = 220.0
    a1 = 70.0
    x1 = 2 * r * cos(a0 / 180.0 * pi)
    y1 = y + 2 * r * sin(a0 / 180.0 * pi)
    x2 = x - r * cos(a0 / 180.0 * pi)
    y2 = y + r * sin(a0 / 180.0 * pi)
    s = (y2 - y) / (x2 - x)
    c.stroke(path.line(r * 2, -r * 2.25, r * 2, -r * 11.4), [style.linestyle.dashed, color.gray(0.7)])
    c.stroke(path.line(r * 6, -r * 2.25, r * 6, -r * 11.4), [style.linestyle.dashed, color.gray(0.7)])
    c.stroke(path.line(r * 10, -r * 2.25, r * 10, -r * 11.4), [style.linestyle.dashed, color.gray(0.7)])
    c.stroke(path.line(r * 14, -r * 2.25, r * 14, -r * 11.4), [style.linestyle.dashed, color.gray(0.7)])
    c.stroke(path.rect(-2 * r, -11.65 * r, 20 * r, 12.65 * r), [color.gray(1.0)])
    ## Untreated cells
    cx = canvas.canvas()
    cx.stroke(path.path(path.arc(x, y, r, a0, 270 * 2 - a0)), [color.rgb.red, style.linewidth.THick])
    cx.stroke(path.path(path.arc(x1, y1, r, a0 - 180, a1)), [color.gray(0.25), style.linewidth.THick])
    cx.stroke(path.path(path.arcn(-x1, y1, r, 360 - a0, 180 - a1)), [color.gray(0.25), style.linewidth.THick])
    cx.stroke(path.path(path.arc(x, -y, r, a0 - 180, 360 - a0)), [color.rgb.blue, style.linewidth.THick])
    cx.stroke(path.path(path.arcn(x1, -y1, r, 270 * 2 - a0, 360 - a1)), [color.gray(0.25), style.linewidth.THick])
    cx.stroke(path.path(path.arc(-x1, -y1, r, a0, 180 + a1)), [color.gray(0.25), style.linewidth.THick])
    cx.stroke(path.circle(0, 0, (y - r) * 1.2),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    cx.stroke(path.line(x2 - 0.1 * s, y2 - 0.1, x2 + 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    cx.stroke(path.line(x2 - 0.1 * s, -y2 + 0.1, x2 + 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    cx.stroke(path.line(-x2 + 0.1 * s, y2 - 0.1, -x2 - 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    cx.stroke(path.line(-x2 + 0.1 * s, -y2 + 0.1, -x2 - 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c.insert(cx, [trafo.translate(0, 0)])
    c.stroke(path.line(r * 2 - w, 0, r * 2 + w, 0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 2, w * 0.5, "Cross-link", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 2, -w * 0.5, "Cells", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Cross-linked cells
    c0 = canvas.canvas()
    c0.stroke(path.path(path.arc(x, y, r, a0, 270 * 2 - a0)), [color.rgb.red, style.linewidth.THick])
    c0.stroke(path.path(path.arc(x1, y1, r, a0 - 180, a1)), [color.gray(0.25), style.linewidth.THick])
    c0.stroke(path.path(path.arcn(-x1, y1, r, 360 - a0, 180 - a1)), [color.gray(0.25), style.linewidth.THick])
    c0.stroke(path.path(path.arc(x, -y, r, a0 - 180, 360 - a0)), [color.rgb.blue, style.linewidth.THick])
    c0.stroke(path.path(path.arcn(x1, -y1, r, 270 * 2 - a0, 360 - a1)), [color.gray(0.25), style.linewidth.THick])
    c0.stroke(path.path(path.arc(-x1, -y1, r, a0, 180 + a1)), [color.gray(0.25), style.linewidth.THick])
    c0.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c0.stroke(path.line(x2 - 0.1 * s, y2 - 0.1, x2 + 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c0.stroke(path.line(x2 - 0.1 * s, -y2 + 0.1, x2 + 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c0.stroke(path.line(-x2 + 0.1 * s, y2 - 0.1, -x2 - 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c0.stroke(path.line(-x2 + 0.1 * s, -y2 + 0.1, -x2 - 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c.insert(c0, [trafo.translate(4 * r, 0)])
    c.stroke(path.line(r * 6 - w, 0, r * 6 + w, 0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 6, w * 0.5, "Restriction Enzyme", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 6, -w * 0.5, "Digestion", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Digested cells
    c1 = canvas.canvas()
    c1.stroke(path.path(path.arc(x, y, r, a0, 270 * 2 - a0)), [color.rgb.red, style.linewidth.THick])
    c1.stroke(path.path(path.arc(x, -y, r, a0 - 180, 360 - a0)), [color.rgb.blue, style.linewidth.THick])
    c1.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c1.stroke(path.line(x2 - 0.1 * s, y2 - 0.1, x2 + 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c1.stroke(path.line(x2 - 0.1 * s, -y2 + 0.1, x2 + 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c1.stroke(path.line(-x2 + 0.1 * s, y2 - 0.1, -x2 - 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c1.stroke(path.line(-x2 + 0.1 * s, -y2 + 0.1, -x2 - 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c.insert(c1, [trafo.translate(r * 8.0, 0)])
    c.stroke(path.line(r * 10 - w, 0, r * 10 + w, 0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 10, w * 0.5, "Ligation", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 10, -w * 0.5, "Reaction", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Ligated cells
    ligated = canvas.canvas()
    ligated.stroke(path.circle(0, 0, 0.035), [style.linewidth.THIN, deco.filled([color.cmyk.Yellow])])
    c2 = canvas.canvas()
    a2 = 295
    x3 = x + r * cos(a2 / 180.0 * pi)
    y3 = y + r * sin(a2 / 180.0 * pi)
    s3 = (y3 - y) / (x3 - x)
    x4 = x3 - y3 / s3
    r3 = (y3 ** 2.0 + (x4 - x3) ** 2.0) ** 0.5
    r4 = (0.1 ** 2.0 + (0.1 * s) ** 2.0) ** 0.5
    c2.stroke(path.path(path.arc(x, y, r, a0, a2), path.arcn(x4, 0, r3, a2 - 180, 0)),
              [color.rgb.red, style.linewidth.THick])
    c2.stroke(path.path(path.arcn(x4, 0, r3, 0, 560 - a2), path.arc(x, -y, r, 360 - a2, 360 - a0)),
              [color.rgb.blue, style.linewidth.THick])
    c2.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c2.stroke(path.line(x4 + r3 - r4, 0, x4 + r3 + r4, 0), [style.linewidth.Thick])
    c2.stroke(path.line(-x2 + 0.1 * s, y2 - 0.1, -x2 - 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c2.stroke(path.line(-x2 + 0.1 * s, -y2 + 0.1, -x2 - 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c2.insert(ligated, [trafo.translate(x4 + r3, 0)])
    c.insert(c2, [trafo.translate(r * 12.0, 0)])
    c.stroke(path.line(r * 14 - w, 0, r * 14 + w, 0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 14, w * 0.5, "Reverse", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 14, -w * 0.5, "Cross-linking", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Ligated uncross-linked
    c3 = canvas.canvas()
    c3.stroke(path.line(-r, 0, 0, 0), [color.rgb.red, style.linewidth.THick])
    c3.stroke(path.line(r, 0, 0, 0), [color.rgb.blue, style.linewidth.THick])
    c3.stroke(path.line(0, -r4, 0, r4), [style.linewidth.Thick])
    c3.stroke(path.line(-r, -r4, -r, r4), [style.linewidth.Thick])
    c3.stroke(path.line(r, -r4, r, r4), [style.linewidth.Thick])
    c3.stroke(path.circle(0, r * 0.45, (y - r) * 1.2), [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c3.fill(path.path(path.moveto(-r * 0.05, w * 0.125), path.lineto(-r * 0.3, w * 0.125),
            path.lineto(-r * 0.3, w * 0.25), path.lineto(-r * 0.175, w * 0.25),
            path.lineto(-r * 0.175, w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.RedOrange])
    c3.fill(path.path(path.moveto(r * 0.05, -w * 0.125), path.lineto(r * 0.3, -w * 0.125),
            path.lineto(r * 0.3, -w * 0.25), path.lineto(r * 0.175, -w * 0.25),
            path.lineto(r * 0.175, -w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.CornflowerBlue])
    c3.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c3, [trafo.translate(r * 16.0, 0)])
    ## 3C
    c.stroke(path.line(r * 16, -r * 0.5, r * 16, -r * 2), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 16, -r * 2.25, "Quantitative PCR", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Circularization of cross-linked fragments
    c5 = canvas.canvas()
    c5.stroke(path.path(path.arcn(-x4, 0, r3, 180, 360 - a2), path.arc(x, y, r, 540 - a2, a2),
              path.arcn(x4, 0, r3, a2 - 180, 0)), [color.rgb.red, style.linewidth.THick])
    c5.stroke(path.path(path.arcn(x4, 0, r3, 0, 560 - a2), path.arc(x, -y, r, 360 - a2, a2 - 180),
              path.arcn(-x4, 0, r3, a2, 180)), [color.rgb.blue, style.linewidth.THick])
    c5.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c5.stroke(path.line(x4 + r3 - r4, 0, x4 + r3 + r4, 0), [style.linewidth.Thick])
    c5.stroke(path.line(-x4 - r3 + r4, 0, -x4 - r3 - r4, 0), [style.linewidth.Thick])
    c5.insert(ligated, [trafo.translate(x4 + r3, 0)])
    c5.insert(ligated, [trafo.translate(-x4 - r3, 0)])
    c.insert(c5, [trafo.translate(r * 8.0, -r * 2.75)])
    c.stroke(path.line(r * 8, -r * 0.5, r * 8, -r * 2), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 8.2, -r * 1.75 + w * 0.25, "Ligation", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 8.2, -r * 1.75 - w * 0.25, "Reaction", [text.halign.left, text.valign.middle, text.size(-4)])
    ## 4C uncross-linked
    c6 = canvas.canvas()
    c6.stroke(path.path(path.arcn(0, 0, x4, 180, 0)), [color.rgb.red, style.linewidth.THick])
    c6.stroke(path.path(path.arc(0, 0, x4, 180, 0)), [color.rgb.blue, style.linewidth.THick])
    c6.stroke(path.line(-x4 - r4, 0, -x4 + r4, 0), [style.linewidth.Thick])
    c6.stroke(path.line(x4 - r4, 0, x4 + r4, 0), [style.linewidth.Thick])
    c6.stroke(path.circle(0, 0, (y - r) * 1.2), [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c6.fill(path.path(path.arcn(0, 0, x4 + w * 0.25, 163, 150),
            path.arc(0, 0, x4 + w * 0.125, 150, 175),
            path.lineto((x4 + w * 0.35) * cos(163 / 180. * pi), (x4 + w * 0.35) * sin(163 / 180. * pi)),
            path.closepath()), [style.linewidth.THIN, color.cmyk.RedOrange])
    c6.fill(path.path(path.arc(0, 0, x4 - w * 0.25, 17, 40),
            path.arcn(0, 0, x4 - w * 0.125, 40, 5),
            path.lineto((x4 - w * 0.35) * cos(17 / 180. * pi), (x4 - w * 0.35) * sin(17 / 180. * pi)),
            path.closepath()), [style.linewidth.THIN, color.cmyk.RedOrange])
    c6.insert(ligated, [trafo.translate(x4, 0)])
    c6.insert(ligated, [trafo.translate(-x4, 0)])
    c.insert(c6, [trafo.translate(r * 8.0, -r * 4.75)])
    c.stroke(path.line(r * 8, -r * 3.5, r * 8, -r * 4.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 8.2, -r * 3.75 + w * 0.25, "Reverse", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 8.2, -r * 3.75 - w * 0.25, "Cross-linking", [text.halign.left, text.valign.middle, text.size(-4)])
    ## 4C
    c.stroke(path.line(r * 8, -r * 5.5, r * 8, -r * 6), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 8, -r * 6.25, "Cloning and Sequencing", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 8, -r * 6.25 - w * 0.5, "or", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 8, -r * 6.25 - w * 1.0, "PCR and Microarray", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Biotin ligation
    c7 = canvas.canvas()
    c7.stroke(path.path(path.arc(x, y, r, a0, a2), path.arcn(x4, 0, r3, a2 - 180, 0)),
              [color.rgb.red, style.linewidth.THick])
    c7.stroke(path.path(path.arcn(x4, 0, r3, 0, 560 - a2), path.arc(x, -y, r, 360 - a2, 360 - a0)),
              [color.rgb.blue, style.linewidth.THick])
    c7.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c7.stroke(path.line(x4 + r3 - r4, 0, x4 + r3 + r4 * 1.5, 0), [style.linewidth.THick])
    c7.stroke(path.line(-x2 + 0.1 * s, y2 - 0.1, -x2 - 0.1 * s, y2 + 0.1), [style.linewidth.Thick])
    c7.stroke(path.line(-x2 + 0.1 * s, -y2 + 0.1, -x2 - 0.1 * s, -y2 - 0.1), [style.linewidth.Thick])
    c7.stroke(path.circle(x2 + r4 * 1.5 + w * 0.08, 0, w * 0.16),
              [style.linewidth.thick, deco.filled([color.rgb.green])])
    c7.insert(ligated, [trafo.translate(x4 + r3, 0)])
    c.insert(c7, [trafo.translate(r * 4.0, -r * 2.75)])
    c.stroke(path.path(path.moveto(r * 8, -r * 1), path.lineto(r * 6.33, -r * 1), path.lineto(r * 5.67, -r * 1.25),
             path.lineto(r * 4, -r * 1.25), path.lineto(r * 4, -r * 2)), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 3.8, -r * 1.75 + w * 0.25, "Biotin", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 3.8, -r * 1.75 - w * 0.25, "Labeling", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 1.75 + w * 0.25, "and Reaction", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 1.75 - w * 0.25, "Reaction", [text.halign.left, text.valign.middle, text.size(-4)])
    ## Uncross-linked biotin
    c8 = canvas.canvas()
    c8.stroke(path.line(-r, 0, 0, 0), [color.rgb.red, style.linewidth.THick])
    c8.stroke(path.line(r, 0, 0, 0), [color.rgb.blue, style.linewidth.THick])
    c8.stroke(path.line(-r, r4, -r, -r4), [style.linewidth.Thick])
    c8.stroke(path.line(r, r4, r, -r4), [style.linewidth.Thick])
    c8.stroke(path.line(0, r4, 0, -r4 * 1.5), [style.linewidth.THick])
    c8.stroke(path.circle(0, -r4 * 1.5 - w * 0.08, w * 0.16),
              [style.linewidth.thick, deco.filled([color.rgb.green])])
    c8.stroke(path.circle(0, r * 0.45, (y - r) * 1.2), [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c8.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c8, [trafo.translate(r * 4, -r * 4.875)])
    c.stroke(path.line(r * 4, -r * 3.5, r * 4, -r * 4.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 4.2, -r * 3.75 + w * 0.25, "Reverse", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 3.75 - w * 0.25, "Cross-linking", [text.halign.left, text.valign.middle, text.size(-4)])
    ## HiC sonicate and precipitate
    fragment = canvas.canvas()
    fragment.stroke(path.path(path.moveto(0, -r4), path.lineto(-r4 * 0.25, 0), path.lineto(r4 * 0.25, 0),
                    path.lineto(0, r4)), [style.linewidth.Thick])
    c9 = canvas.canvas()
    c9.stroke(path.line(-r * 0.75, 0, 0, 0), [color.rgb.red, style.linewidth.THick])
    c9.stroke(path.line(r * 0.75, 0, 0, 0), [color.rgb.blue, style.linewidth.THick])
    c9.stroke(path.line(0, r4, 0, -r4 * 1.5), [style.linewidth.THick])
    c9.stroke(path.circle(0, -r4 * 1.5 - w * 0.08, w * 0.16),
              [style.linewidth.thick, deco.filled([color.rgb.green])])
    strep = canvas.canvas()
    strep.stroke(path.line(0, 0, 0, -r4 * 1.5), [style.linewidth.THick, color.cmyk.DarkOrchid])
    strep.stroke(path.path(path.moveto(-r4, r4), path.lineto(0, 0),
              path.lineto(r4, r4)), [style.linewidth.THick, color.cmyk.DarkOrchid])
    c9.insert(strep, [trafo.translate(0, -r4 * 2.0 - w * 0.24)])
    c9.insert(fragment, [trafo.translate(-r * 0.75, 0)])
    c9.insert(fragment, [trafo.translate(r * 0.75, 0)])
    c9.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c9, [trafo.translate(r * 4, -r * 6.6)])
    c.text(r * 3.8, -r * 5.75 + w * 0.25, "Sonication",
           [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 3.8, -r * 5.75 - w * 0.25, "Fragmentation",
           [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 5.75 + w * 0.25, "and Streptavidin",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 5.75 - w * 0.25, "Precipitation",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.stroke(path.line(r * 4, -r * 5.5, r * 4, -r * 6), [style.linewidth.Thick, deco.earrow.normal])
    ## HiC linker addition
    linker = canvas.canvas()
    linker.stroke(path.rect(-r4, -r4 * 0.5, r4 * 2, r4), [deco.filled([color.cmyk.CadetBlue])])
    c10 = canvas.canvas()
    c10.stroke(path.line(-r * 0.75, 0, 0, 0), [color.rgb.red, style.linewidth.THick])
    c10.stroke(path.line(r * 0.75, 0, 0, 0), [color.rgb.blue, style.linewidth.THick])
    c10.stroke(path.line(0, r4, 0, -r4), [style.linewidth.Thick])
    c10.insert(linker, [trafo.translate(-r * 0.75 - r4, 0)])
    c10.insert(linker, [trafo.translate(r * 0.75 + r4, 0)])
    c10.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c10, [trafo.translate(r * 4, -r * 8.5)])
    c.text(r * 4.2, -r * 7.75 + w * 0.5, "Addition of",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 7.75, "Sequencing",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 4.2, -r * 7.75 - w * 0.5, "Linkers",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.stroke(path.line(r * 4, -r * 7.5, r * 4, -r * 8), [style.linewidth.Thick, deco.earrow.normal])
    ## HiC
    c.stroke(path.line(r * 4, -r * 9, r * 4, -r * 9.5), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 4, -r * 9.75, "High-throughput", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 4, -r * 9.75 - w * 0.5, "Sequencing", [text.halign.center, text.valign.middle, text.size(-4)])
    ## ChIA-PET fragmentation and precipitation
    c.stroke(path.path(path.moveto(r * 4, -r * 0.5), path.lineto(r * 4, -r * 1), path.lineto(r * 2.33, -r * 1),
             path.lineto(r * 1.67, -r * 1.25), path.lineto(0, -r * 1.25), path.lineto(0, -r * 2)),
             [style.linewidth.Thick, deco.earrow.normal])
    c.text(-r * 0.2, -r * 1.75 + w * 0.25, "Sonication",
           [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(-r * 0.2, -r * 1.75 - w * 0.25, "Fragmentation",
           [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 0.2, -r * 1.75 + w * 0.25, "and Antibody",
           [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 0.2, -r * 1.75 - w * 0.25, "Precipitation",
           [text.halign.left, text.valign.middle, text.size(-4)])
    ## ChIA-PET sonication and precipitation    antibody = canvas.canvas()
    antibody = canvas.canvas()
    antibody.stroke(path.line(0, 0, 0, -r4 * 1.5), [style.linewidth.THick, color.cmyk.PineGreen])
    antibody.stroke(path.path(path.moveto(-r4, r4), path.lineto(0, 0),
              path.lineto(r4, r4)), [style.linewidth.THick, color.cmyk.PineGreen])
    c11 = canvas.canvas()
    c11.stroke(path.path(path.arc(x, y, r, a0, 270 * 2 - a0)), [color.rgb.red, style.linewidth.THick])
    c11.stroke(path.path(path.arc(x, -y, r, a0 - 180, 360 - a0)), [color.rgb.blue, style.linewidth.THick])
    c11.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c11.insert(fragment, [trafo.rotate(a0 - 180), trafo.translate(x2, y2)])
    c11.insert(fragment, [trafo.rotate(180 - a0), trafo.translate(-x2, y2)])
    c11.insert(fragment, [trafo.rotate(180 - a0), trafo.translate(x2, -y2)])
    c11.insert(fragment, [trafo.rotate(a0 - 180), trafo.translate(-x2, -y2)])
    c11.insert(antibody, [trafo.translate(0, -(y - r) * 2.6)])
    c.insert(c11, [trafo.translate(0, -r * 2.75)])
    c.stroke(path.line(0, -r * 3.5, 0, -r * 4.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(0.2, -r * 3.75 + w * 0.25, "Addition of", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(0.2, -r * 3.75 - w * 0.25, "RE Linkers", [text.halign.left, text.valign.middle, text.size(-4)])
    ## ChIA-PET linkered
    linker2 = canvas.canvas()
    linker2.stroke(path.rect(-r4 * 0.5, -r4 * 0.5, r4, r4), [deco.filled([color.cmyk.OliveGreen])])
    c12 = canvas.canvas()
    c12.stroke(path.path(path.arc(x, y, r, a0, 270 * 2 - a0)), [color.rgb.red, style.linewidth.THick])
    c12.stroke(path.path(path.arc(x, -y, r, a0 - 180, 360 - a0)), [color.rgb.blue, style.linewidth.THick])
    c12.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c12.insert(linker2, [trafo.translate(r4 * 0.5, 0), trafo.rotate(a0 - 180), trafo.translate(x2, y2)])
    c12.insert(linker2, [trafo.translate(-r4 * 0.5, 0), trafo.rotate(180 - a0), trafo.translate(-x2, y2)])
    c12.insert(linker2, [trafo.translate(r4 * 0.5, 0), trafo.rotate(180 - a0), trafo.translate(x2, -y2)])
    c12.insert(linker2, [trafo.translate(-r4 * 0.5, 0), trafo.rotate(a0 - 180), trafo.translate(-x2, -y2)])
    c.insert(c12, [trafo.translate(0, -r * 4.75)])
    c.stroke(path.line(0, -r * 5.5, 0, -r * 6), [style.linewidth.Thick, deco.earrow.normal])
    c.text(0.2, -r * 5.75 + w * 0.25, "Ligation", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(0.2, -r * 5.75 - w * 0.25, "Reaction", [text.halign.left, text.valign.middle, text.size(-4)])
    ## ChIA-PET ligated
    c13 = canvas.canvas()
    c13.stroke(path.path(path.moveto(-r * 0.75, r4), path.curveto(-r * 0.75, r4 * 4, -0.4 * r, y - r, 0, y - r),
               path.curveto(0.4 * r, y - r, r * 0.75, r4 * 4, r * 0.75, r4)),
               [style.linewidth.THick, color.rgb.red])
    c13.stroke(path.path(path.moveto(-r * 0.75, -r4), path.curveto(-r * 0.75, -r4 * 4, -0.4 * r, -y + r, 0, -y + r),
               path.curveto(0.4 * r, -y + r, r * 0.75, -r4 * 4, r * 0.75, -r4)),
               [style.linewidth.THick, color.rgb.blue])
    c13.stroke(path.rect(-0.1 * r, -(y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c13.insert(linker2, [trafo.translate(0.75 * r, r4 * 0.5)])
    c13.insert(linker2, [trafo.translate(0.75 * r, -r4 * 0.5)])
    c13.insert(linker2, [trafo.translate(-0.75 * r, r4 * 0.5)])
    c13.insert(linker2, [trafo.translate(-0.75 * r, -r4 * 0.5)])
    c13.stroke(path.line(0.75 * r - 0.15 * s, r4 * 1.5, 0.75 * r + 0.15 * s, r4 * 1.5), [style.linewidth.Thick])
    c13.stroke(path.line(0.75 * r - 0.15 * s, -r4 * 1.5, 0.75 * r + 0.15 * s, -r4 * 1.5), [style.linewidth.Thick])
    c13.stroke(path.line(-0.75 * r - 0.15 * s, r4 * 1.5, -0.75 * r + 0.15 * s, r4 * 1.5), [style.linewidth.Thick])
    c13.stroke(path.line(-0.75 * r - 0.15 * s, -r4 * 1.5, -0.75 * r + 0.15 * s, -r4 * 1.5), [style.linewidth.Thick])
    c13.insert(ligated, [trafo.translate(0.75 * r, 0)])
    c13.insert(ligated, [trafo.translate(-0.75 * r, 0)])
    c.insert(c13, [trafo.translate(0, -r * 6.75)])
    c.stroke(path.line(0, -r * 7.5, 0, -r * 8), [style.linewidth.Thick, deco.earrow.normal])
    c.text(-0.2, -r * 7.75 + w * 0.25, "Restriction", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(-0.2, -r * 7.75 - w * 0.25, "Enzyme Digest", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(0.2, -r * 7.75 + w * 0.25, "and Sequencing", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(0.2, -r * 7.75 - w * 0.25, "Linker Addition", [text.halign.left, text.valign.middle, text.size(-4)])
    ## ChIA-PET linker addition
    c14 = canvas.canvas()
    c14.stroke(path.line(-r4, 0, -r4 * 2, 0), [color.rgb.red, style.linewidth.THick])
    c14.stroke(path.line(r4, 0, r4 * 2, 0), [color.rgb.blue, style.linewidth.THick])
    c14.insert(linker2, [trafo.translate(-r4 * 0.5, 0)])
    c14.insert(linker2, [trafo.translate(r4 * 0.5, 0)])
    c14.insert(linker, [trafo.translate(-r4 * 3, 0)])
    c14.insert(linker, [trafo.translate(r4 * 3, 0)])
    c14.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c14, [trafo.translate(0, -r * 8.5)])
    ## ChIA-PET
    c.stroke(path.line(0, -r * 9, 0, -r * 9.5), [style.linewidth.Thick, deco.earrow.normal])
    c.text(0, -r * 9.75, "High-throughput", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(0, -r * 9.75 - w * 0.5, "Sequencing", [text.halign.center, text.valign.middle, text.size(-4)])
    ## 5C annealing
    c.stroke(path.line(r * 12, -r * 0.5, r * 12, -r * 2), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 11.8, -r * 1.75 + w * 0.25, "Reverse", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 11.8, -r * 1.75 - w * 0.25, "Cross-linking", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 1.75 + w * 0.25, "and Primer", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 1.75 - w * 0.25, "Annealing", [text.halign.left, text.valign.middle, text.size(-4)])
    c15 = canvas.canvas()
    c15.stroke(path.line(-r, 0, 0, 0), [color.rgb.red, style.linewidth.THick])
    c15.stroke(path.line(r, 0, 0, 0), [color.rgb.blue, style.linewidth.THick])
    c15.stroke(path.line(0, -r4, 0, r4), [style.linewidth.Thick])
    c15.stroke(path.line(-r, -r4, -r, r4), [style.linewidth.Thick])
    c15.stroke(path.line(r, -r4, r, r4), [style.linewidth.Thick])
    c15.stroke(path.circle(0, r * 0.45, (y - r) * 1.2), [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c15.fill(path.path(path.moveto(0, w * 0.125), path.lineto(-r * 0.35, w * 0.125),
            path.lineto(-r * 0.35, w * 0.25), path.lineto(-r * 0.125, w * 0.25),
            path.lineto(-r * 0.125, w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.RedOrange])
    a = 0.15 * (2.0 ** 0.5) * r
    c15.fill(path.path(path.moveto(-r * 0.35, w * 0.125), path.lineto(-r * 0.35 - a, w * 0.125 + a),
             path.lineto(-r * 0.35 - a + w * 0.125, w * 0.25 + a), path.lineto(-r * 0.35 + w * 0.125, w * 0.25),
             path.closepath()), [style.linewidth.THIN, color.cmyk.CadetBlue])
    c15.fill(path.path(path.moveto(0, w * 0.125), path.lineto(r * 0.35, w * 0.125),
            path.lineto(r * 0.35, w * 0.25), path.lineto(r * 0.125, w * 0.25),
            path.lineto(r * 0.125, w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.CornflowerBlue])
    c15.fill(path.path(path.moveto(r * 0.35, w * 0.125), path.lineto(r * 0.35 + a, w * 0.125 + a),
             path.lineto(r * 0.35 + a - w * 0.125, w * 0.25 + a), path.lineto(r * 0.35 - w * 0.125, w * 0.25),
             path.closepath()), [style.linewidth.THIN, color.cmyk.CadetBlue])
    c.insert(c15, [trafo.translate(r * 12.0, -r * 3.0)])
    ## 5C PCR
    c.stroke(path.line(r * 12, -r * 3.5, r * 12, -r * 4.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 11.8, -r * 3.75, "Nick Repair", [text.halign.right, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 3.75 + w * 0.25, "and PCR", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 3.75 - w * 0.25, "Enrichment", [text.halign.left, text.valign.middle, text.size(-4)])
    c16 = canvas.canvas()
    c16.stroke(path.line(0, 0, -0.25, 0), [color.cmyk.RedOrange, style.linewidth.THick])
    c16.stroke(path.line(-0.25, 0, -0.5, 0), [color.cmyk.CadetBlue, style.linewidth.THick])
    c16.stroke(path.line(0, 0, 0.25, 0), [color.cmyk.CornflowerBlue, style.linewidth.THick])
    c16.stroke(path.line(0.25, 0, 0.5, 0), [color.cmyk.CadetBlue, style.linewidth.THick])
    c16.fill(path.path(path.moveto(-0.25, w * 0.125), path.lineto(-0.5, w * 0.125),
            path.lineto(-0.5, w * 0.25), path.lineto(-0.375, w * 0.25),
            path.lineto(-0.375, w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.CadetBlue])
    c16.fill(path.path(path.moveto(0.25, -w * 0.125), path.lineto(0.5, -w * 0.125),
            path.lineto(0.5, -w * 0.25), path.lineto(0.375, -w * 0.25),
            path.lineto(0.375, -w * 0.325), path.closepath()), [style.linewidth.THIN, color.cmyk.CadetBlue])
    c16.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c16, [trafo.translate(r * 12., -r * 4.75)])
    ## 5C Linkered
    c.stroke(path.line(r * 12, -r * 5.5, r * 12, -r * 6.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 12.2, -r * 5.75 + w * 0.5, "Addition of", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 5.75, "Sequencing", [text.halign.left, text.valign.middle, text.size(-4)])
    c.text(r * 12.2, -r * 5.75 - w * 0.5, "Linkers", [text.halign.left, text.valign.middle, text.size(-4)])
    c17 = canvas.canvas()
    c17.stroke(path.line(0, 0, -0.25, 0), [color.cmyk.RedOrange, style.linewidth.THick])
    c17.stroke(path.line(-0.25, 0, -0.5, 0), [color.cmyk.CadetBlue, style.linewidth.THick])
    c17.stroke(path.line(0, 0, 0.25, 0), [color.cmyk.CornflowerBlue, style.linewidth.THick])
    c17.stroke(path.line(0.25, 0, 0.5, 0), [color.cmyk.CadetBlue, style.linewidth.THick])
    c17.insert(linker, [trafo.translate(-0.5 - r4, 0)])
    c17.insert(linker, [trafo.translate(+0.5 + r4, 0)])
    c17.insert(ligated, [trafo.translate(0, 0)])
    c.insert(c17, [trafo.translate(r * 12., -r * 6.75)])
    ## 5C
    c.stroke(path.line(r * 12, -r * 7.5, r * 12, -r * 8.0), [style.linewidth.Thick, deco.earrow.normal])
    c.text(r * 12, -r * 8.25, "High-throughput", [text.halign.center, text.valign.middle, text.size(-4)])
    c.text(r * 12, -r * 8.25 - w * 0.5, "Sequencing", [text.halign.center, text.valign.middle, text.size(-4)])
    ## Labels
    c.text(0, -r * 10.65, "ChIA-PET", [text.halign.center, text.valign.top, text.size(0)])
    c.text(0, -r * 11.15, "many to many", [text.halign.center, text.valign.top, text.size(-3)])
    c.text(4 * r, -r * 10.65, "HiC", [text.halign.center, text.valign.top, text.size(0)])
    c.text(4 * r, -r * 11.15, "all to all", [text.halign.center, text.valign.top, text.size(-3)])
    c.text(8 * r, -r * 10.65, "4C", [text.halign.center, text.valign.top, text.size(0)])
    c.text(8 * r, -r * 11.15, "one to all", [text.halign.center, text.valign.top, text.size(-3)])
    c.text(12 * r, -r * 10.65, "5C", [text.halign.center, text.valign.top, text.size(0)])
    c.text(12 * r, -r * 11.15, "many to many", [text.halign.center, text.valign.top, text.size(-3)])
    c.text(16 * r, -r * 10.65, "3C", [text.halign.center, text.valign.top, text.size(0)])
    c.text(16 * r, -r * 11.15, "one to one", [text.halign.center, text.valign.top, text.size(-3)])
    ## Key
    c18 = canvas.canvas()
    c18.stroke(path.rect(0.05, -7.8 * r, 3.5 * r, 7.8 * r))
    c18.stroke(path.line(0.25 * r, -0.2 * r, 0.75 * r, -0.2 * r), [style.linewidth.THick, color.rgb.red])
    c18.stroke(path.line(0.25 * r, -0.4 * r, 0.75 * r, -0.4 * r), [style.linewidth.THick, color.rgb.blue])
    c18.text(1.0 * r, -0.3 * r, "RE Fragments", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.circle(0.5 * r, -0.9 * r, (y - r) * 1.2),
              [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c18.text(1.0 * r, -0.9 * r, "DNA Binding Protein", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.rect(0.4 * r, -1.5 * r - (y - r) * 1.7, 0.2 * r, (y - r) * 3.4),
               [style.linewidth.thick, deco.filled([color.gray(0.8)])])
    c18.text(1.0 * r, -1.5 * r, "Cross-linked Protein", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.insert(linker, [trafo.translate(0.5 * r, -2.1 * r)])
    c18.text(1.0 * r, -2.1 * r, "Sequencing Linker", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.insert(linker2, [trafo.translate(0.5 * r, -2.7 * r)])
    c18.text(1.0 * r, -2.7 * r, "MmeI-containing Linker", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.insert(antibody, [trafo.translate(0.5 * r, -3.3 * r)])
    c18.text(1.0 * r, -3.3 * r, "Antibody", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.insert(strep, [trafo.translate(0.5 * r, -3.9 * r)])
    c18.text(1.0 * r, -3.9 * r, "Streptavidin", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.line(0.5 * r, -4.5 * r, 0.5 * r, -4.7 * r), [style.linewidth.THick])
    c18.stroke(path.circle(0.5 * r, -4.5 * r + w * 0.16, w * 0.16),
              [style.linewidth.thick, deco.filled([color.rgb.green])])
    c18.text(1.0 * r, -4.5 * r, "Biotin", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.line(0.35 * r, -5.1 * r, 0.65 * r, -5.1 * r), [style.linewidth.THick, color.rgb.red])
    c18.stroke(path.line(0.5 * r, -5.1 * r - 0.14 * s, 0.5 * r, -5.1 * r + 0.14 * s),
               [style.linewidth.Thick])
    c18.text(1.0 * r, -5.1 * r, "RE Site", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.line(0.35 * r, -5.7 * r, 0.65 * r, -5.7 * r), [style.linewidth.THick, color.rgb.red])
    c18.insert(fragment, [trafo.translate(0.5 * r, -5.7 * r)])
    c18.text(1.0 * r, -5.7 * r, "Sonication Break Site", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.stroke(path.line(0.35 * r, -6.3 * r, 0.5 * r, -6.3 * r), [style.linewidth.THick, color.rgb.red])
    c18.stroke(path.line(0.65 * r, -6.3 * r, 0.5 * r, -6.3 * r), [style.linewidth.THick, color.rgb.blue])
    c18.insert(ligated, [trafo.translate(0.5 * r, -6.3 * r)])
    c18.text(1.0 * r, -6.3 * r, "Ligation Event", [text.halign.left, text.valign.middle, text.size(-4)])
    primer = path.path(path.moveto(r * 0.175, -w * 0.1), path.lineto(-r * 0.175, -w * 0.1),
            path.lineto(-r * 0.175, w * 0.025), path.lineto(0, w * 0.025),
            path.lineto(0, w * 0.1), path.closepath())
    c18.fill(primer, [style.linewidth.THIN, color.cmyk.RedOrange, trafo.translate(0.4 * r, -6.8 * r)])
    c18.fill(primer, [style.linewidth.THIN, color.cmyk.CornflowerBlue, trafo.scale(-1, 1),
             trafo.translate(0.6 * r, -7.0 * r)])
    c18.text(1.0 * r, -6.9 * r, "Targeted Primers", [text.halign.left, text.valign.middle, text.size(-4)])
    c18.fill(primer, [style.linewidth.THIN, color.cmyk.CadetBlue, trafo.translate(0.4 * r, -7.4 * r)])
    c18.fill(primer, [style.linewidth.THIN, color.cmyk.CadetBlue, trafo.scale(-1, 1),
             trafo.translate(0.6 * r, -7.6 * r)])
    c18.text(1.0 * r, -7.5 * r, "Universal Primers", [text.halign.left, text.valign.middle, text.size(-4)])
    c.insert(c18, [trafo.translate(r * 14.2, -r * 2.6)])
    return c


if __name__ == "__main__":
    main()
