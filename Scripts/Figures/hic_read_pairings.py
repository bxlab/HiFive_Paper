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

def sub_figA():
    c = canvas.canvas()
    subprotein1 = canvas.canvas()
    subprotein1.stroke(path.circle(0, 0, 0.3), [style.linewidth.Thick, deco.filled([color.cmyk.Mahogany])])
    subprotein2 = canvas.canvas()
    subprotein2.fill(path.circle(0, 0, 0.15), [color.cmyk.White])
    subprotein3 = canvas.canvas([canvas.clip(path.path(path.arc(0, 0.02, 0.16, 140, 400), path.lineto(-0.4, 0.1),
        path.lineto(-0.4, -0.4), path.lineto(0.4, -0.4), path.lineto(0.4, 0.1), path.closepath()))])
    subprotein3.stroke(path.path(path.arc(0, 0, 0.15, 180, 360)), [style.linewidth.Thick])
    subprotein4 = canvas.canvas([canvas.clip(path.path(path.arc(0, 0.1, 0.15, 160, 380), path.lineto(0.4, 0),
        path.lineto(0.4, -0.4), path.lineto(-0.4, -0.4), path.lineto(-0.4, 0), path.closepath()))])
    subprotein4.stroke(path.circle(0, 0, 0.3), [style.linewidth.Thick, deco.filled([color.cmyk.Mahogany])])
    protein1 = canvas.canvas()
    protein1.insert(subprotein1, [trafo.scale(1, 0.5)])
    protein1.insert(subprotein2, [trafo.translate(0, -0.1), trafo.scale(1, 0.5)])
    protein1.insert(subprotein3, [trafo.translate(0, 0.1), trafo.scale(1, 0.5)])
    protein2 = canvas.canvas()
    protein2.insert(subprotein4, [trafo.scale(1, 0.5)])
    protein2.insert(subprotein3, [trafo.translate(0, 0.1), trafo.scale(1, 0.5)])

    ligation = canvas.canvas()
    ligation.stroke(path.circle(0, 0, 0.1), [deco.filled([color.cmyk.MidnightBlue])])

    sequenced = canvas.canvas()
    sequenced.fill(path.circle(0, 0, 0.35), [color.gray(0.825)])

    top = canvas.canvas()
    top.stroke(path.line(0, 0, 0.1, 0), [style.linewidth.THIck, color.rgb.red])
    top.stroke(path.line(0.1, 0, 1.0, 0), [style.linewidth.THIck, color.rgb.blue])
    top.stroke(path.line(1.0, 0, 1.9, 0), [style.linewidth.THIck, color.cmyk.Dandelion])
    top.stroke(path.line(1.9, 0, 2.0, 0), [style.linewidth.THIck, color.rgb.red])
    top.stroke(path.line(1.09, 0, 1.1, 0), [deco.earrow.normal, color.rgb.blue])
    top.text(0, 0.7, r"Upstream fend", [text.halign.center, text.valign.top, text.size(-3)])
    top.stroke(path.line(0, 0.45, 0.5, 0.1), [deco.earrow.normal])
    top.text(2.0, 0.7, r"Downstream fend", [text.halign.center, text.valign.top, text.size(-3)])
    top.stroke(path.line(2.0, 0.45, 1.5, 0.1), [deco.earrow.normal])
    top.text(1.0, -0.5, r"RE sites", [text.halign.center, text.valign.top, text.size(-3)])
    top.stroke(path.path(path.moveto(0.15, -0.1), path.lineto(1.0, -0.45), path.lineto(1.85, -0.1)),
        [deco.earrow.normal, deco.barrow.normal])
    top.insert(protein1, [trafo.translate(4, 0)])
    top.insert(protein2, [trafo.translate(4, 0)])
    top.text(4.0, -0.2, r"Cross-linked", [text.halign.center, text.valign.top, text.size(-3)])
    top.text(4.0, -0.5, r"protein", [text.halign.center, text.valign.top, text.size(-3)])
    top.insert(ligation, [trafo.translate(5.5, 0)])
    top.text(5.5, -0.2, r"Ligation", [text.halign.center, text.valign.top, text.size(-3)])
    top.text(5.5, -0.5, r"event", [text.halign.center, text.valign.top, text.size(-3)])
    top.insert(sequenced, [trafo.translate(7.0, 0.2)])
    top.text(7.0, -0.2, r"Sequenced", [text.halign.center, text.valign.top, text.size(-3)])
    top.text(7.0, -0.5, r"DNA", [text.halign.center, text.valign.top, text.size(-3)])

    selfI = canvas.canvas()
    selfI.insert(sequenced, [trafo.translate(0, 0.6)])
    selfI.insert(protein1)
    selfI.stroke(path.path(path.moveto(0, -1.2), path.curveto(-0.5, -1.2, -0.075, -0.6, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    selfI.stroke(path.path(path.moveto(0, -1.2), path.curveto(0.5, -1.2, 0.075, -0.6, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    selfI.stroke(path.line(0.09, -1.2, 0.1, -1.2), [deco.earrow.normal, color.rgb.blue])
    selfI.insert(ligation, [trafo.translate(0, 0.6)])
    selfI.insert(protein2)

    adjacentI = canvas.canvas()
    adjacentI.insert(protein1)
    adjacentI.insert(sequenced, [trafo.translate(0, 0.6)])
    adjacentI.stroke(path.path(path.moveto(0, -1.2), path.curveto(-0.4, -1.2, -0.8, -0.8, -0.5, -0.5)),
        [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentI.stroke(path.path(path.moveto(0, -1.2), path.curveto(0.4, -1.2, 0.8, -0.8, 0.5, -0.5)),
        [style.linewidth.THIck, color.rgb.blue])
    adjacentI.stroke(path.line(-0.05, -1.2, 0.05, -1.2), [style.linewidth.THIck, color.rgb.red])
    adjacentI.stroke(path.path(path.moveto(-0.5, -0.5), path.curveto(-0.2, -0.2, -0.075, -0.3, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    adjacentI.stroke(path.path(path.moveto(0.5, -0.5), path.curveto(0.2, -0.2, 0.075, -0.3, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentI.stroke(path.line(-0.54, -0.54, -0.55, -0.55), [deco.earrow.normal, color.rgb.blue])
    adjacentI.stroke(path.line(0.46, -0.46, 0.45, -0.45), [deco.earrow.normal, color.rgb.blue])
    adjacentI.insert(ligation, [trafo.translate(0, 0.6)])
    adjacentI.insert(protein2)

    uncutI = canvas.canvas()
    uncutI.insert(protein1)
    uncutI.insert(sequenced, [trafo.translate(0, 0.6)])
    uncutI.stroke(path.path(path.moveto(-0.7, -1.2), path.curveto(-0.5, -1, -0.8, -0.8, -0.5, -0.5)),
        [style.linewidth.THIck, color.cmyk.Dandelion])
    uncutI.stroke(path.path(path.moveto(0.7, -1.2), path.curveto(0.5, -1, 0.8, -0.8, 0.5, -0.5)),
        [style.linewidth.THIck, color.rgb.blue])
    uncutI.stroke(path.path(path.moveto(-0.5, -0.5), path.curveto(-0.2, -0.2, -0.075, -0.3, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    uncutI.stroke(path.path(path.moveto(0.5, -0.5), path.curveto(0.2, -0.2, 0.075, -0.3, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    uncutI.stroke(path.line(-0.54, -0.54, -0.55, -0.55), [deco.earrow.normal, color.rgb.blue])
    uncutI.stroke(path.line(0.46, -0.46, 0.45, -0.45), [deco.earrow.normal, color.rgb.blue])
    uncutI.stroke(path.line(-0.05, 0.6, 0.05, 0.6), [style.linewidth.THIck, color.rgb.red])
    uncutI.insert(protein2)

    adjacentopp1I = canvas.canvas()
    adjacentopp1I.insert(protein1)
    adjacentopp1I.insert(sequenced, [trafo.translate(0, 0.6)])
    adjacentopp1I.stroke(path.path(path.moveto(-0.7, -1.2), path.curveto(-0.5, -1, -0.8, -0.8, -0.5, -0.5)),
        [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentopp1I.stroke(path.path(path.moveto(0.7, -1.2), path.curveto(0.5, -1, 0.8, -0.8, 0.5, -0.5)),
        [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentopp1I.stroke(path.path(path.moveto(-0.5, -0.5), path.curveto(-0.2, -0.2, -0.075, -0.3, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    adjacentopp1I.stroke(path.path(path.moveto(0.5, -0.5), path.curveto(0.2, -0.2, 0.075, -0.3, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    adjacentopp1I.stroke(path.line(-0.54, -0.54, -0.55, -0.55), [deco.earrow.normal, color.rgb.blue])
    adjacentopp1I.stroke(path.line(0.54, -0.54, 0.55, -0.55), [deco.earrow.normal, color.rgb.blue])
    adjacentopp1I.insert(ligation, [trafo.translate(0, 0.6)])
    adjacentopp1I.insert(protein2)

    adjacentopp2I = canvas.canvas()
    adjacentopp2I.insert(protein1)
    adjacentopp2I.insert(sequenced, [trafo.translate(0, 0.6)])
    adjacentopp2I.stroke(path.path(path.moveto(-0.7, -1.2), path.curveto(-0.5, -1, -0.8, -0.8, -0.5, -0.5)),
        [style.linewidth.THIck, color.rgb.blue])
    adjacentopp2I.stroke(path.path(path.moveto(0.7, -1.2), path.curveto(0.5, -1, 0.8, -0.8, 0.5, -0.5)),
        [style.linewidth.THIck, color.rgb.blue])
    adjacentopp2I.stroke(path.path(path.moveto(-0.5, -0.5), path.curveto(-0.2, -0.2, -0.075, -0.3, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentopp2I.stroke(path.path(path.moveto(0.5, -0.5), path.curveto(0.2, -0.2, 0.075, -0.3, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    adjacentopp2I.stroke(path.line(-0.46, -0.46, -0.45, -0.45), [deco.earrow.normal, color.rgb.blue])
    adjacentopp2I.stroke(path.line(0.46, -0.46, 0.45, -0.45), [deco.earrow.normal, color.rgb.blue])
    adjacentopp2I.insert(ligation, [trafo.translate(0, 0.6)])
    adjacentopp2I.insert(protein2)

    validI = canvas.canvas()
    validI.insert(protein1)
    validI.insert(sequenced, [trafo.translate(0, 0.6)])
    validI.stroke(path.path(path.moveto(-0.7, -1.2), path.curveto(-0.5, -1, -0.8, -0.8, -0.5, -0.5)),
        [style.linewidth.THIck, color.cmyk.Dandelion])
    validI.stroke(path.path(path.moveto(0.7, -1.2), path.curveto(0.5, -1, 0.8, -0.8, 0.5, -0.5)),
        [style.linewidth.THIck, color.rgb.blue])
    validI.stroke(path.path(path.moveto(-0.5, -0.5), path.curveto(-0.2, -0.2, -0.075, -0.3, -0.075, 0),
        path.curveto(-0.075, 0.4, -0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.rgb.blue])
    validI.stroke(path.path(path.moveto(0.5, -0.5), path.curveto(0.2, -0.2, 0.075, -0.3, 0.075, 0),
        path.curveto(0.075, 0.4, 0.4, 0.6, 0, 0.6)), [style.linewidth.THIck, color.cmyk.Dandelion])
    validI.stroke(path.line(-0.54, -0.54, -0.55, -0.55), [deco.earrow.normal, color.rgb.blue])
    validI.stroke(path.line(0.46, -0.46, 0.45, -0.45), [deco.earrow.normal, color.rgb.blue])
    validI.insert(ligation, [trafo.translate(0, 0.6)])
    validI.insert(protein2)

    bracket = canvas.canvas()
    bracket.stroke(path.path(path.moveto(-3.5, 0), path.arc(-3.2, -0.1, 0.3, 180, 270),
        path.arcn(-0.3, -0.7, 0.3, 90, 0), path.lineto(0, -0.8)),
        [deco.barrow.Small, deco.earrow.Small, style.linewidth.THick])
    bracket.stroke(path.path(path.moveto(3.5, 0), path.arcn(3.2, -0.1, 0.3, 0, -90),
        path.arc(0.3, -0.7, 0.3, 90, 180), path.lineto(0, -0.8)),
        [deco.barrow.Small, deco.earrow.Small, style.linewidth.THick])

    c.insert(top, [trafo.translate(3.0, 0.5)])
    c.insert(selfI, [trafo.translate(0, -2.2)])
    c.insert(uncutI, [trafo.translate(3, -2.2)])
    c.insert(adjacentI, [trafo.translate(5, -2.2)])
    c.insert(adjacentopp1I, [trafo.translate(8, -2.2)])
    c.insert(adjacentopp2I, [trafo.translate(10, -2.2)])
    c.insert(validI, [trafo.translate(13, -2.2)])
    c.insert(bracket, [trafo.translate(2.5, -4.7)])
    c.insert(bracket, [trafo.translate(10.5, -4.7)])

    c.text(0, -3.9, r"Self-interacting", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(0, -4.2, r"fragment", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(4, -3.9, r"Opposite strands,", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(4, -4.2, r"adjacent fragments", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(9, -3.9, r"Same strand,", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(9, -4.2, r"adjacent fragments", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(9, -4.5, r"or", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(9, -4.8, r"non-adjacent fragments", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(13, -3.9, r"Opposite strands,", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(13, -4.2, r"non-adjacent", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(13, -4.5, r"fragments", [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(2.5, -5.7, r"Excluded", [text.halign.center, text.valign.middle, text.size(-2)])
    c.text(10.5, -5.7, r"Included", [text.halign.center, text.valign.middle, text.size(-2)])
    return c



if __name__ == "__main__":
    main()
