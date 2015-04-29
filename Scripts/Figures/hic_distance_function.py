#!/usr/bin/env python

import sys
from pyx import *
import hifive
import numpy

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular(labelattrs=[text.size(-3)], titleattrs=[text.size(-2)])


basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
out_fname = sys.argv[1]
hifive_fname = "%s/data/HiC/HiFive/mm9_ESC_HindIII.hcp" % basedir

hic = hifive.HiC(hic_fname, 'r')
params = numpy.copy(hic.distance_parameters)
points = numpy.zeros((params.shape[0] + 1, 2), dtype=numpy.float32)
points[1:-1, 0] = params[:-1, 0]
points[0, 0] = 2 * points[1, 0] - points[2, 0]
points[-1, 0] = 2 * points[-2, 0] - points[-3, 0]
points[:-1, 1] = points[:-1, 0] * params[:, 1] + params[:, 2]
points[-1, 1] = points[-1, 0] * params[-1, 1] + params[-1, 2]
points = numpy.exp(points)

params2 = numpy.copy(hic.bin_distance_parameters)
points2 = numpy.zeros((params2.shape[0] + 1, 2), dtype=numpy.float32)
points2[1:-1, 0] = params2[:-1, 0]
points2[0, 0] = 2 * points2[1, 0] - points2[2, 0]
points2[-1, 0] = 2 * points2[-2, 0] - points2[-3, 0]
points2[:-1, 1] = points2[:-1, 0] * params2[:, 1] + params2[:, 2]
points2[-1, 1] = points2[-1, 0] * params2[-1, 1] + params2[-1, 2]
points2 = numpy.exp(points2)

xmin = min(numpy.amin(points[:, 0]), numpy.amin(points2[:, 0]))
xmax = max(numpy.amax(points[:, 0]), numpy.amax(points2[:, 0]))
ymin = min(numpy.amin(points[:, 1]), numpy.amin(points2[:, 1]))
ymax = max(numpy.amax(points[:, 1]), numpy.amax(points2[:, 1]))
xmin2 = numpy.exp(numpy.log(xmin) + (numpy.log(xmax) - numpy.log(xmin)) * (4.5 / 6.0))
xmax2 = numpy.exp(numpy.log(xmin) + (numpy.log(xmax) - numpy.log(xmin)) * (5.1 / 6.0))
ymin2 = numpy.exp(numpy.log(ymin) + (numpy.log(ymax) - numpy.log(ymin)) * (1.3 / 6.0))
ymax2 = numpy.exp(numpy.log(ymin) + (numpy.log(ymax) - numpy.log(ymin)) * (1.9 / 6.0))
c = canvas.canvas()
g = graph.graphxy( width=6.0, height=6.0, 
    x=graph.axis.log(min=xmin, max=xmax, title='', painter=painter), 
    y=graph.axis.log(min=ymin, max=ymax, title='', painter=painter), 
    x2=graph.axis.lin(min=0, max=1, parter=None),
    y2=graph.axis.lin(min=0, max=1, parter=None))
for i in range(points.shape[0] - 1):
    g.plot(graph.data.function("y(x)=exp(log(x)*%f+%f)" % (params[i, 1], params[i, 2]), min=points[i, 0], max=points[i + 1, 0]), [graph.style.line([style.linewidth.Thick])])
    g.plot(graph.data.function("y(x)=exp(log(x)*%f+%f)" % (params2[i, 1], params2[i, 2]), min=points2[i, 0], max=points2[i + 1, 0]), [graph.style.line([style.linewidth.Thick, color.rgb.red])])
g.stroke(path.rect(4.5, 1.3, 0.6, 0.6), [color.rgb.red])
g.text(3.0, -0.8, "Inter-fragment Distance (Kb)", [text.size(-2), text.halign.center, text.valign.top])
g.text(-1.0, 3.0, "Mean Interaction Count", [text.size(-2), text.halign.center, text.valign.bottom, trafo.rotate(90)])
c.insert(g)
c.stroke(path.line(5.1, 1.9, 8, 6), [color.rgb.red])
c.stroke(path.line(5.1, 1.3, 8, 0), [color.rgb.red])
g = graph.graphxy( width=6.0, height=6.0, 
    x=graph.axis.log(min=xmin2, max=xmax2, title='', painter=painter), 
    y=graph.axis.log(min=ymin2, max=ymax2, title='', painter=painter), 
    x2=graph.axis.lin(min=0, max=1, parter=None),
    y2=graph.axis.lin(min=0, max=1, parter=None))
line_colors = [color.rgb.black, color.rgb.red]
for i in range(points.shape[0] - 1):
    g.plot(graph.data.function("y(x)=exp(log(x)*%f+%f)" % (params[i, 1], params[i, 2]), min=points[i, 0], max=points[i + 1, 0]), [graph.style.line([style.linewidth.Thick, line_colors[i%2]])])
    g.plot(graph.data.function("y(x)=exp(log(x)*%f+%f)" % (params2[i, 1], params2[i, 2]), min=points2[i, 0], max=points2[i + 1, 0]), [graph.style.line([style.linewidth.Thick, line_colors[i%2]])])
g.text(3.0, -0.8, "Inter-fragment Distance (Kb)", [text.size(-2), text.halign.center, text.valign.top])
g.text(-1.0, 3.0, "Mean Interaction Count", [text.size(-2), text.halign.center, text.valign.bottom, trafo.rotate(90)])
c.insert(g, [trafo.translate(8, 0)])
c.writePDFfile(out_fname)
