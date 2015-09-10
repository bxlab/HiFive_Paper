#!/usr/bin/env python

import sys
import os

import numpy
from math import log
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-2)] )

def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    bin_file = '%s/Analysis/hifive_mm9_ESC_prob_correlations.txt' % basedir
    pois_file = '%s/Analysis/hifive_mm9_ESC_probpois_correlations.txt' % basedir
    bin_data = load_data(bin_file)
    pois_data = load_data(pois_file)
    width = 16.8
    spacer = 0.4
    range_width = (width - spacer) / 5.0
    range_height = range_width * 1.4
    c = canvas.canvas()
    ranges_img = plot_dataset_ranges(bin_data, pois_data, range_width * 4, range_height)
    c.insert(ranges_img)
    overall_img = plot_overall(bin_data, pois_data, range_width, range_height)
    c.insert(overall_img, [trafo.translate(range_width * 4 + spacer, 0)])
    c.writePDFfile(out_fname)

def load_data(fname):
    data = []
    for line in open(fname):
        temp = line[:-1].split('\t')
        data.append((int(temp[0]), int(temp[1].replace('NA', '-1')), temp[2], int(temp[3]), float(temp[4])))
    all_data = numpy.array(data, dtype=numpy.dtype([('binsize', numpy.int32), ('range', numpy.int32),
                                 ('interaction', 'a5'), ('count', numpy.int32), ('correlation', numpy.float32)]))
    return all_data

def plot_overall(data0, data1, width, height):
    vo = 1.15
    ho = 1.15
    plot_width = width - ho
    plot_height = height - vo
    c = canvas.canvas()
    cis_binsizes = numpy.unique(data0['binsize'][numpy.where(data0['interaction'] == 'cis')])
    trans_binsizes = numpy.unique(data0['binsize'][numpy.where(data0['interaction'] == 'trans')])
    ymin = numpy.inf
    ymax = -numpy.inf
    Y = numpy.zeros((cis_binsizes.shape[0] + trans_binsizes.shape[0]), dtype=numpy.float32)    
    for j, binsize in enumerate(cis_binsizes):
        where = numpy.where((data0['binsize'] == binsize) *
                            (data0['interaction'] == 'cis') *
                            (data0['range'] < 0))
        where1 = numpy.where((data1['binsize'] == binsize) *
                             (data1['interaction'] == 'cis') *
                             (data1['range'] < 0))
        if where[0].shape[0] > 0:
            Y[j] = (data1['correlation'][where1] - data0['correlation'][where])
            ymin = min(ymin, Y[j])
            ymax = max(ymax, Y[j])
    for j, binsize in enumerate(trans_binsizes):
        where = numpy.where((data0['binsize'] == binsize) *
                            (data0['interaction'] == 'trans') *
                            (data0['range'] < 0))
        where1 = numpy.where((data1['binsize'] == binsize) *
                             (data1['interaction'] == 'trans') *
                             (data1['range'] < 0))
        if where[0].shape[0] > 0:
            Y[j + cis_binsizes.shape[0]] = (data1['correlation'][where1] - data0['correlation'][where])
            ymin = min(ymin, Y[j + cis_binsizes.shape[0]])
            ymax = max(ymax, Y[j + cis_binsizes.shape[0]])
    yspan = ymax - ymin
    ymin -= yspan * 0.05
    ymax += yspan * 0.05
    Y = numpy.array(Y)
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.bar(painter=graph.axis.painter.bar(nameattrs=None)),
                      y=graph.axis.lin(painter=painter, min=ymin, max=ymax),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    y0 = plot_height * (-ymin) / (ymax - ymin)
    g.stroke(path.line(0, plot_height * (-ymin) / (ymax - ymin), plot_width, plot_height * (-ymin) / (ymax - ymin)),
        [style.linestyle.dotted, style.linewidth.THin])
    w0 = plot_width / Y.shape[0]
    w1 = w0 / 1.5
    for j in range(Y.shape[0]):
        x = j * w0 + 0.25 * w1
        y = plot_height * (Y[j] - ymin) / (ymax - ymin)
        g.fill( path.rect(x, y0, w1, y - y0) )
    c.insert(g, [trafo.translate(ho, vo)])
    for i, label in enumerate(["10Kb", "50Kb", "250Kb", "1Mb", "250Kb", "1Mb"]):
        c.text(ho + plot_width * (i + 0.5) / 6.0, vo - 0.05, "%s" % label,
               [text.halign.right, text.valign.middle, text.size(-3), trafo.rotate(90)])
    c.text(ho + plot_width * 2.0 / 6.0, 0, "cis",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2, 0.3, ho - 0.2 + plot_width * 4.0 / 6.0, 0.3), [style.linewidth.THin])
    c.text(ho + plot_width * 5.0 / 6.0, 0, "trans",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2 + plot_width * 4.0 / 6.0, 0.3, ho - 0.2 + plot_width, 0.3), [style.linewidth.THin])
    c.text(0, plot_height * 0.5 + vo, r"$r_{Poisson} - r_{binomial}$",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
    c.text(0, height, 'b', [text.halign.left, text.valign.top, text.size(-1)])
    return c

def plot_dataset_ranges(data0, data1, width, height):
    binsizes = numpy.unique(data0['binsize'])
    ho = 1.1
    ho2 = 0.2
    vo = 1.0
    plot_width = (width - ho - ho2 * 3) / 4.0
    plot_height = height - vo
    c = canvas.canvas()
    for i, binsize in enumerate(binsizes):
        if i == 0:
            ylabel = True
        else:
            ylabel = False
        img = plot_single_range(data0, data1, binsize, plot_width, plot_height, ylabel)
        c.insert(img, [trafo.translate(ho + i * (plot_width + ho2), vo - 0.3)])
    c.text(0, plot_height * 0.5 + vo - 0.3, r"$r_{Poisson} - r_{binomial}$",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
    c.text(plot_width * 2 + ho + 1.5 * ho2, 0, "Interaction Range (bp)",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.text(0, height, 'a', [text.halign.left, text.valign.top, text.size(-1)])
    return c

def plot_single_range(data0, data1, binsize, width, height, ylabel):
    plot_width = width
    plot_height = height
    c = canvas.canvas()
    xmax = 0.0
    where = numpy.where((data0['binsize'] == binsize) *
                        (data0['interaction'] == 'cis') *
                        (data0['range'] >= 0))
    if where[0].shape[0] > 0:
        xmax = max(xmax, numpy.amax(data0['range'][where]))
        X = data0['range'][where]
    X = numpy.r_[0, X]
    X[0] = X[1] ** 2.0 / X[2]
    xmin = X[0]
    ymin = numpy.inf
    ymax = -numpy.inf
    binsizes = numpy.unique(data0['binsize'])
    for b in binsizes:
        where = numpy.where((data0['binsize'] == b) *
                            (data0['interaction'] == 'cis') *
                            (data0['range'] >= 0))
        where1 = numpy.where((data1['binsize'] == b) *
                             (data1['interaction'] == 'cis') *
                             (data1['range'] >= 0))
        if where[0].shape[0] > 0:
            Y = (data1['correlation'][where1] -
                             data0['correlation'][where] )
            ymin = min(ymin, numpy.amin(Y))
            ymax = max(ymax, numpy.amax(Y))
    yspan = ymax - ymin
    ymin -= yspan * 0.05
    ymax += yspan * 0.05
    if ylabel:
        yaxis = graph.axis.lin(painter=painter, min=ymin, max=ymax)
    else:
        yaxis = graph.axis.lin(painter=None, min=ymin, max=ymax)
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.log(painter=painter, min=X[0], max=xmax),
                      y=yaxis,
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    g.stroke( path.line(0, 0, 0, plot_height) )
    for x in X[1:-1]:
        pos = ((log(x) - log(xmin)) / (log(xmax) - log(xmin)) * plot_width)
        g.stroke(path.line(pos, 0, pos, plot_height), [style.linestyle.dotted, style.linewidth.THin])
    X = (X[1:] ** 0.5) * (X[:-1] ** 0.5)
    where = numpy.where((data0['binsize'] == binsize) *
                        (data0['interaction'] == 'cis') *
                        (data0['range'] > 0))
    where1 = numpy.where((data1['binsize'] == binsize) *
                         (data1['interaction'] == 'cis') *
                         (data1['range'] > 0))
    if where[0].shape[0] > 0:
        Y = ( data1['correlation'][where1] - data0['correlation'][where] )
        g.plot(graph.data.points(zip(X, Y), x=1, y=2),
               [graph.style.line(lineattrs=[style.linewidth.Thick])])
    g.stroke(path.line(0, plot_height * (-ymin) / (ymax - ymin), plot_width, plot_height * (-ymin) / (ymax - ymin)),
        [style.linestyle.dotted, style.linewidth.THin])
    if binsize / 1000000 > 0:
        binstring = "%iMb" % (binsize / 1000000)
    elif binsize / 1000 > 0:
        binstring = "%iKb" % (binsize / 1000)
    else:
        binstring = str(binsize)
    g.text(plot_width / 2, plot_height + 0.3, "%s binning" % (binstring),
           [text.halign.center, text.valign.top, text.size(-2)])
    c.insert(g)
    return c

if __name__ == "__main__":
    main()
