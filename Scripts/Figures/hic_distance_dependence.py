#!/usr/bin/env python

import sys
import os

import numpy
from math import log, ceil
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style

import hifive

unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )

method_colors = {
        'raw':color.cmyk.Black,
        'distance-corrected':color.rgb.red,
}


def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    mm9_methods = {
        'distance-corrected':'%s/Analysis/hifive_mm9_ESC_exp_correlations.txt' % basedir, 
        'raw':'%s/Analysis/hifive_mm9_ESC_expdist_correlations.txt' % basedir, 
    }
    mm9_data = load_data(mm9_methods)
    width = 16.8  
    spacer = 0.4
    plot_width = (width - spacer * 4) / 5
    c = canvas.canvas()
    mm9_ranges_img, mm9_ranges_height = plot_dataset_ranges(mm9_data, plot_width, spacer)
    mm9_ranges_img.text(0, mm9_ranges_height, 'b',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(mm9_ranges_img, [trafo.translate(plot_width + spacer, 0)])
    mm9_overall_img = plot_overall(mm9_data, plot_width - spacer, (mm9_ranges_height - spacer) - 0.3)
    c.text(0,  mm9_ranges_height, 'a',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(mm9_overall_img, [trafo.translate(0, 0.4)])
    c.writePDFfile(out_fname)

def load_data(fdict):
    all_data = {}
    for name in fdict.keys():
        if not os.path.exists(fdict[name]):
            continue
        data = []
        for line in open(fdict[name]):
            temp = line[:-1].split('\t')
            data.append((int(temp[0]), int(temp[1].replace('NA', '0')), temp[2], int(temp[3]), float(temp[4])))
        all_data[name] = numpy.array(data, dtype=numpy.dtype([('binsize', numpy.int32), ('range', numpy.int32),
                                     ('interaction', 'a5'), ('count', numpy.int32), ('correlation', numpy.float32)]))
    return all_data

def plot_overall(data, width, height):
    plot_width = width - 0.4
    plot_height = height - 0.4
    c = canvas.canvas()
    methods = data.keys()
    methods.sort()
    bar_colors = []
    cis_binsizes = numpy.unique(data[methods[0]]['binsize'][numpy.where(data[methods[0]]['interaction'] == 'cis')])
    trans_binsizes = numpy.unique(data[methods[0]]['binsize'][numpy.where(data[methods[0]]['interaction'] == 'trans')])
    Y = numpy.zeros((len(methods), cis_binsizes.shape[0] + trans_binsizes.shape[0]), dtype=numpy.float32)    
    for i, method in enumerate(methods):
        for j, binsize in enumerate(cis_binsizes):
            where = numpy.where((data[method]['binsize'] == binsize) *
                                (data[method]['interaction'] == 'cis') *
                                (data[method]['range'] == 0))
            if where[0].shape[0] > 0:
                Y[i, j] = data[method]['correlation'][where]
        for j, binsize in enumerate(trans_binsizes):
            where = numpy.where((data[method]['binsize'] == binsize) *
                                (data[method]['interaction'] == 'trans') *
                                (data[method]['range'] == 0))
            if where[0].shape[0] > 0:
                Y[i, j + cis_binsizes.shape[0]] = data[method]['correlation'][where]
        bar_colors.append(method_colors[method])
    Y = numpy.array(Y)
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.nestedbar(painter=graph.axis.painter.bar(nameattrs=None)),
                      y=graph.axis.lin(painter=painter),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    for i in range(len(methods)):
        g.plot(graph.data.points(zip(zip(range(Y.shape[1]), [i] * Y.shape[1]), Y[i, :]), xname=1, y=2),
               [graph.style.changebar([method_colors[methods[i]]])])
    step = plot_width / (cis_binsizes.shape[0] + trans_binsizes.shape[0])
    for i, binsize in enumerate(cis_binsizes):
        g.text(step * (0.5 + i), -0.05, "%s cis" % (str(binsize/1000) + 'Kb').replace('000Kb', 'Mb'),
               [text.halign.right, text.valign.middle, text.size(-4), trafo.rotate(45)])
    for i, binsize in enumerate(trans_binsizes):
        g.text(step * (0.5 + i + cis_binsizes.shape[0]), -0.05, "%s trans" % (str(binsize/1000) + 'Kb').replace('000Kb', 'Mb'),
               [text.halign.right, text.valign.middle, text.size(-4), trafo.rotate(45)])
    c.insert(g, [trafo.translate(0.7, 0.4)])
    c.text(0, plot_height / 2.0 + 0.4, "Dataset Correlation",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    return c

def plot_dataset_ranges(data, width, spacer):
    methods = data.keys()
    binsizes = numpy.unique(data[methods[0]]['binsize'])
    plot_width = width - 0.4
    plot_height = plot_width + 0.2
    c = canvas.canvas()
    for i, binsize in enumerate(binsizes):
        img = plot_single_range(data, binsize, plot_width, plot_height)
        c.insert(img, [trafo.translate((plot_width + spacer) * i + 0.4, 0.4 + spacer * 0.5)])
    c.text(0, plot_height * 0.5 + 0.2 + spacer * 0.75, "Dataset correlation",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    c.text(plot_width * 2 + spacer * 1.5 + 0.4, 0.35, "Interaction range (bp)",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.fill(path.rect(0.4, 0.025, 0.2, 0.2), [color.rgb.black])
    c.text(0.7, 0.125, "HiFive-Express using raw reads", [text.halign.left, text.valign.middle, text.size(-2)])
    c.fill(path.rect(plot_width * 2 + spacer * 1.5, 0.025, 0.2, 0.2), [color.rgb.red])
    c.text(plot_width * 2 + spacer * 1.5 + 0.3, 0.125, "HiFive-Express using distance-corrected reads",
           [text.halign.left, text.valign.middle, text.size(-2)])
    return c, plot_height + 0.4 + spacer

def plot_single_range(data, binsize, width, height):
    plot_width = width - 0.4
    plot_height = height - 0.8
    c = canvas.canvas()
    xmax = 0.0
    methods = data.keys()
    methods.sort()
    for method in methods:
        where = numpy.where((data[method]['binsize'] == binsize) *
                            (data[method]['interaction'] == 'cis') *
                            (data[method]['range'] > 0))
        if where[0].shape[0] > 0:
            xmax = max(xmax, numpy.amax(data[method]['range'][where]))
            X = data[method]['range'][where]
    X = numpy.r_[0, X]
    X[0] = X[1] ** 2.0 / X[2]
    xmin = X[0]
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.log(painter=painter, min=X[0], max=xmax),
                      y=graph.axis.lin(painter=painter),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    for x in X[1:-1]:
        pos = ((log(x) - log(xmin)) / (log(xmax) - log(xmin)) * plot_width)
        g.stroke(path.line(pos, 0, pos, plot_height), [style.linestyle.dotted, style.linewidth.THin])

    X = (X[1:] ** 0.5) * (X[:-1] ** 0.5)
    for method in methods:
        where = numpy.where((data[method]['binsize'] == binsize) *
                            (data[method]['interaction'] == 'cis') *
                            (data[method]['range'] > 0))
        if where[0].shape[0] > 0:
            Y = data[method]['correlation'][where]
            g.plot(graph.data.points(zip(X, Y), x=1, y=2),
                   [graph.style.line(lineattrs=[method_colors[method], style.linewidth.Thick])])
    if binsize / 1000000 > 0:
        binstring = "%iMb" % (binsize / 1000000)
    elif binsize / 1000 > 0:
        binstring = "%iKb" % (binsize / 1000)
    else:
        binstring = str(binsize)
    g.text(plot_width / 2, plot_height + 0.05, "%s Binning" % (binstring),
           [text.halign.center, text.valign.bottom, text.size(-2)])
    c.insert(g, [trafo.translate(0.4, 0.4)])
    return c

if __name__ == "__main__":
    main()
