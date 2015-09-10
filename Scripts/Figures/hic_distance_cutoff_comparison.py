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

method_colors = {
        'HiFive-Probability':color.cmyk.Black,
        'HiFive-Express':color.cmyk.CadetBlue,
        'HiFive-Binning':color.cmyk.MidnightBlue
}


def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    mm9_methods = {
        'HiFive-Probability':'%s/Analysis/hifive_mm9_ESC_prob_correlations.txt' % basedir, 
        'HiFive-Express':'%s/Analysis/hifive_mm9_ESC_exp_correlations.txt' % basedir, 
        'HiFive-Binning':'%s/Analysis/hifive_mm9_ESC_bin_correlations.txt' % basedir, 
    }
    dist_methods = {
        'HiFive-Probability':'%s/Analysis/hifive_mm9_ESC_probnodist_correlations.txt' % basedir, 
        'HiFive-Express':'%s/Analysis/hifive_mm9_ESC_expnodist_correlations.txt' % basedir, 
        'HiFive-Binning':'%s/Analysis/hifive_mm9_ESC_binnodist_correlations.txt' % basedir, 
    }
    mm9_data = load_data(mm9_methods)
    dist_data = load_data(dist_methods)
    width = 16.8
    spacer = 0.4
    range_width = (width - spacer) / 3.0
    c = canvas.canvas()
    mm9_ranges_img = plot_dataset_ranges(mm9_data, dist_data, range_width)
    c.insert(mm9_ranges_img)
    mm9_overall_img = plot_overall(mm9_data, dist_data, range_width, range_width)
    mm9_overall_img.text(0, range_width, 'b',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(mm9_overall_img, [trafo.translate(range_width * 2 + spacer, range_width)])
    c.insert(plot_key(range_width / 2, range_width / 2),
             [trafo.translate(range_width * 2.25 + spacer, range_width * 0.25)])
    c.writePDFfile(out_fname)

def load_data(fdict):
    all_data = {}
    for name in fdict.keys():
        if not os.path.exists(fdict[name]):
            continue
        data = []
        for line in open(fdict[name]):
            temp = line[:-1].split('\t')
            data.append((int(temp[0]), int(temp[1].replace('NA', '-1')), temp[2], int(temp[3]), float(temp[4])))
        all_data[name] = numpy.array(data, dtype=numpy.dtype([('binsize', numpy.int32), ('range', numpy.int32),
                                     ('interaction', 'a5'), ('count', numpy.int32), ('correlation', numpy.float32)]))
    return all_data

def plot_overall(data0, data1, width, height):
    vo = 0.55
    ho = 1.1
    plot_width = width - ho
    plot_height = height - vo
    c = canvas.canvas()
    methods = data0.keys()
    methods.sort()
    bar_colors = []
    cis_binsizes = numpy.unique(data0[methods[0]]['binsize'][numpy.where(data0[methods[0]]['interaction'] == 'cis')])
    trans_binsizes = numpy.unique(data0[methods[0]]['binsize'][numpy.where(data0[methods[0]]['interaction'] == 'trans')])
    ymin = numpy.inf
    ymax = -numpy.inf
    Y = numpy.zeros((len(methods), cis_binsizes.shape[0] + trans_binsizes.shape[0]), dtype=numpy.float32)    
    for i, method in enumerate(methods):
        for j, binsize in enumerate(cis_binsizes):
            where = numpy.where((data0[method]['binsize'] == binsize) *
                                (data0[method]['interaction'] == 'cis') *
                                (data0[method]['range'] < 0))
            where1 = numpy.where((data1[method]['binsize'] == binsize) *
                                 (data1[method]['interaction'] == 'cis') *
                                 (data1[method]['range'] < 0))
            if where[0].shape[0] > 0:
                Y[i, j] = (data1[method]['correlation'][where1] - data0[method]['correlation'][where])
                ymin = min(ymin, Y[i, j])
                ymax = max(ymax, Y[i, j])
        for j, binsize in enumerate(trans_binsizes):
            where = numpy.where((data0[method]['binsize'] == binsize) *
                                (data0[method]['interaction'] == 'trans') *
                                (data0[method]['range'] < 0))
            where1 = numpy.where((data1[method]['binsize'] == binsize) *
                                 (data1[method]['interaction'] == 'trans') *
                                 (data1[method]['range'] < 0))
            if where[0].shape[0] > 0:
                Y[i, j + cis_binsizes.shape[0]] = (data1[method]['correlation'][where1] -
                                                              data0[method]['correlation'][where])
                ymin = min(ymin, Y[i, j + cis_binsizes.shape[0]])
                ymax = max(ymax, Y[i, j + cis_binsizes.shape[0]])
        bar_colors.append(method_colors[method])
    yspan = ymax - ymin
    ymin -= yspan * 0.05
    ymax += yspan * 0.05
    Y = numpy.array(Y)
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.nestedbar(painter=graph.axis.painter.bar(nameattrs=None)),
                      y=graph.axis.lin(painter=painter, min=ymin, max=ymax),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    y0 = plot_height * (-ymin) / (ymax - ymin)
    g.stroke(path.line(0, plot_height * (-ymin) / (ymax - ymin), plot_width, plot_height * (-ymin) / (ymax - ymin)),
        [style.linestyle.dotted, style.linewidth.THin])
    w0 = plot_width / Y.shape[1]
    w1 = w0 / (len(methods) + 0.5)
    for i in range(len(methods)):
        for j in range(Y.shape[1]):
            x = j * w0 + (i + 0.25) * w1
            y = plot_height * (Y[i, j] - ymin) / (ymax - ymin)
            g.stroke( path.rect(x, y0, w1, y - y0), [deco.filled([method_colors[methods[i]]])])
    c.insert(g, [trafo.translate(ho, vo)])
    for i, label in enumerate(["10Kb", "50Kb", "250Kb", "1Mb", "250Kb", "1Mb"]):
        c.text(ho + plot_width * (i + 0.5) / 6.0, vo - 0.05, "%s" % label,
               [text.halign.center, text.valign.top, text.size(-3)])
    c.text(ho + plot_width * 2.0 / 6.0, 0.05, "cis",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2, vo * 0.5, ho - 0.2 + plot_width * 4.0 / 6.0, vo * 0.5), [style.linewidth.THin])
    c.text(ho + plot_width * 5.0 / 6.0, 0.05, "trans",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2 + plot_width * 4.0 / 6.0, vo * 0.5, ho - 0.2 + plot_width, vo * 0.5), [style.linewidth.THin])
    c.text(0, plot_height * 0.5 + vo, r"$r_{0K} - r_{500K}$",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
    return c

def plot_dataset_ranges(data0, data1, width):
    methods = data0.keys()
    binsizes = numpy.unique(data0[methods[0]]['binsize'])
    ho = 0.8
    ho2 = 0.25
    vo = 0.7
    vo2 = 0.4
    spacer = 0.0
    plot_width = (2 * width - ho * 2 - ho2 - spacer) / 2
    plot_height = width - vo - vo2
    c = canvas.canvas()
    for i, binsize in enumerate(binsizes):
        img = plot_single_range(data0, data1, binsize, plot_width, plot_width, vo, vo2)
        c.insert(img, [trafo.translate((plot_width + spacer + ho) * (i % 2) + ho2 + ho,
                                        (1 - i /
                                         2) * (plot_height + spacer + vo + vo2) + vo)])
    c.text(0, plot_height + 0.5 * spacer + vo + vo2, r"$r_{0K} - r_{500K}$",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
    c.text(plot_width + ho2 + spacer * 0.5 + ho, 0, "Interaction Range (bp)",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.text(0, width * 2 + spacer, 'a',
                        [text.halign.left, text.valign.top, text.size(-1)])
    return c

def plot_single_range(data0, data1, binsize, width, height, vo, vo2):
    plot_width = width
    plot_height = height
    c = canvas.canvas()
    xmax = 0.0
    methods = data0.keys()
    methods.sort()
    for method in methods:
        where = numpy.where((data0[method]['binsize'] == binsize) *
                            (data0[method]['interaction'] == 'cis') *
                            (data0[method]['range'] >= 0))
        if where[0].shape[0] > 0:
            xmax = max(xmax, numpy.amax(data0[method]['range'][where]))
            X = data0[method]['range'][where]
    X = numpy.r_[0, X]
    X[0] = X[1] ** 2.0 / X[2]
    xmin = X[0]
    ymin = numpy.inf
    ymax = -numpy.inf
    binsizes = numpy.unique(data0[methods[0]]['binsize'])
    for method in methods:
        for b in binsizes:
            where = numpy.where((data0[method]['binsize'] == b) *
                                (data0[method]['interaction'] == 'cis') *
                                (data0[method]['range'] >= 0))
            where1 = numpy.where((data1[method]['binsize'] == b) *
                                 (data1[method]['interaction'] == 'cis') *
                                 (data1[method]['range'] >= 0))
            if where[0].shape[0] > 0:
                Y = (data1[method]['correlation'][where1] -
                                 data0[method]['correlation'][where] )
                ymin = min(ymin, numpy.amin(Y))
                ymax = max(ymax, numpy.amax(Y))
    yspan = ymax - ymin
    ymin -= yspan * 0.05
    ymax += yspan * 0.05
    g = graph.graphxy(width=plot_width, height=plot_height,
                      x=graph.axis.log(painter=painter, min=X[0], max=xmax),
                      y=graph.axis.lin(painter=painter, min=ymin, max=ymax),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    for x in X[1:-1]:
        pos = ((log(x) - log(xmin)) / (log(xmax) - log(xmin)) * plot_width)
        g.stroke(path.line(pos, 0, pos, plot_height), [style.linestyle.dotted, style.linewidth.THin])
    X = (X[1:] ** 0.5) * (X[:-1] ** 0.5)
    for method in methods:
        where = numpy.where((data0[method]['binsize'] == binsize) *
                            (data0[method]['interaction'] == 'cis') *
                            (data0[method]['range'] > 0))
        where1 = numpy.where((data1[method]['binsize'] == binsize) *
                             (data1[method]['interaction'] == 'cis') *
                             (data1[method]['range'] > 0))
        if where[0].shape[0] > 0:
            Y = ( data1[method]['correlation'][where1] -
                             data0[method]['correlation'][where] )
            g.plot(graph.data.points(zip(X, Y), x=1, y=2),
                   [graph.style.line(lineattrs=[method_colors[method], style.linewidth.Thick])])
    g.stroke(path.line(0, plot_height * (-ymin) / (ymax - ymin), plot_width, plot_height * (-ymin) / (ymax - ymin)),
        [style.linestyle.dotted, style.linewidth.THin])
    if binsize / 1000000 > 0:
        binstring = "%iMb" % (binsize / 1000000)
    elif binsize / 1000 > 0:
        binstring = "%iKb" % (binsize / 1000)
    else:
        binstring = str(binsize)
    g.text(plot_width / 2, plot_height + vo2, "%s binning" % (binstring),
           [text.halign.center, text.valign.top, text.size(-2)])
    c.insert(g)
    return c

def plot_key(width, height):
    c = canvas.canvas()
    step = height / float(len(method_colors))
    for i, meth in enumerate(['HiFive-Probability', 'HiFive-Express', 'HiFive-Binning']):
        c.fill(path.rect(0.2, height - step * (i + 0.5) - 0.1, 0.2, 0.2),
               [method_colors[meth]])
        c.text(0.5, height - step * (i + 0.5), meth, [text.halign.left, text.valign.middle, text.size(-2)])
    return c
if __name__ == "__main__":
    main()
