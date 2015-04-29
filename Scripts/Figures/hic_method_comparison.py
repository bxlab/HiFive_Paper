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
        'HiFive-Binning':color.cmyk.MidnightBlue, 
        'HiCNorm':color.cmyk.Dandelion,
        'HiCPipe':color.cmyk.Mahogany,
        'Matrix-Balancing':color.cmyk.OliveGreen,  
}


def main():
    out_fname = sys.argv[1]
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    mm9_methods = {
        'HiFive-Probability':'%s/Analysis/hifive_mm9_ESC_prob_correlations.txt' % basedir, 
        'HiFive-Express':'%s/Analysis/hifive_mm9_ESC_exp_correlations.txt' % basedir, 
        'HiFive-Binning':'%s/Analysis/hifive_mm9_ESC_bin_correlations.txt' % basedir, 
        'HiCNorm':'%s/Analysis/hicnorm_mm9_ESC_correlations.txt' % basedir, 
        'HiCPipe':'%s/Analysis/hicpipe_mm9_ESC_correlations.txt' % basedir, 
        'Matrix-Balancing':'%s/Analysis/mb_mm9_ESC_correlations.txt' % basedir,
    }
    hg19_methods = {
        'HiFive-Probability':'%s/Analysis/hifive_hg19_GM12878_prob_correlations.txt' % basedir, 
        'HiFive-Express':'%s/Analysis/hifive_hg19_GM12878_exp_correlations.txt' % basedir, 
        'HiFive-Binning':'%s/Analysis/hifive_hg19_GM12878_bin_correlations.txt' % basedir, 
        'HiCNorm':'%s/Analysis/hicnorm_hg19_GM12878_correlations.txt' % basedir, 
        'HiCPipe':'%s/Analysis/hicpipe_hg19_GM12878_correlations.txt' % basedir, 
        'Matrix-Balancing':'%s/Analysis/mb_hg19_GM12878_correlations.txt' % basedir,
    }
    mm9_data = load_data(mm9_methods)
    hg19_data = load_data(hg19_methods)
    width = 16.8
    spacer = 0.4
    overall_width = (width - spacer * 2) / 2.6
    c = canvas.canvas()
    mm9_ranges_img, mm9_ranges_height = plot_dataset_ranges(mm9_data, width, "MM9 ESC")
    mm9_ranges_img.text(0, mm9_ranges_height, 'a',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(mm9_ranges_img)
    hg19_ranges_img, hg19_ranges_height = plot_dataset_ranges(hg19_data, width, "HG19 GM12878")
    hg19_ranges_img.text(0, hg19_ranges_height, 'b',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(hg19_ranges_img, [trafo.translate(0, -hg19_ranges_height - spacer)])
    overall_height = mm9_ranges_height * 0.6
    mm9_overall_img = plot_overall(mm9_data, overall_width, overall_height, "MM9 ESC")
    mm9_overall_img.text(0, overall_height + 0.1, 'c',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(mm9_overall_img, [trafo.translate(0, -hg19_ranges_height - overall_height - spacer * 2)])
    hg19_overall_img = plot_overall(hg19_data, overall_width, overall_height, "HG19 GM12878")
    hg19_overall_img.text(0, overall_height + 0.1, 'd',
                        [text.halign.left, text.valign.top, text.size(-1)])
    c.insert(hg19_overall_img, [trafo.translate(overall_width * 1.6 + spacer * 2,
                                -hg19_ranges_height - overall_height - spacer * 2)])
    c.insert(plot_key(overall_width * 0.6 + 0.4, overall_height),
             [trafo.translate(overall_width + spacer + 0.6, -hg19_ranges_height - overall_height - spacer * 2)])
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

def plot_overall(data, width, height, name):
    vo = 0.55
    ho = 0.7
    plot_width = width - ho
    plot_height = height - vo - 0.3
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
    c.insert(g, [trafo.translate(ho, vo)])
    for i, label in enumerate(["10K", "50K", "250K", "1M", "250K", "1M"]):
        c.text(ho + plot_width * (i + 0.5) / 6.0, vo - 0.05, "%s" % label,
               [text.halign.center, text.valign.top, text.size(-3)])
    c.text(ho + plot_width * 2.0 / 6.0, 0.05, "cis",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2, vo * 0.5, ho - 0.2 + plot_width * 4.0 / 6.0, vo * 0.5), [style.linewidth.THin])
    c.text(ho + plot_width * 5.0 / 6.0, 0.05, "trans",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    c.stroke(path.line(ho + 0.2 + plot_width * 4.0 / 6.0, vo * 0.5, ho - 0.2 + plot_width, vo * 0.5), [style.linewidth.THin])
    c.text(0, plot_height * 0.5 + vo, "Correlation",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    c.text(plot_width * 0.5 + ho, height, name,
           [text.halign.center, text.valign.top, text.size(-3)])
    return c

def plot_dataset_ranges(data, width, label):
    methods = data.keys()
    binsizes = numpy.unique(data[methods[0]]['binsize'])
    ho = 0.4
    ho2 = 0.4
    vo = 0.6
    spacer = 0.25
    plot_width = (width - ho * 2 - (binsizes.shape[0] - 1) * spacer) / binsizes.shape[0] - ho2
    plot_height = plot_width
    c = canvas.canvas()
    for i, binsize in enumerate(binsizes):
        img = plot_single_range(data, binsize, plot_width, plot_width)
        c.insert(img, [trafo.translate((plot_width + spacer) * i + ho2 * (i + 1) + ho, vo)])
    c.text(0, plot_height * 0.5 + vo, "Correlation",
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(90)])
    c.text(width, plot_height * 0.5 + vo, label,
           [text.halign.center, text.valign.top, text.size(-3), trafo.rotate(-90)])
    c.text((plot_width + ho2) * 2 + spacer * 1.5 + ho, 0, "Interaction Range",
           [text.halign.center, text.valign.bottom, text.size(-3)])
    return c, plot_height + vo + 0.3

def plot_single_range(data, binsize, width, height):
    plot_width = width
    plot_height = height
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
        binstring = "%iM" % (binsize / 1000000)
    elif binsize / 1000 > 0:
        binstring = "%iK" % (binsize / 1000)
    else:
        binstring = str(binsize)
    g.text(plot_width / 2, plot_height + 0.3, "%s binning" % (binstring),
           [text.halign.center, text.valign.top, text.size(-2)])
    c.insert(g)
    return c

def plot_key(width, height):
    c = canvas.canvas()
    step = height / float(len(method_colors))
    for i, meth in enumerate(['HiFive-Probability', 'HiFive-Express', 'HiFive-Binning',
                              'HiCNorm', 'HiCPipe', 'Matrix-Balancing']):
        c.fill(path.rect(0.2, height - step * (i + 0.5) - 0.1, 0.2, 0.2),
               [method_colors[meth]])
        c.text(0.5, height - step * (i + 0.5), meth, [text.halign.left, text.valign.middle, text.size(-2)])
    return c
if __name__ == "__main__":
    main()
