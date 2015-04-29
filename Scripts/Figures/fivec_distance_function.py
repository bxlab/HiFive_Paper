#!/usr/bin/env python

import sys
import os

from pyx import *
from PIL import Image
import numpy
import scipy.stats
from scipy.stats import linregress
from math import log, exp

import hifive


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )


def main():
    basedir = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    out_fname = sys.argv[1]
    fivec_fname = "%s/Data/FiveC/HiFive/Nora_ESC_male_E14_prob.fcp" % basedir
    fivec = hifive.FiveC(fivec_fname)
    num_bins = 200
    width = 16.92
    plot_width = (width) / 2.0 - 2
    c = canvas.canvas()
    hist1, hist_min, hist_max, hist_ranges1 = get_histograms(fivec, num_bins, 'raw')
    temp = get_histograms(fivec, num_bins, 'fragment')
    hist2 = temp[0]
    hist_ranges2 = temp[3]
    hist_min = min(hist_min, temp[1])
    hist_max = max(hist_max, temp[2])
    c.insert(plot_distances(fivec, plot_width, hist1, hist_min, hist_max, hist_ranges1))
    c.insert(plot_distances(fivec, plot_width, hist2, hist_min, hist_max, hist_ranges2),
            [trafo.translate(plot_width + 0.8, 0)])
    c.text(-0.6, plot_width * 0.5, "Log count",
           [text.halign.center, text.valign.bottom, text.size(-3), trafo.rotate(90)])
    c.text(plot_width * 0.5, -0.4, "Log distance",
           [text.halign.center, text.valign.top, text.size(-3)])
    c.text(plot_width * 1.5 + 0.8, -0.4, "Log distance",
           [text.halign.center, text.valign.top, text.size(-3)])
    c.text(plot_width * 0.5, plot_width + 0.1, "Raw",
           [text.halign.center, text.valign.bottom, text.size(-2)])
    c.text(plot_width * 1.5 + 0.8, plot_width + 0.1, "Corrected",
           [text.halign.center, text.valign.bottom, text.size(-2)])
    c.insert(hifive.plotting.plot_key(hist_min, hist_max, 0.4, plot_width * 2 + 0.8, "%.0f", orientation='top',
             min_color="000000", max_color="0000ff", mid_color=None, labelattrs=[text.size(-3)]),
             [trafo.translate(0, plot_width + 0.4)])
    c.text(plot_width * 1.0 + 0.4, plot_width + 1.3, "Interactions per bin",
           [text.halign.center, text.valign.bottom, text.size(-2)])
    c.writePDFfile(out_fname)


def get_histograms(fivec, num_bins, datatype):
    hm, xmap, ymap = fivec.cis_heatmap(0, datatype=datatype, arraytype='compact', skipfiltered=True, returnmapping=True)
    where = numpy.where(hm[:, :, 0] > 0)
    Y = hm[where[0], where[1], 0] / hm[where[0], where[1], 1]
    X = numpy.abs(fivec.frags['fragments']['mid'][xmap[:, 2]].reshape(-1, 1) -
                  fivec.frags['fragments']['mid'][ymap[:, 2]].reshape(1, -1))[where]    
    xmin = numpy.amin(X)
    xmax = numpy.amax(X)
    ymin = numpy.amin(Y)
    ymax = numpy.amax(Y)
    hist_ranges = numpy.log(numpy.array([[xmin, xmax], [ymin, ymax]], dtype=numpy.float32))
    log_X = numpy.log(X)
    log_Y = numpy.log(Y)
    hist = numpy.histogram2d(log_X, log_Y, bins=num_bins, range=hist_ranges)[0]
    where = numpy.where(hist > 0)
    hist[where] = numpy.log(hist[where])
    hist_min = numpy.amin(hist[where])
    hist_max = numpy.amax(hist[where])
    valid = numpy.zeros(hist.shape, dtype=numpy.float32)
    valid[where] = 1.0
    hist = numpy.dstack((hist, valid))
    return hist, hist_min, hist_max, numpy.exp(hist_ranges)

def plot_distances(fivec, width, hist, hist_min, hist_max, hist_ranges):
    gamma = fivec.gamma
    mu = fivec.region_means[0]
    hist_bm = numpy.zeros((hist.shape[0], hist.shape[1]), dtype=numpy.uint32)
    hist_bm.fill(int('ffffffff', 16))
    hist1 = numpy.zeros(hist_bm.shape, dtype=numpy.float32)
    where = numpy.where(hist[:, :, 1])
    hist1[where] = hist[where[0], where[1], 0] / hist_max
    hist_bm[-1 - where[1], where[0]] = 256**3 * 255 + numpy.round(255 * hist1[where]).astype(numpy.int32) * 256 ** 2
    hist_img = Image.frombuffer('RGBA', (hist.shape[0], hist.shape[1]), hist_bm, 'raw', 'RGBA', 0, 1)
    c = canvas.canvas()
    c.insert(bitmap.bitmap(0, 0, hist_img, width=width))
    xmin = hist_ranges[0, 0] / 1000
    xmax = hist_ranges[0, 1] / 1000
    ymin = hist_ranges[1, 0]
    ymax = hist_ranges[1, 1]
    g = graph.graphxy( width=width, height=width, 
    x=graph.axis.log(min=xmin, max=xmax, title='', painter=painter), 
    y=graph.axis.log(min=ymin, max=ymax, title='', painter=painter), 
    x2=graph.axis.lin(min=0, max=1, parter=None),
    y2=graph.axis.lin(min=0, max=1, parter=None) )
    g.plot(graph.data.function("y(x) = exp(-%f * log(x*1000) + %f)" % (gamma, mu), min=xmin, max=xmax),
           [graph.style.line([color.cmyk.Red, style.linewidth.THick])])
    c.insert(g)
    return c


def CreateKey(min_val, max_val, height, width, log_display=True, orientation='left', num_ticks=9):
    height = float(height)
    width = float(width)
    min_val = float(min_val)
    max_val = float(max_val)
    c = canvas.canvas()
    if orientation in ['left', 'right']:
        img = numpy.zeros( (int(round(256.0*width/height)), 256), dtype=numpy.uint32 )
        img.shape = (img.shape[1],img.shape[0])
        img[:, :] = 256**3 * 255 + (255 - numpy.arange(256).reshape(-1, 1)) * 256**2
        pilImage = Image.frombuffer( 'RGBA',(img.shape[1],img.shape[0]),img,'raw','RGBA',0,1)
        c.insert( bitmap.bitmap( 0, 0, pilImage, height=height ) )
    else:
        img = numpy.zeros( (256, int(round(256.0*height/width))), dtype=numpy.uint32 )
        img.shape = (img.shape[1],img.shape[0])
        img[:, :] = 256**3 * 255 + (255 - (numpy.arange(256)[::-1]).reshape(1, -1)) * 256**2
        pilImage = Image.frombuffer( 'RGBA',(img.shape[1],img.shape[0]),img,'raw','RGBA',0,1)
        c.insert( bitmap.bitmap( 0, 0, pilImage, width=width ) )
    c.stroke( path.rect( 0, 0, width, height ), [style.linewidth.THin] )
    if orientation in ['left', 'right']:
        tick_step = height/(num_ticks-1)
    else:
        tick_step = width/(num_ticks-1)
    if log_display == True:
        labels = numpy.exp( min_val + numpy.arange(num_ticks) * (max_val - min_val) / (num_ticks-1) )
    else:
        labels = min_val + numpy.arange(num_ticks) * (max_val - min_val) / (num_ticks-1)
    if orientation == 'left':
        for i in range(num_ticks):
            c.stroke( path.line( -width * 0.4, tick_step * i, 0.0, tick_step * i ), [style.linewidth.THin] )
            c.text(-width * 0.5, tick_step*i, r"%i" % (numpy.round(labels[i])),
                   [text.halign.right, text.valign.middle, text.size(-3)] )
    elif orientation == 'right':
        for i in range(num_ticks):
            c.stroke( path.line( width * 1.4, tick_step * i, width, tick_step * i ), [style.linewidth.THin] )
            c.text(width * 1.5, tick_step * i, r"%i" % (numpy.round(labels[i])),
                   [text.halign.left, text.valign.middle, text.size(-3)] )
    elif orientation == 'top':
        for i in range(num_ticks):
            c.stroke( path.line(tick_step * i, height, tick_step * i, height * 1.4 ), [style.linewidth.THin] )
            c.text(tick_step*i, height * 1.5, r"%i" % (numpy.round(labels[i])),
                   [text.halign.center, text.valign.bottom, text.size(-3)] )
    else:
        for i in range(num_ticks):
            c.stroke( path.line(tick_step * i, 0, tick_step * i, -height * 0.4 ), [style.linewidth.THin] )
            c.text(tick_step*i, -height * 0.5, r"%i" % (numpy.round(labels[i])),
                   [text.halign.center, text.valign.top, text.size(-3)] )
    return c


def CreateFigureB(hic_data, width, height):
    results = {}
    plot_width = (width - 1.6) / 3.0
    input = open(hic_data, 'r')
    for line in input:
        temp = line.strip('\n').split('\t')
        key = temp[0]
        temp = temp[1:]
        for i in range(len(temp)):
            temp[i] = float(temp[i])
        results[key] = numpy.array(temp, dtype=numpy.float32)
    input.close()
    c = canvas.canvas()
    xmin = 1.0
    xmax = 1.0
    for i, dataset in enumerate(['uncorrected', 'hifive', 'express']):
        X = results['%s_bins' % dataset] / 1000.0
        xmin = min(xmin, numpy.amin(X))
        xmax = max(xmax, numpy.amax(X))
    for i, dataset in enumerate(['uncorrected', 'hifive', 'express']):
        X = results['%s_bins' % dataset] / 1000.0   
        ymin = numpy.inf
        ymax = -numpy.inf
        for key in results:
            if key.count('bins') > 0 or key.count(dataset) == 0:
                continue
            where = numpy.where(results[key] > 0.0)
            ymin = min(numpy.amin(results[key][where]), ymin)
            ymax = max(numpy.amax(results[key][where]), ymax)
        if i == 0:
            yaxis = graph.axis.log(min=ymin, max=ymax, title='', painter=painter)
        else:
            yaxis = graph.axis.log(min=ymin, max=ymax, title='', parter=None)
        g = graph.graphxy( width=plot_width, height=plot_width, 
        x=graph.axis.log(min=xmin, max=xmax, title='', painter=painter), 
        y=yaxis, 
        x2=graph.axis.lin(min=0, max=1, parter=None),
        y2=graph.axis.lin(min=0, max=1, parter=None) )
        for key in results.keys():
            if key.count('bins') > 0 or key.count(dataset) == 0:
                continue
            temp = key.split('_')
            if temp[0] == 'raw' or temp[0] == 'smoothed':
                continue
            where = numpy.where(results[key] > 0.0)
            g.plot(graph.data.values(x=X[where], y=results[key][where]),
                   [graph.style.line([style.linestyle.dotted])] )
        where = numpy.where(results['raw_%s' % dataset] > 0.0)
        g.plot(graph.data.values(x=X[where], y=results['raw_%s' % dataset][where]),
               [graph.style.line([])] )
        where = numpy.where(results['smoothed_%s' % dataset] > 0.0)
        g.plot(graph.data.values(x=X[where], y=results['smoothed_%s' % dataset][where]),
               [graph.style.line([color.cmyk.Red])] )
        c.insert(g, [trafo.translate(1.2 + plot_width * i, 1.0)])
    c.text(0, 1.0 + plot_width * 0.5, r"Mean Interaction Count",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(90)])
    c.text(1.2 + plot_width * 1.5, 0, r"Inter-fragment Distance (Kb)",
           [text.halign.center, text.valign.bottom, text.size(-2)])
    c.text(0, height, r"b", [text.halign.left, text.valign.top, text.size(0)])
    c.text(width, 1.0 + plot_width * 0.5, r"HiC",
           [text.halign.center, text.valign.top, text.size(-2), trafo.rotate(-90)])
    c.text(1.2 + plot_width * 0.5, 1.2 + plot_width, r"Raw",
           [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(1.2 + plot_width * 1.5, 1.2 + plot_width, r"HiFive corrected",
           [text.halign.center, text.valign.middle, text.size(-3)])
    c.text(1.2 + plot_width * 2.5, 1.2 + plot_width, r"HiFive-Express corrected",
           [text.halign.center, text.valign.middle, text.size(-3)])
    return c


if __name__ == "__main__":
    main()
